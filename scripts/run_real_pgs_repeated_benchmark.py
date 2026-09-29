from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
import statistics
import sys
from time import perf_counter


ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from genome_he_pilot import prs_scores  # noqa: E402
from genome_he_pilot.memory import process_memory, stop_process_memory_measurement  # noqa: E402
from genome_he_pilot.pgs import read_pgs_scoring_file  # noqa: E402
from genome_he_pilot.pgs_vcf import load_pgs_matched_vcf_dataset  # noqa: E402
from genome_he_pilot.tenseal_backend import TenSEALUnavailable, encrypted_prs_benchmark  # noqa: E402
from genome_he_pilot.vcf_cache import load_pgs_matched_cached_dataset  # noqa: E402


DEFAULT_VCF = (
    ROOT
    / "data"
    / "1000genomes"
    / "phase3"
    / "ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"
)
DEFAULT_PGS = ROOT / "data" / "pgs_catalog" / "PGS004941" / "PGS004941.txt.gz"
DEFAULT_CACHE = ROOT / "data" / "cache" / "1000g_chr22_samples32.sqlite"
DEFAULT_RAW_OUTPUT = ROOT / "results" / "real_pgs_chr22_pgs004941_repeated_raw.csv"
DEFAULT_SUMMARY_OUTPUT = ROOT / "results" / "real_pgs_chr22_pgs004941_repeated_summary.csv"
DEFAULT_REPORT = ROOT / "submission" / "real_pgs_repeated_benchmark_audit.md"

SUMMARY_METRICS = [
    "match_vcf_seconds",
    "plain_prs_seconds",
    "context_seconds",
    "encryption_seconds",
    "evaluation_seconds",
    "decryption_seconds",
    "serialization_seconds",
    "total_ckks_compute_seconds",
    "max_abs_error",
    "eval_overhead_ratio",
    "total_overhead_ratio",
    "ckks_samples_per_second",
    "ckks_snp_products_per_second",
    "plain_samples_per_second",
    "plain_snp_products_per_second",
    "input_ciphertext_mean_bytes",
    "result_ciphertext_mean_bytes",
    "rss_delta_bytes",
    "rss_delta_per_sample_bytes",
]


@dataclass(frozen=True)
class Check:
    name: str
    passed: bool
    detail: str


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--pgs", type=Path, default=DEFAULT_PGS)
    parser.add_argument("--vcf", type=Path, default=DEFAULT_VCF)
    parser.add_argument("--chrom", default="22")
    parser.add_argument("--samples", type=int, default=32)
    parser.add_argument("--cache", type=Path, default=DEFAULT_CACHE)
    parser.add_argument("--no-cache", action="store_true")
    parser.add_argument("--repeats", type=int, default=5)
    parser.add_argument("--max-snps", type=int, nargs="+", default=[512, 1024, 2048, 4096])
    parser.add_argument("--raw-output", type=Path, default=DEFAULT_RAW_OUTPUT)
    parser.add_argument("--summary-output", type=Path, default=DEFAULT_SUMMARY_OUTPUT)
    parser.add_argument("--report", type=Path, default=DEFAULT_REPORT)
    return parser.parse_args()


def format_float(value: float) -> str:
    return f"{value:.12g}"


def mean(values: list[float]) -> float:
    return sum(values) / len(values)


def sample_sd(values: list[float]) -> float:
    if len(values) < 2:
        return 0.0
    return statistics.stdev(values)


def float_or_none(value: str) -> float | None:
    if value == "":
        return None
    try:
        return float(value)
    except ValueError:
        return None


def run_once(args: argparse.Namespace, scoring, max_snps: int, repeat: int) -> dict[str, str]:
    memory_before = process_memory()

    start = perf_counter()
    if args.cache.exists() and not args.no_cache:
        match = load_pgs_matched_cached_dataset(
            cache_path=args.cache,
            scoring=scoring,
            chrom=args.chrom,
            max_snps=max_snps,
        )
        match_mode = "sqlite_cache"
    else:
        match = load_pgs_matched_vcf_dataset(
            vcf_path=args.vcf,
            scoring=scoring,
            chrom=args.chrom,
            n_samples=args.samples,
            max_snps=max_snps,
        )
        match_mode = "sequential_vcf_scan"
    match_seconds = perf_counter() - start
    dataset = match.dataset

    start = perf_counter()
    plain_scores = prs_scores(dataset.genotypes, dataset.weights)
    plain_seconds = perf_counter() - start

    status = "tenseal_ok"
    context_seconds = encryption_seconds = evaluation_seconds = 0.0
    decryption_seconds = serialization_seconds = 0.0
    max_abs_error = 0.0
    public_context_bytes = input_ciphertext_bytes = result_ciphertext_bytes = 0
    poly_modulus_degree = ""
    scale = ""
    coeff_mod_bit_sizes = ""

    try:
        encrypted = encrypted_prs_benchmark(dataset.genotypes, dataset.weights)
        max_abs_error = max(abs(plain - secure) for plain, secure in zip(plain_scores, encrypted.scores))
        context_seconds = encrypted.context_seconds
        encryption_seconds = encrypted.encryption_seconds
        evaluation_seconds = encrypted.evaluation_seconds
        decryption_seconds = encrypted.decryption_seconds
        serialization_seconds = encrypted.serialization_seconds
        public_context_bytes = encrypted.public_context_bytes
        input_ciphertext_bytes = encrypted.input_ciphertext_bytes
        result_ciphertext_bytes = encrypted.result_ciphertext_bytes
        poly_modulus_degree = str(encrypted.poly_modulus_degree)
        scale = str(encrypted.scale)
        coeff_mod_bit_sizes = "-".join(str(value) for value in encrypted.coeff_mod_bit_sizes)
    except TenSEALUnavailable as exc:
        status = f"tenseal_unavailable: {exc}"

    memory = stop_process_memory_measurement(memory_before)
    total_ckks_seconds = encryption_seconds + evaluation_seconds + decryption_seconds
    eval_overhead = evaluation_seconds / plain_seconds if plain_seconds > 0 else 0.0
    total_overhead = total_ckks_seconds / plain_seconds if plain_seconds > 0 else 0.0
    ckks_samples_per_second = dataset.n_samples / evaluation_seconds if evaluation_seconds > 0 else 0.0
    ckks_snp_products_per_second = (
        dataset.n_samples * dataset.n_snps / evaluation_seconds if evaluation_seconds > 0 else 0.0
    )
    plain_samples_per_second = dataset.n_samples / plain_seconds if plain_seconds > 0 else 0.0
    plain_snp_products_per_second = (
        dataset.n_samples * dataset.n_snps / plain_seconds if plain_seconds > 0 else 0.0
    )
    input_ciphertext_mean_bytes = input_ciphertext_bytes / dataset.n_samples if dataset.n_samples else 0.0
    result_ciphertext_mean_bytes = result_ciphertext_bytes / dataset.n_samples if dataset.n_samples else 0.0
    rss_delta_per_sample = memory.rss_delta_bytes / dataset.n_samples if dataset.n_samples else 0.0

    return {
        "pgs_id": scoring.metadata.get("pgs_id", ""),
        "trait_reported": scoring.metadata.get("trait_reported", ""),
        "genome_build": scoring.metadata.get("genome_build", ""),
        "chrom": args.chrom,
        "repeat": str(repeat),
        "samples": str(dataset.n_samples),
        "matched_snps": str(dataset.n_snps),
        "max_snps": str(max_snps),
        "requested_weights_on_chrom": str(match.requested_weights_on_chrom),
        "skipped_allele_mismatch": str(match.skipped_allele_mismatch),
        "backend": "tenseal",
        "status": status,
        "match_mode": match_mode,
        "match_vcf_seconds": format_float(match_seconds),
        "plain_prs_seconds": format_float(plain_seconds),
        "context_seconds": format_float(context_seconds),
        "encryption_seconds": format_float(encryption_seconds),
        "evaluation_seconds": format_float(evaluation_seconds),
        "decryption_seconds": format_float(decryption_seconds),
        "serialization_seconds": format_float(serialization_seconds),
        "total_ckks_compute_seconds": format_float(total_ckks_seconds),
        "max_abs_error": format_float(max_abs_error),
        "eval_overhead_ratio": format_float(eval_overhead),
        "total_overhead_ratio": format_float(total_overhead),
        "ckks_samples_per_second": format_float(ckks_samples_per_second),
        "ckks_snp_products_per_second": format_float(ckks_snp_products_per_second),
        "plain_samples_per_second": format_float(plain_samples_per_second),
        "plain_snp_products_per_second": format_float(plain_snp_products_per_second),
        "input_ciphertext_bytes": str(input_ciphertext_bytes),
        "input_ciphertext_mean_bytes": format_float(input_ciphertext_mean_bytes),
        "result_ciphertext_bytes": str(result_ciphertext_bytes),
        "result_ciphertext_mean_bytes": format_float(result_ciphertext_mean_bytes),
        "public_context_bytes": str(public_context_bytes),
        "poly_modulus_degree": poly_modulus_degree,
        "scale": scale,
        "coeff_mod_bit_sizes": coeff_mod_bit_sizes,
        "rss_before_bytes": str(memory.rss_before_bytes),
        "rss_after_bytes": str(memory.rss_after_bytes),
        "rss_delta_bytes": str(memory.rss_delta_bytes),
        "rss_delta_per_sample_bytes": format_float(rss_delta_per_sample),
        "peak_rss_bytes": str(memory.peak_rss_bytes),
    }


def write_csv(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        raise ValueError("cannot write an empty benchmark CSV")
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def summarize_rows(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    groups: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        groups[row["max_snps"]].append(row)

    summary_rows: list[dict[str, str]] = []
    for max_snps in sorted(groups, key=lambda value: int(value)):
        group = groups[max_snps]
        first = group[0]
        summary: dict[str, str] = {
            "pgs_id": first["pgs_id"],
            "chrom": first["chrom"],
            "samples": first["samples"],
            "max_snps": first["max_snps"],
            "matched_snps": first["matched_snps"],
            "repeats": str(len(group)),
            "ok_repeats": str(sum(1 for row in group if row["status"] == "tenseal_ok")),
            "match_mode": first["match_mode"],
        }
        for metric in SUMMARY_METRICS:
            values = [float(value) for row in group if (value := float_or_none(row.get(metric, ""))) is not None]
            if values:
                summary[f"{metric}_mean"] = format_float(mean(values))
                summary[f"{metric}_sd"] = format_float(sample_sd(values))
                summary[f"{metric}_min"] = format_float(min(values))
                summary[f"{metric}_max"] = format_float(max(values))
            else:
                summary[f"{metric}_mean"] = ""
                summary[f"{metric}_sd"] = ""
                summary[f"{metric}_min"] = ""
                summary[f"{metric}_max"] = ""
        summary_rows.append(summary)
    return summary_rows


def collect_checks(args: argparse.Namespace, rows: list[dict[str, str]], summary_rows: list[dict[str, str]]) -> list[Check]:
    expected_raw_rows = args.repeats * len(args.max_snps)
    observed_snps = sorted({int(row["max_snps"]) for row in rows})
    checks = [
        Check("Raw repeated benchmark rows exist", len(rows) == expected_raw_rows, f"rows={len(rows)} expected={expected_raw_rows}"),
        Check("Summary rows cover all SNP sizes", observed_snps == sorted(args.max_snps), f"snps={observed_snps}"),
        Check("All repeats completed with TenSEAL", all(row["status"] == "tenseal_ok" for row in rows), "all ok"),
        Check("Plaintext baseline time recorded", all(float(row["plain_prs_seconds"]) > 0 for row in rows), "present"),
        Check("CKKS evaluation time recorded", all(float(row["evaluation_seconds"]) > 0 for row in rows), "present"),
        Check("Overhead ratios recorded", all(float(row["eval_overhead_ratio"]) > 0 for row in rows), "present"),
        Check("Throughput metrics recorded", all(float(row["ckks_snp_products_per_second"]) > 0 for row in rows), "present"),
        Check("Ciphertext size per sample recorded", all(float(row["input_ciphertext_mean_bytes"]) > 0 for row in rows), "present"),
        Check("Memory per sample recorded", all(row["rss_delta_per_sample_bytes"] != "" for row in rows), "present"),
        Check("CKKS max error remains below 1e-5", max(float(row["max_abs_error"]) for row in rows) < 1e-5, f"max={max(float(row['max_abs_error']) for row in rows):.3g}"),
        Check("Summary includes standard deviation", all(row["evaluation_seconds_sd"] != "" for row in summary_rows), "present"),
    ]
    return checks


def passfail(value: bool) -> str:
    return "PASS" if value else "FAIL"


def write_report(path: Path, checks: list[Check], summary_rows: list[dict[str, str]]) -> None:
    failures = [check for check in checks if not check.passed]
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = [
        "# Real PGS Repeated Benchmark Audit",
        "",
        "## Verdict",
        "",
        "PASS: repeated PGS004941 chr22 benchmark completed and summary metrics are available."
        if not failures
        else "FAIL: repeated benchmark has unresolved issues.",
        "",
        "## Summary Table",
        "",
        "| SNPs | Repeats | Plain sec mean | CKKS eval sec mean | Eval overhead mean | CKKS SNP products/sec mean | Max error mean |",
        "| ---: | ---: | ---: | ---: | ---: | ---: | ---: |",
    ]
    for row in summary_rows:
        lines.append(
            "| "
            + " | ".join(
                [
                    row["matched_snps"],
                    row["repeats"],
                    row["plain_prs_seconds_mean"],
                    row["evaluation_seconds_mean"],
                    row["eval_overhead_ratio_mean"],
                    row["ckks_snp_products_per_second_mean"],
                    row["max_abs_error_mean"],
                ]
            )
            + " |"
        )
    lines.extend(
        [
            "",
            "## Checks",
            "",
            "| Check | Status | Detail |",
            "| --- | --- | --- |",
        ]
    )
    lines.extend(f"| {check.name} | {passfail(check.passed)} | `{check.detail}` |" for check in checks)
    if failures:
        lines.extend(["", "## Failures", ""])
        lines.extend(f"- {check.name}: `{check.detail}`" for check in failures)
    lines.append("")
    path.write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    args = parse_args()
    if args.repeats < 2:
        raise ValueError("--repeats must be at least 2 to report standard deviation")

    scoring = read_pgs_scoring_file(args.pgs, chrom=args.chrom)
    rows: list[dict[str, str]] = []
    for max_snps in args.max_snps:
        for repeat in range(1, args.repeats + 1):
            print(f"Repeated PGS benchmark: snps={max_snps} repeat={repeat}/{args.repeats}", flush=True)
            row = run_once(args, scoring, max_snps, repeat)
            rows.append(row)
            print(
                "  ",
                f"plain={row['plain_prs_seconds']}",
                f"eval={row['evaluation_seconds']}",
                f"overhead={row['eval_overhead_ratio']}",
                f"error={row['max_abs_error']}",
                flush=True,
            )

    write_csv(args.raw_output, rows)
    summary_rows = summarize_rows(rows)
    write_csv(args.summary_output, summary_rows)
    checks = collect_checks(args, rows, summary_rows)
    write_report(args.report, checks, summary_rows)

    print(f"Wrote {args.raw_output}")
    print(f"Wrote {args.summary_output}")
    print(f"Wrote {args.report}")
    failures = [check for check in checks if not check.passed]
    if failures:
        print("Repeated benchmark audit failures:")
        for check in failures:
            print(f"- {check.name}: {check.detail}")
        return 1
    print("Repeated benchmark audit passed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
