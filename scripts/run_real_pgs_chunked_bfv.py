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

from genome_he_pilot.chunked_fixed_point import (  # noqa: E402
    DEFAULT_FIXED_POINT_WINDOW_SNPS,
    chunked_bfv_fixed_point_prs_benchmark,
)
from genome_he_pilot.bfv_backend import DEFAULT_BFV_FIXED_POINT_SCALE  # noqa: E402
from genome_he_pilot.pgs import read_pgs_scoring_file  # noqa: E402
from genome_he_pilot.pgs_vcf import load_pgs_matched_vcf_dataset  # noqa: E402
from genome_he_pilot.plaintext import prs_scores  # noqa: E402
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
DEFAULT_RAW_OUTPUT = ROOT / "results" / "real_pgs_chr22_pgs004941_chunked_bfv_raw.csv"
DEFAULT_SUMMARY_OUTPUT = ROOT / "results" / "real_pgs_chr22_pgs004941_chunked_bfv_summary.csv"
DEFAULT_REPORT = ROOT / "submission" / "real_pgs_chunked_bfv_scaling_audit.md"

SUMMARY_METRICS = [
    "match_vcf_seconds",
    "plain_prs_seconds",
    "context_seconds",
    "encryption_seconds",
    "evaluation_seconds",
    "decryption_seconds",
    "serialization_seconds",
    "total_he_compute_seconds",
    "max_abs_error",
    "eval_overhead_ratio",
    "total_overhead_ratio",
    "he_samples_per_second",
    "he_snp_products_per_second",
    "input_ciphertext_total_bytes",
    "input_ciphertext_mean_bytes_per_sample",
    "input_ciphertext_mean_bytes_per_chunk",
    "result_ciphertext_mean_bytes_per_sample",
    "fixed_point_scale",
    "plain_modulus",
    "integer_dot_abs_bound",
    "safe_modulus_margin",
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
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--window-snps", type=int, default=DEFAULT_FIXED_POINT_WINDOW_SNPS)
    parser.add_argument("--max-snps", type=int, nargs="+", default=[8192, 16384, 32768, 51227])
    parser.add_argument("--fixed-point-scale", type=int, default=DEFAULT_BFV_FIXED_POINT_SCALE)
    parser.add_argument("--raw-output", type=Path, default=DEFAULT_RAW_OUTPUT)
    parser.add_argument("--summary-output", type=Path, default=DEFAULT_SUMMARY_OUTPUT)
    parser.add_argument("--report", type=Path, default=DEFAULT_REPORT)
    return parser.parse_args()


def format_float(value: float) -> str:
    return f"{value:.12g}"


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


def load_dataset(args: argparse.Namespace, scoring, max_snps: int):
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
    return match, match_mode, perf_counter() - start


def run_once(args: argparse.Namespace, scoring, max_snps: int, repeat: int) -> dict[str, str]:
    match, match_mode, match_seconds = load_dataset(args, scoring, max_snps)
    dataset = match.dataset

    start = perf_counter()
    plain_scores = prs_scores(dataset.genotypes, dataset.weights)
    plain_seconds = perf_counter() - start

    result = chunked_bfv_fixed_point_prs_benchmark(
        dataset.genotypes,
        dataset.weights,
        window_snps=args.window_snps,
        fixed_point_scale=args.fixed_point_scale,
    )
    max_error = max(abs(actual - expected) for actual, expected in zip(result.scores, plain_scores))
    total_he_seconds = result.encryption_seconds + result.evaluation_seconds + result.decryption_seconds
    samples = dataset.n_samples
    snps = dataset.n_snps

    return {
        "pgs_id": scoring.metadata.get("pgs_id", ""),
        "trait_reported": scoring.metadata.get("trait_reported", ""),
        "genome_build": scoring.metadata.get("genome_build", ""),
        "chrom": args.chrom,
        "repeat": str(repeat),
        "samples": str(samples),
        "matched_snps": str(snps),
        "max_snps": str(max_snps),
        "window_snps": str(result.window_snps),
        "windows_per_sample": str(result.windows_per_sample),
        "input_ciphertexts": str(result.input_ciphertexts),
        "result_ciphertexts": str(result.result_ciphertexts),
        "requested_weights_on_chrom": str(match.requested_weights_on_chrom),
        "skipped_allele_mismatch": str(match.skipped_allele_mismatch),
        "match_mode": match_mode,
        "match_vcf_seconds": format_float(match_seconds),
        "scheme": "BFV",
        "backend": "TenSEAL BFVVector chunked",
        "status": "ok",
        "plain_prs_seconds": format_float(plain_seconds),
        "context_seconds": format_float(result.context_seconds),
        "encryption_seconds": format_float(result.encryption_seconds),
        "evaluation_seconds": format_float(result.evaluation_seconds),
        "decryption_seconds": format_float(result.decryption_seconds),
        "serialization_seconds": format_float(result.serialization_seconds),
        "total_he_compute_seconds": format_float(total_he_seconds),
        "max_abs_error": format_float(max_error),
        "eval_overhead_ratio": format_float(result.evaluation_seconds / plain_seconds if plain_seconds > 0 else 0.0),
        "total_overhead_ratio": format_float(total_he_seconds / plain_seconds if plain_seconds > 0 else 0.0),
        "he_samples_per_second": format_float(samples / result.evaluation_seconds if result.evaluation_seconds > 0 else 0.0),
        "he_snp_products_per_second": format_float(samples * snps / result.evaluation_seconds if result.evaluation_seconds > 0 else 0.0),
        "input_ciphertext_total_bytes": str(result.input_ciphertext_bytes),
        "input_ciphertext_mean_bytes_per_sample": format_float(result.input_ciphertext_bytes / samples if samples else 0.0),
        "input_ciphertext_mean_bytes_per_chunk": format_float(result.input_ciphertext_bytes / result.input_ciphertexts if result.input_ciphertexts else 0.0),
        "result_ciphertext_total_bytes": str(result.result_ciphertext_bytes),
        "result_ciphertext_mean_bytes_per_sample": format_float(result.result_ciphertext_bytes / samples if samples else 0.0),
        "public_context_bytes": str(result.public_context_bytes),
        "poly_modulus_degree": str(result.poly_modulus_degree),
        "plain_modulus": str(result.plain_modulus),
        "fixed_point_scale": str(result.fixed_point_scale),
        "integer_dot_abs_bound": str(result.integer_dot_abs_bound),
        "safe_modulus_margin": str(result.safe_modulus_margin),
        "note": "Chunked TenSEAL BFV fixed-point PRS dot product",
    }


def write_csv(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        raise ValueError("cannot write an empty CSV")
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def summarize_rows(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    groups: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        groups[row["max_snps"]].append(row)
    summary_rows: list[dict[str, str]] = []
    for max_snps, group in sorted(groups.items(), key=lambda item: int(item[0])):
        first = group[0]
        summary = {
            "pgs_id": first["pgs_id"],
            "chrom": first["chrom"],
            "samples": first["samples"],
            "max_snps": max_snps,
            "matched_snps": first["matched_snps"],
            "window_snps": first["window_snps"],
            "windows_per_sample": first["windows_per_sample"],
            "input_ciphertexts": first["input_ciphertexts"],
            "result_ciphertexts": first["result_ciphertexts"],
            "scheme": "BFV",
            "backend": first["backend"],
            "repeats": str(len(group)),
            "ok_repeats": str(sum(1 for row in group if row["status"] == "ok")),
            "status": "ok" if all(row["status"] == "ok" for row in group) else "incomplete",
            "match_mode": first["match_mode"],
            "note": first["note"],
        }
        for metric in SUMMARY_METRICS:
            present = [value for row in group if (value := float_or_none(row.get(metric, ""))) is not None]
            if present:
                summary[f"{metric}_mean"] = format_float(sum(present) / len(present))
                summary[f"{metric}_sd"] = format_float(sample_sd(present))
                summary[f"{metric}_min"] = format_float(min(present))
                summary[f"{metric}_max"] = format_float(max(present))
            else:
                summary[f"{metric}_mean"] = ""
                summary[f"{metric}_sd"] = ""
                summary[f"{metric}_min"] = ""
                summary[f"{metric}_max"] = ""
        summary_rows.append(summary)
    return summary_rows


def passfail(value: bool) -> str:
    return "PASS" if value else "FAIL"


def write_report(path: Path, rows: list[dict[str, str]], summary_rows: list[dict[str, str]], expected_rows: int) -> list[Check]:
    checks = [
        Check("BFV raw rows completed", len(rows) == expected_rows, f"rows={len(rows)} expected={expected_rows}"),
        Check("All BFV rows are ok", all(row["status"] == "ok" for row in rows), "all ok"),
        Check(
            "BFV fixed-point error remains below 1e-3",
            bool(rows) and max(float(row["max_abs_error"]) for row in rows) < 1e-3,
            f"max={max(float(row['max_abs_error']) for row in rows):.3g}" if rows else "missing",
        ),
        Check("BFV modulus margins are positive", all(float(row["safe_modulus_margin"]) > 0 for row in rows), "positive"),
    ]
    failures = [check for check in checks if not check.passed]
    lines = [
        "# Chunked BFV Real PGS Audit",
        "",
        "## Verdict",
        "",
        "PASS: Chunked BFV fixed-point PRS benchmark completed."
        if not failures
        else "FAIL: chunked BFV benchmark has unresolved issues.",
        "",
        "## Summary Table",
        "",
        "| SNPs | Scheme | Windows/sample | Eval sec mean | Max abs err mean | Input MB/sample mean |",
        "| ---: | --- | ---: | ---: | ---: | ---: |",
    ]
    for row in summary_rows:
        input_mb = ""
        if row["input_ciphertext_mean_bytes_per_sample_mean"]:
            input_mb = format_float(float(row["input_ciphertext_mean_bytes_per_sample_mean"]) / (1024 * 1024))
        lines.append(
            f"| {row['matched_snps']} | {row['scheme']} | {row['windows_per_sample']} | "
            f"{row['evaluation_seconds_mean']} | {row['max_abs_error_mean']} | {input_mb} |"
        )
    lines.extend(["", "## Checks", "", "| Check | Status | Detail |", "| --- | --- | --- |"])
    lines.extend(f"| {check.name} | {passfail(check.passed)} | `{check.detail}` |" for check in checks)
    if failures:
        lines.extend(["", "## Failures", ""])
        lines.extend(f"- {check.name}: `{check.detail}`" for check in failures)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return checks


def main() -> int:
    args = parse_args()
    scoring = read_pgs_scoring_file(args.pgs, chrom=args.chrom)
    rows: list[dict[str, str]] = []
    for max_snps in args.max_snps:
        for repeat in range(1, args.repeats + 1):
            print(f"Chunked BFV PGS benchmark: snps={max_snps} repeat={repeat}/{args.repeats}", flush=True)
            row = run_once(args, scoring, max_snps, repeat)
            rows.append(row)
            print(f"  eval={row['evaluation_seconds']} err={row['max_abs_error']}", flush=True)
    write_csv(args.raw_output, rows)
    summary_rows = summarize_rows(rows)
    write_csv(args.summary_output, summary_rows)
    checks = write_report(args.report, rows, summary_rows, args.repeats * len(args.max_snps))
    print(f"Wrote {args.raw_output}")
    print(f"Wrote {args.summary_output}")
    print(f"Wrote {args.report}")
    failures = [check for check in checks if not check.passed]
    if failures:
        print("Chunked BFV audit failures:")
        for check in failures:
            print(f"- {check.name}: {check.detail}")
        return 1
    print("Chunked BFV audit passed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
