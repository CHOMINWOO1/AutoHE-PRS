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

from genome_he_pilot.openfhe_bgv_backend import (  # noqa: E402
    DEFAULT_BGV_FIXED_POINT_SCALE,
    DEFAULT_BGV_PLAIN_MODULUS,
    DEFAULT_BGV_POLY_MODULUS_DEGREE,
    openfhe_bgv_fixed_point_prs_benchmark,
)
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
DEFAULT_RAW_OUTPUT = ROOT / "results" / "real_pgs_chr22_pgs004941_openfhe_bgv_raw.csv"
DEFAULT_SUMMARY_OUTPUT = ROOT / "results" / "real_pgs_chr22_pgs004941_openfhe_bgv_summary.csv"
DEFAULT_REPORT = ROOT / "submission" / "openfhe_bgv_real_pgs_audit.md"

SUMMARY_METRICS = [
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
    "input_ciphertext_mean_bytes",
    "result_ciphertext_mean_bytes",
    "public_context_bytes",
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
    parser.add_argument("--max-snps", type=int, nargs="+", default=[512, 1024, 2048, 4096])
    parser.add_argument("--poly-modulus-degree", type=int, default=DEFAULT_BGV_POLY_MODULUS_DEGREE)
    parser.add_argument("--plain-modulus", type=int, default=DEFAULT_BGV_PLAIN_MODULUS)
    parser.add_argument("--fixed-point-scale", type=int, default=DEFAULT_BGV_FIXED_POINT_SCALE)
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


def bgv_row(
    args: argparse.Namespace,
    scoring,
    match,
    match_mode: str,
    match_seconds: float,
    max_snps: int,
    repeat: int,
) -> dict[str, str]:
    dataset = match.dataset

    start = perf_counter()
    plain_scores = prs_scores(dataset.genotypes, dataset.weights)
    plain_seconds = perf_counter() - start

    result = openfhe_bgv_fixed_point_prs_benchmark(
        dataset.genotypes,
        dataset.weights,
        poly_modulus_degree=args.poly_modulus_degree,
        plain_modulus=args.plain_modulus,
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
        "requested_weights_on_chrom": str(match.requested_weights_on_chrom),
        "skipped_allele_mismatch": str(match.skipped_allele_mismatch),
        "match_mode": match_mode,
        "match_vcf_seconds": format_float(match_seconds),
        "scheme": "BGV",
        "backend": "OpenFHE BGVRNS",
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
        "input_ciphertext_bytes": str(result.input_ciphertext_bytes),
        "input_ciphertext_mean_bytes": format_float(result.input_ciphertext_bytes / samples if samples else 0.0),
        "result_ciphertext_bytes": str(result.result_ciphertext_bytes),
        "result_ciphertext_mean_bytes": format_float(result.result_ciphertext_bytes / samples if samples else 0.0),
        "public_context_bytes": str(result.public_context_bytes),
        "poly_modulus_degree": str(result.poly_modulus_degree),
        "plain_modulus": str(result.plain_modulus),
        "fixed_point_scale": str(result.fixed_point_scale),
        "integer_dot_abs_bound": str(result.integer_dot_abs_bound),
        "safe_modulus_margin": str(result.safe_modulus_margin),
        "note": "OpenFHE BGV fixed-point PRS dot product; error includes weight quantization",
    }


def summarize_rows(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    groups: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        groups[(row["max_snps"], row["scheme"])].append(row)

    summary_rows: list[dict[str, str]] = []
    for (max_snps, scheme), group_rows in sorted(groups.items(), key=lambda item: (int(item[0][0]), item[0][1])):
        first = group_rows[0]
        summary = {
            "pgs_id": first["pgs_id"],
            "chrom": first["chrom"],
            "samples": first["samples"],
            "max_snps": max_snps,
            "matched_snps": first["matched_snps"],
            "scheme": scheme,
            "backend": first["backend"],
            "repeats": str(len(group_rows)),
            "ok_repeats": str(sum(1 for row in group_rows if row["status"] == "ok")),
            "status": "ok" if all(row["status"] == "ok" for row in group_rows) else "incomplete",
            "note": first["note"],
        }
        for metric in SUMMARY_METRICS:
            values = [float_or_none(row.get(metric, "")) for row in group_rows]
            present = [value for value in values if value is not None]
            if not present:
                for suffix in ("mean", "sd", "min", "max"):
                    summary[f"{metric}_{suffix}"] = ""
                continue
            summary[f"{metric}_mean"] = format_float(sum(present) / len(present))
            summary[f"{metric}_sd"] = format_float(sample_sd(present))
            summary[f"{metric}_min"] = format_float(min(present))
            summary[f"{metric}_max"] = format_float(max(present))
        summary_rows.append(summary)
    return summary_rows


def write_csv(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def passfail(value: bool) -> str:
    return "PASS" if value else "FAIL"


def write_audit(path: Path, rows: list[dict[str, str]], summary_rows: list[dict[str, str]], expected_rows: int) -> list[Check]:
    checks = [
        Check("BGV raw rows completed", len(rows) == expected_rows, f"rows={len(rows)} expected={expected_rows}"),
        Check("All BGV rows are ok", all(row["status"] == "ok" for row in rows), "all ok"),
        Check(
            "BGV fixed-point error remains below 1e-3",
            bool(rows) and max(float(row["max_abs_error"]) for row in rows) < 1e-3,
            f"max={max(float(row['max_abs_error']) for row in rows):.3g}" if rows else "missing",
        ),
        Check(
            "BGV modulus margins are positive",
            all(float(row["safe_modulus_margin"]) > 0 for row in rows),
            "positive" if rows else "missing",
        ),
        Check(
            "Summary rows include standard deviations",
            all(row["evaluation_seconds_sd"] != "" for row in summary_rows),
            "present" if summary_rows else "missing",
        ),
    ]
    failures = [check for check in checks if not check.passed]
    lines = [
        "# OpenFHE BGV Real PGS Audit",
        "",
        "## Verdict",
        "",
        "PASS: OpenFHE BGV fixed-point PRS benchmark completed."
        if not failures
        else "FAIL: OpenFHE BGV benchmark has unresolved issues.",
        "",
        "## Summary Table",
        "",
        "| SNPs | Scheme | Status | Eval sec mean | Max abs err mean | Input KB/sample mean | Note |",
        "| ---: | --- | --- | ---: | ---: | ---: | --- |",
    ]
    for row in summary_rows:
        input_kb = ""
        if row["input_ciphertext_mean_bytes_mean"]:
            input_kb = format_float(float(row["input_ciphertext_mean_bytes_mean"]) / 1024)
        lines.append(
            f"| {row['matched_snps']} | {row['scheme']} | {row['status']} | "
            f"{row['evaluation_seconds_mean']} | {row['max_abs_error_mean']} | {input_kb} | {row['note']} |"
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
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return checks


def main() -> int:
    args = parse_args()
    scoring = read_pgs_scoring_file(args.pgs)
    rows: list[dict[str, str]] = []
    for max_snps in args.max_snps:
        for repeat in range(1, args.repeats + 1):
            print(f"OpenFHE BGV PGS benchmark: snps={max_snps} repeat={repeat}/{args.repeats}", flush=True)
            match, match_mode, match_seconds = load_dataset(args, scoring, max_snps)
            rows.append(bgv_row(args, scoring, match, match_mode, match_seconds, max_snps, repeat))

    summary_rows = summarize_rows(rows)
    write_csv(args.raw_output, rows)
    write_csv(args.summary_output, summary_rows)
    checks = write_audit(args.report, rows, summary_rows, len(args.max_snps) * args.repeats)
    print(f"Wrote {args.raw_output}")
    print(f"Wrote {args.summary_output}")
    print(f"Wrote {args.report}")
    failures = [check for check in checks if not check.passed]
    if failures:
        print("OpenFHE BGV audit failures:")
        for check in failures:
            print(f"- {check.name}: {check.detail}")
        return 1
    print("OpenFHE BGV audit passed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
