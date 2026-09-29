from __future__ import annotations

import argparse
import csv
from collections import defaultdict
from dataclasses import dataclass
from importlib.metadata import version
from pathlib import Path
import platform
import statistics
import sys
from time import perf_counter


ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from genome_he_pilot.openfhe_controlled_backend import (  # noqa: E402
    DEFAULT_OPENFHE_CKKS_SCALING_MOD_SIZE,
    DEFAULT_OPENFHE_FIXED_POINT_SCALE,
    DEFAULT_OPENFHE_PLAIN_MODULUS,
    DEFAULT_OPENFHE_RING_DIM,
    DEFAULT_OPENFHE_WINDOW_SNPS,
    openfhe_ckks_prs_benchmark,
    openfhe_fixed_point_prs_benchmark,
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
DEFAULT_RAW_OUTPUT = ROOT / "results" / "real_pgs_chr22_pgs004941_openfhe_controlled_raw.csv"
DEFAULT_SUMMARY_OUTPUT = ROOT / "results" / "real_pgs_chr22_pgs004941_openfhe_controlled_summary.csv"
DEFAULT_REPORT = ROOT / "submission" / "openfhe_controlled_scheme_comparison_audit.md"

SCHEMES = ("CKKS", "BFV", "BGV")
SCHEME_ORDER = {scheme: index for index, scheme in enumerate(SCHEMES)}
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
    parser.add_argument("--max-snps", type=int, nargs="+", default=[512, 1024, 2048, 4096, 8192, 16384, 32768, 51227])
    parser.add_argument("--window-snps", type=int, default=DEFAULT_OPENFHE_WINDOW_SNPS)
    parser.add_argument("--ring-dim", type=int, default=DEFAULT_OPENFHE_RING_DIM)
    parser.add_argument("--plain-modulus", type=int, default=DEFAULT_OPENFHE_PLAIN_MODULUS)
    parser.add_argument("--fixed-point-scale", type=int, default=DEFAULT_OPENFHE_FIXED_POINT_SCALE)
    parser.add_argument("--ckks-scaling-mod-size", type=int, default=DEFAULT_OPENFHE_CKKS_SCALING_MOD_SIZE)
    parser.add_argument("--raw-output", type=Path, default=DEFAULT_RAW_OUTPUT)
    parser.add_argument("--summary-output", type=Path, default=DEFAULT_SUMMARY_OUTPUT)
    parser.add_argument("--report", type=Path, default=DEFAULT_REPORT)
    return parser.parse_args()


def format_float(value: float | int | None) -> str:
    if value is None:
        return ""
    return f"{float(value):.12g}"


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


def run_scheme(args: argparse.Namespace, dataset, scheme: str):
    if scheme == "CKKS":
        return openfhe_ckks_prs_benchmark(
            dataset.genotypes,
            dataset.weights,
            window_snps=args.window_snps,
            ring_dim=args.ring_dim,
            scaling_mod_size=args.ckks_scaling_mod_size,
        )
    return openfhe_fixed_point_prs_benchmark(
        dataset.genotypes,
        dataset.weights,
        scheme=scheme,
        window_snps=args.window_snps,
        ring_dim=args.ring_dim,
        plain_modulus=args.plain_modulus,
        fixed_point_scale=args.fixed_point_scale,
    )


def row_for_result(
    args: argparse.Namespace,
    scoring,
    match,
    match_mode: str,
    match_seconds: float,
    max_snps: int,
    repeat: int,
    scheme: str,
) -> dict[str, str]:
    dataset = match.dataset
    start = perf_counter()
    plain_scores = prs_scores(dataset.genotypes, dataset.weights)
    plain_seconds = perf_counter() - start

    result = run_scheme(args, dataset, scheme)
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
        "scheme": scheme,
        "backend": f"OpenFHE {scheme}RNS" if scheme != "CKKS" else "OpenFHE CKKSRNS",
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
        "ring_dim": str(result.ring_dim),
        "batch_size": str(result.batch_size),
        "plain_modulus": "" if result.plain_modulus is None else str(result.plain_modulus),
        "fixed_point_scale": "" if result.fixed_point_scale is None else str(result.fixed_point_scale),
        "ckks_scaling_mod_size": "" if result.ckks_scaling_mod_size is None else str(result.ckks_scaling_mod_size),
        "integer_weight_abs_sum": "" if result.integer_weight_abs_sum is None else str(result.integer_weight_abs_sum),
        "integer_dot_abs_bound": "" if result.integer_dot_abs_bound is None else str(result.integer_dot_abs_bound),
        "safe_modulus_margin": "" if result.safe_modulus_margin is None else str(result.safe_modulus_margin),
        "openfhe_version": version("openfhe"),
        "python_version": platform.python_version(),
        "platform": platform.platform(),
        "note": "Same-backend OpenFHE controlled PRS comparison",
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
    groups: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        groups[(row["max_snps"], row["scheme"])].append(row)
    summary_rows: list[dict[str, str]] = []
    for (max_snps, scheme), group in sorted(groups.items(), key=lambda item: (int(item[0][0]), SCHEME_ORDER[item[0][1]])):
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
            "scheme": scheme,
            "backend": first["backend"],
            "ring_dim": first["ring_dim"],
            "batch_size": first["batch_size"],
            "openfhe_version": first["openfhe_version"],
            "python_version": first["python_version"],
            "platform": first["platform"],
            "repeats": str(len(group)),
            "ok_repeats": str(sum(1 for row in group if row["status"] == "ok")),
            "status": "ok" if all(row["status"] == "ok" for row in group) else "incomplete",
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
                for suffix in ("mean", "sd", "min", "max"):
                    summary[f"{metric}_{suffix}"] = ""
        summary_rows.append(summary)
    return summary_rows


def passfail(value: bool) -> str:
    return "PASS" if value else "FAIL"


def scheme_sets(rows: list[dict[str, str]]) -> dict[str, set[str]]:
    found: dict[str, set[str]] = defaultdict(set)
    for row in rows:
        found[row["max_snps"]].add(row["scheme"])
    return dict(found)


def write_report(path: Path, rows: list[dict[str, str]], summary_rows: list[dict[str, str]], args: argparse.Namespace) -> list[Check]:
    sets = scheme_sets(summary_rows)
    expected_sizes = {str(size) for size in args.max_snps}
    checks = [
        Check("All requested SNP sizes are summarized", set(sets) == expected_sizes, ",".join(sorted(sets, key=int))),
        Check("Each SNP size has CKKS, BFV, and BGV", all(sets.get(str(size)) == set(SCHEMES) for size in args.max_snps), str(sets)),
        Check("All controlled rows are ok", all(row["status"] == "ok" for row in summary_rows), "all ok"),
        Check("All rows use OpenFHE 1.5.1 package", all(row["openfhe_version"].startswith("1.5.1") for row in summary_rows), summary_rows[0]["openfhe_version"] if summary_rows else "missing"),
        Check("All rows use the same ring dimension", {row["ring_dim"] for row in summary_rows} == {str(args.ring_dim)}, str(sorted({row["ring_dim"] for row in summary_rows}))),
        Check("All errors remain below 1e-3", all(float(row["max_abs_error_mean"]) < 1e-3 for row in summary_rows), f"max={max(float(row['max_abs_error_mean']) for row in summary_rows):.3g}"),
    ]
    failures = [check for check in checks if not check.passed]
    lines = [
        "# OpenFHE Controlled CKKS/BFV/BGV Scheme Comparison Audit",
        "",
        "## Verdict",
        "",
        "PASS: CKKS, BFV, and BGV were benchmarked in the same OpenFHE Docker environment."
        if not failures
        else "FAIL: controlled OpenFHE scheme comparison has unresolved issues.",
        "",
        "## Controlled Environment",
        "",
        f"- Docker image: `genome-he-openfhe:1.5.1`",
        f"- OpenFHE-Python: `{summary_rows[0]['openfhe_version'] if summary_rows else 'missing'}`",
        f"- Python: `{summary_rows[0]['python_version'] if summary_rows else 'missing'}`",
        f"- Platform: `{summary_rows[0]['platform'] if summary_rows else 'missing'}`",
        f"- Ring dimension: `{args.ring_dim}`",
        f"- Window SNPs: `{args.window_snps}`",
        f"- Fixed-point scale for BFV/BGV: `{args.fixed_point_scale}`",
        f"- Plain modulus for BFV/BGV: `{args.plain_modulus}`",
        f"- CKKS scaling modulus size: `{args.ckks_scaling_mod_size}`",
        "",
        "## Summary Table",
        "",
        "| SNPs | Scheme | Windows/sample | Eval sec mean | Max abs err mean | Input MB/sample mean |",
        "| ---: | --- | ---: | ---: | ---: | ---: |",
    ]
    for row in summary_rows:
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
            match, match_mode, match_seconds = load_dataset(args, scoring, max_snps)
            for scheme in SCHEMES:
                print(f"OpenFHE controlled {scheme}: snps={max_snps} repeat={repeat}/{args.repeats}", flush=True)
                row = row_for_result(args, scoring, match, match_mode, match_seconds, max_snps, repeat, scheme)
                rows.append(row)
                print(f"  eval={row['evaluation_seconds']} err={row['max_abs_error']}", flush=True)
    write_csv(args.raw_output, rows)
    summary_rows = summarize_rows(rows)
    write_csv(args.summary_output, summary_rows)
    checks = write_report(args.report, rows, summary_rows, args)
    print(f"Wrote {args.raw_output}")
    print(f"Wrote {args.summary_output}")
    print(f"Wrote {args.report}")
    failures = [check for check in checks if not check.passed]
    if failures:
        print("OpenFHE controlled audit failures:")
        for check in failures:
            print(f"- {check.name}: {check.detail}")
        return 1
    print("OpenFHE controlled scheme comparison audit passed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
