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

from genome_he_pilot.bfv_backend import (  # noqa: E402
    DEFAULT_BFV_FIXED_POINT_SCALE,
    bfv_fixed_point_prs_benchmark,
)
from genome_he_pilot.plaintext import prs_scores  # noqa: E402
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
DEFAULT_RAW_OUTPUT = ROOT / "results" / "real_pgs_chr22_pgs004941_scheme_comparison_raw.csv"
DEFAULT_SUMMARY_OUTPUT = ROOT / "results" / "real_pgs_chr22_pgs004941_scheme_comparison_summary.csv"
DEFAULT_REPORT = ROOT / "submission" / "real_pgs_scheme_comparison_audit.md"

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
    parser.add_argument("--max-snps", type=int, nargs="+", default=[512, 1024, 2048])
    parser.add_argument("--bfv-fixed-point-scale", type=int, default=DEFAULT_BFV_FIXED_POINT_SCALE)
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
    match_seconds = perf_counter() - start
    return match, match_mode, match_seconds


def base_row(args: argparse.Namespace, scoring, match, match_mode: str, match_seconds: float, max_snps: int, repeat: int) -> dict[str, str]:
    dataset = match.dataset
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
        "match_mode": match_mode,
        "match_vcf_seconds": format_float(match_seconds),
    }


def scheme_row(
    base: dict[str, str],
    *,
    scheme: str,
    backend: str,
    status: str,
    plain_seconds: float,
    context_seconds: float = 0.0,
    encryption_seconds: float = 0.0,
    evaluation_seconds: float = 0.0,
    decryption_seconds: float = 0.0,
    serialization_seconds: float = 0.0,
    max_abs_error: float = 0.0,
    public_context_bytes: int = 0,
    input_ciphertext_bytes: int = 0,
    result_ciphertext_bytes: int = 0,
    poly_modulus_degree: str = "",
    scale: str = "",
    coeff_mod_bit_sizes: str = "",
    plain_modulus: str = "",
    fixed_point_scale: str = "",
    integer_dot_abs_bound: str = "",
    safe_modulus_margin: str = "",
    note: str = "",
) -> dict[str, str]:
    samples = int(base["samples"])
    snps = int(base["matched_snps"])
    total_he_seconds = encryption_seconds + evaluation_seconds + decryption_seconds
    eval_overhead = evaluation_seconds / plain_seconds if plain_seconds > 0 else 0.0
    total_overhead = total_he_seconds / plain_seconds if plain_seconds > 0 else 0.0
    he_samples_per_second = samples / evaluation_seconds if evaluation_seconds > 0 else 0.0
    he_snp_products_per_second = samples * snps / evaluation_seconds if evaluation_seconds > 0 else 0.0
    input_mean = input_ciphertext_bytes / samples if samples else 0.0
    result_mean = result_ciphertext_bytes / samples if samples else 0.0

    row = dict(base)
    row.update(
        {
            "scheme": scheme,
            "backend": backend,
            "status": status,
            "plain_prs_seconds": format_float(plain_seconds),
            "context_seconds": format_float(context_seconds),
            "encryption_seconds": format_float(encryption_seconds),
            "evaluation_seconds": format_float(evaluation_seconds),
            "decryption_seconds": format_float(decryption_seconds),
            "serialization_seconds": format_float(serialization_seconds),
            "total_he_compute_seconds": format_float(total_he_seconds),
            "max_abs_error": format_float(max_abs_error),
            "eval_overhead_ratio": format_float(eval_overhead),
            "total_overhead_ratio": format_float(total_overhead),
            "he_samples_per_second": format_float(he_samples_per_second),
            "he_snp_products_per_second": format_float(he_snp_products_per_second),
            "input_ciphertext_bytes": str(input_ciphertext_bytes),
            "input_ciphertext_mean_bytes": format_float(input_mean),
            "result_ciphertext_bytes": str(result_ciphertext_bytes),
            "result_ciphertext_mean_bytes": format_float(result_mean),
            "public_context_bytes": str(public_context_bytes),
            "poly_modulus_degree": poly_modulus_degree,
            "scale": scale,
            "coeff_mod_bit_sizes": coeff_mod_bit_sizes,
            "plain_modulus": plain_modulus,
            "fixed_point_scale": fixed_point_scale,
            "integer_dot_abs_bound": integer_dot_abs_bound,
            "safe_modulus_margin": safe_modulus_margin,
            "note": note,
        }
    )
    return row


def run_once(args: argparse.Namespace, scoring, max_snps: int, repeat: int) -> list[dict[str, str]]:
    match, match_mode, match_seconds = load_dataset(args, scoring, max_snps)
    dataset = match.dataset
    base = base_row(args, scoring, match, match_mode, match_seconds, max_snps, repeat)

    start = perf_counter()
    plain_scores = prs_scores(dataset.genotypes, dataset.weights)
    plain_seconds = perf_counter() - start

    rows: list[dict[str, str]] = []

    try:
        ckks = encrypted_prs_benchmark(dataset.genotypes, dataset.weights)
        ckks_error = max(abs(plain - encrypted) for plain, encrypted in zip(plain_scores, ckks.scores))
        rows.append(
            scheme_row(
                base,
                scheme="CKKS",
                backend="TenSEAL CKKSVector",
                status="ok",
                plain_seconds=plain_seconds,
                context_seconds=ckks.context_seconds,
                encryption_seconds=ckks.encryption_seconds,
                evaluation_seconds=ckks.evaluation_seconds,
                decryption_seconds=ckks.decryption_seconds,
                serialization_seconds=ckks.serialization_seconds,
                max_abs_error=ckks_error,
                public_context_bytes=ckks.public_context_bytes,
                input_ciphertext_bytes=ckks.input_ciphertext_bytes,
                result_ciphertext_bytes=ckks.result_ciphertext_bytes,
                poly_modulus_degree=str(ckks.poly_modulus_degree),
                scale=str(ckks.scale),
                coeff_mod_bit_sizes="-".join(str(value) for value in ckks.coeff_mod_bit_sizes),
                note="Approximate real-valued PRS dot product",
            )
        )
    except TenSEALUnavailable as exc:
        rows.append(
            scheme_row(
                base,
                scheme="CKKS",
                backend="TenSEAL CKKSVector",
                status=f"unavailable: {exc}",
                plain_seconds=plain_seconds,
            )
        )

    try:
        bfv = bfv_fixed_point_prs_benchmark(
            dataset.genotypes,
            dataset.weights,
            fixed_point_scale=args.bfv_fixed_point_scale,
        )
        bfv_error = max(abs(plain - encrypted) for plain, encrypted in zip(plain_scores, bfv.scores))
        rows.append(
            scheme_row(
                base,
                scheme="BFV",
                backend="TenSEAL BFVVector",
                status="ok",
                plain_seconds=plain_seconds,
                context_seconds=bfv.context_seconds,
                encryption_seconds=bfv.encryption_seconds,
                evaluation_seconds=bfv.evaluation_seconds,
                decryption_seconds=bfv.decryption_seconds,
                serialization_seconds=bfv.serialization_seconds,
                max_abs_error=bfv_error,
                public_context_bytes=bfv.public_context_bytes,
                input_ciphertext_bytes=bfv.input_ciphertext_bytes,
                result_ciphertext_bytes=bfv.result_ciphertext_bytes,
                poly_modulus_degree=str(bfv.poly_modulus_degree),
                plain_modulus=str(bfv.plain_modulus),
                fixed_point_scale=str(bfv.fixed_point_scale),
                integer_dot_abs_bound=str(bfv.integer_dot_abs_bound),
                safe_modulus_margin=str(bfv.safe_modulus_margin),
                note="Integer fixed-point PRS dot product; error includes weight quantization",
            )
        )
    except (TenSEALUnavailable, ValueError) as exc:
        rows.append(
            scheme_row(
                base,
                scheme="BFV",
                backend="TenSEAL BFVVector",
                status=f"unavailable: {exc}",
                plain_seconds=plain_seconds,
                fixed_point_scale=str(args.bfv_fixed_point_scale),
            )
        )

    rows.append(
        scheme_row(
            base,
            scheme="BGV",
            backend="not available in TenSEAL 0.3.16",
            status="not_run_backend_unavailable",
            plain_seconds=plain_seconds,
            note="TenSEAL exposes CKKS and BFV only in this environment; BGV requires another validated backend",
        )
    )
    return rows


def write_csv(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        raise ValueError("cannot write an empty CSV")
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def summarize_rows(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    groups: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    for row in rows:
        groups[(row["max_snps"], row["scheme"])].append(row)

    summary_rows: list[dict[str, str]] = []
    for (max_snps, scheme) in sorted(groups, key=lambda item: (int(item[0]), item[1])):
        group = groups[(max_snps, scheme)]
        first = group[0]
        ok_group = [row for row in group if row["status"] == "ok"]
        summary: dict[str, str] = {
            "pgs_id": first["pgs_id"],
            "chrom": first["chrom"],
            "samples": first["samples"],
            "max_snps": first["max_snps"],
            "matched_snps": first["matched_snps"],
            "scheme": scheme,
            "backend": first["backend"],
            "repeats": str(len(group)),
            "ok_repeats": str(len(ok_group)),
            "status": "ok" if ok_group else first["status"],
            "note": first["note"],
        }
        for metric in SUMMARY_METRICS:
            values = [float(value) for row in ok_group if (value := float_or_none(row.get(metric, ""))) is not None]
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
    ckks_rows = [row for row in rows if row["scheme"] == "CKKS"]
    bfv_rows = [row for row in rows if row["scheme"] == "BFV"]
    bgv_rows = [row for row in rows if row["scheme"] == "BGV"]
    ok_summaries = [row for row in summary_rows if row["status"] == "ok"]
    checks = [
        Check(
            "Raw rows include CKKS, BFV, and BGV status for each repeat",
            len(rows) == args.repeats * len(args.max_snps) * 3,
            f"rows={len(rows)} expected={args.repeats * len(args.max_snps) * 3}",
        ),
        Check("CKKS runs completed", all(row["status"] == "ok" for row in ckks_rows), "all ok"),
        Check("BFV fixed-point runs completed", all(row["status"] == "ok" for row in bfv_rows), "all ok"),
        Check(
            "BGV is explicitly marked unavailable rather than reported as fake runtime",
            all(row["status"] == "not_run_backend_unavailable" for row in bgv_rows),
            "not_run_backend_unavailable",
        ),
        Check("CKKS error remains below 1e-4", max(float(row["max_abs_error"]) for row in ckks_rows) < 1e-4, f"max={max(float(row['max_abs_error']) for row in ckks_rows):.3g}"),
        Check("BFV fixed-point error remains below 1e-3", max(float(row["max_abs_error"]) for row in bfv_rows) < 1e-3, f"max={max(float(row['max_abs_error']) for row in bfv_rows):.3g}"),
        Check("BFV modulus margins are positive", all(float(row["safe_modulus_margin"]) > 0 for row in bfv_rows), "positive"),
        Check("Summary rows include standard deviations for completed schemes", all(row["evaluation_seconds_sd"] != "" for row in ok_summaries), "present"),
    ]
    return checks


def passfail(value: bool) -> str:
    return "PASS" if value else "FAIL"


def write_report(path: Path, checks: list[Check], summary_rows: list[dict[str, str]]) -> None:
    failures = [check for check in checks if not check.passed]
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = [
        "# Real PGS HE Scheme Comparison Audit",
        "",
        "## Verdict",
        "",
        "PASS: CKKS/BFV comparison completed and BGV backend unavailability is explicitly recorded."
        if not failures
        else "FAIL: HE scheme comparison has unresolved issues.",
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
            "| "
            + " | ".join(
                [
                    row["matched_snps"],
                    row["scheme"],
                    row["status"],
                    row["evaluation_seconds_mean"],
                    row["max_abs_error_mean"],
                    input_kb,
                    row["note"],
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
            print(f"HE scheme comparison: snps={max_snps} repeat={repeat}/{args.repeats}", flush=True)
            run_rows = run_once(args, scoring, max_snps, repeat)
            rows.extend(run_rows)
            for row in run_rows:
                print(
                    "  ",
                    row["scheme"],
                    row["status"],
                    f"eval={row['evaluation_seconds']}",
                    f"err={row['max_abs_error']}",
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
        print("HE scheme comparison audit failures:")
        for check in failures:
            print(f"- {check.name}: {check.detail}")
        return 1
    print("HE scheme comparison audit passed.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
