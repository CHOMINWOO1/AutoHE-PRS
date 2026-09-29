from __future__ import annotations

import argparse
import csv
from pathlib import Path
import sys
from time import perf_counter

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from genome_he_pilot import prs_scores  # noqa: E402
from genome_he_pilot.packed_ckks import packed_prs_aggregate_benchmark  # noqa: E402
from genome_he_pilot.pgs import read_pgs_scoring_file  # noqa: E402
from genome_he_pilot.pgs_vcf import load_pgs_matched_vcf_dataset  # noqa: E402
from genome_he_pilot.tenseal_backend import (  # noqa: E402
    DEFAULT_POLY_MODULUS_DEGREE,
    TenSEALUnavailable,
    encrypted_prs_benchmark,
)
from genome_he_pilot.vcf_cache import load_pgs_matched_cached_dataset  # noqa: E402
from scripts.run_federated_packed_aggregate import (  # noqa: E402
    format_float_list,
    partition_rows,
    ratio,
    site_sums_from_scores,
)


DEFAULT_VCF = (
    ROOT
    / "data"
    / "1000genomes"
    / "phase3"
    / "ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"
)
DEFAULT_PGS = ROOT / "data" / "pgs_catalog" / "PGS004941" / "PGS004941.txt.gz"
DEFAULT_CACHE = ROOT / "data" / "cache" / "1000g_chr22_samples32.sqlite"
DEFAULT_OUTPUT = ROOT / "results" / "real_pgs_federated_packed_aggregate.csv"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--pgs", type=Path, default=DEFAULT_PGS)
    parser.add_argument("--vcf", type=Path, default=DEFAULT_VCF)
    parser.add_argument("--chrom", default="22")
    parser.add_argument("--samples", type=int, default=32)
    parser.add_argument("--sites", type=int, default=4)
    parser.add_argument("--preferred-block-size", type=int, default=4)
    parser.add_argument("--cache", type=Path, default=DEFAULT_CACHE)
    parser.add_argument("--no-cache", action="store_true")
    parser.add_argument(
        "--max-snps",
        type=int,
        nargs="+",
        default=[512, 1024, 2048, 4096],
        help="Matched SNP subset sizes for real-data federated packed aggregate scaling.",
    )
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    return parser.parse_args()


def choose_block_size(
    snps: int,
    preferred_block_size: int,
    poly_modulus_degree: int = DEFAULT_POLY_MODULUS_DEGREE,
) -> int:
    if snps <= 0:
        raise ValueError("snps must be positive")
    if preferred_block_size <= 0:
        raise ValueError("preferred_block_size must be positive")
    slots = poly_modulus_degree // 2
    return max(1, min(preferred_block_size, slots // snps))


def run_row(args: argparse.Namespace, scoring, max_snps: int) -> dict[str, str | int]:
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
    block_size = choose_block_size(dataset.n_snps, args.preferred_block_size)
    site_genotypes = partition_rows(dataset.genotypes, args.sites)
    site_sizes = [len(site_rows) for site_rows in site_genotypes]

    start = perf_counter()
    plain_scores = prs_scores(dataset.genotypes, dataset.weights)
    plain_seconds = perf_counter() - start
    plain_site_sums = site_sums_from_scores(plain_scores, site_sizes)

    baseline_site_sums: list[float] = []
    aggregate_site_sums: list[float] = []
    baseline_context_seconds = 0.0
    baseline_encryption_seconds = 0.0
    baseline_evaluation_seconds = 0.0
    baseline_decryption_seconds = 0.0
    baseline_serialization_seconds = 0.0
    baseline_public_context_bytes = 0
    baseline_input_ciphertext_bytes = 0
    baseline_result_ciphertext_bytes = 0
    aggregate_context_seconds = 0.0
    aggregate_encryption_seconds = 0.0
    aggregate_evaluation_seconds = 0.0
    aggregate_decryption_seconds = 0.0
    aggregate_serialization_seconds = 0.0
    aggregate_public_context_bytes = 0
    aggregate_input_ciphertext_bytes = 0
    aggregate_result_ciphertext_bytes = 0
    aggregate_blocks_total = 0
    slots_per_ciphertext: int | str = ""
    poly_modulus_degree: int | str = ""
    scale: float | str = ""
    coeff_mod_bit_sizes = ""

    try:
        for site_rows in site_genotypes:
            baseline = encrypted_prs_benchmark(site_rows, dataset.weights)
            baseline_site_sums.append(sum(baseline.scores))
            baseline_context_seconds += baseline.context_seconds
            baseline_encryption_seconds += baseline.encryption_seconds
            baseline_evaluation_seconds += baseline.evaluation_seconds
            baseline_decryption_seconds += baseline.decryption_seconds
            baseline_serialization_seconds += baseline.serialization_seconds
            baseline_public_context_bytes += baseline.public_context_bytes
            baseline_input_ciphertext_bytes += baseline.input_ciphertext_bytes
            baseline_result_ciphertext_bytes += baseline.result_ciphertext_bytes
            poly_modulus_degree = baseline.poly_modulus_degree
            scale = baseline.scale
            coeff_mod_bit_sizes = "-".join(str(value) for value in baseline.coeff_mod_bit_sizes)

        for site_rows in site_genotypes:
            aggregate = packed_prs_aggregate_benchmark(
                site_rows,
                dataset.weights,
                block_size=block_size,
            )
            aggregate_site_sums.append(sum(aggregate.block_sums))
            aggregate_context_seconds += aggregate.context_seconds
            aggregate_encryption_seconds += aggregate.encryption_seconds
            aggregate_evaluation_seconds += aggregate.evaluation_seconds
            aggregate_decryption_seconds += aggregate.decryption_seconds
            aggregate_serialization_seconds += aggregate.serialization_seconds
            aggregate_public_context_bytes += aggregate.public_context_bytes
            aggregate_input_ciphertext_bytes += aggregate.input_ciphertext_bytes
            aggregate_result_ciphertext_bytes += aggregate.result_ciphertext_bytes
            aggregate_blocks_total += aggregate.blocks
            slots_per_ciphertext = aggregate.slots_per_ciphertext

        baseline_site_sum_max_abs_error: float | str = max(
            abs(plain - secure)
            for plain, secure in zip(plain_site_sums, baseline_site_sums)
        )
        aggregate_site_sum_max_abs_error: float | str = max(
            abs(plain - secure)
            for plain, secure in zip(plain_site_sums, aggregate_site_sums)
        )
        aggregate_site_mean_max_abs_error: float | str = max(
            abs(plain - secure) / size
            for plain, secure, size in zip(plain_site_sums, aggregate_site_sums, site_sizes)
        )
        status = "tenseal_ok"
    except TenSEALUnavailable as exc:
        status = f"tenseal_unavailable: {exc}"
        baseline_site_sum_max_abs_error = ""
        aggregate_site_sum_max_abs_error = ""
        aggregate_site_mean_max_abs_error = ""

    return {
        "pgs_id": scoring.metadata.get("pgs_id", ""),
        "trait_reported": scoring.metadata.get("trait_reported", ""),
        "genome_build": scoring.metadata.get("genome_build", ""),
        "chrom": args.chrom,
        "sites": args.sites,
        "samples": dataset.n_samples,
        "site_sample_counts": ";".join(str(size) for size in site_sizes),
        "matched_snps": dataset.n_snps,
        "max_snps": max_snps,
        "preferred_block_size": args.preferred_block_size,
        "block_size": block_size,
        "aggregate_slots_used_per_block": block_size * dataset.n_snps,
        "requested_weights_on_chrom": match.requested_weights_on_chrom,
        "skipped_allele_mismatch": match.skipped_allele_mismatch,
        "status": status,
        "match_mode": match_mode,
        "match_vcf_seconds": f"{match_seconds:.6f}",
        "plain_prs_seconds": f"{plain_seconds:.6f}",
        "plain_site_sums": format_float_list(plain_site_sums) if status == "tenseal_ok" else "",
        "baseline_site_sums": format_float_list(baseline_site_sums) if status == "tenseal_ok" else "",
        "aggregate_site_sums": format_float_list(aggregate_site_sums) if status == "tenseal_ok" else "",
        "baseline_context_seconds_total": f"{baseline_context_seconds:.6f}" if status == "tenseal_ok" else "",
        "baseline_encryption_seconds_total": f"{baseline_encryption_seconds:.6f}" if status == "tenseal_ok" else "",
        "baseline_evaluation_seconds_total": f"{baseline_evaluation_seconds:.6f}" if status == "tenseal_ok" else "",
        "baseline_decryption_seconds_total": f"{baseline_decryption_seconds:.6f}" if status == "tenseal_ok" else "",
        "baseline_serialization_seconds_total": f"{baseline_serialization_seconds:.6f}" if status == "tenseal_ok" else "",
        "baseline_site_sum_max_abs_error": (
            f"{baseline_site_sum_max_abs_error:.12g}"
            if isinstance(baseline_site_sum_max_abs_error, float)
            else ""
        ),
        "baseline_public_context_bytes_total": baseline_public_context_bytes if status == "tenseal_ok" else "",
        "baseline_input_ciphertext_bytes_total": baseline_input_ciphertext_bytes if status == "tenseal_ok" else "",
        "baseline_result_ciphertext_bytes_total": baseline_result_ciphertext_bytes if status == "tenseal_ok" else "",
        "aggregate_blocks_total": aggregate_blocks_total if status == "tenseal_ok" else "",
        "aggregate_context_seconds_total": f"{aggregate_context_seconds:.6f}" if status == "tenseal_ok" else "",
        "aggregate_encryption_seconds_total": f"{aggregate_encryption_seconds:.6f}" if status == "tenseal_ok" else "",
        "aggregate_evaluation_seconds_total": f"{aggregate_evaluation_seconds:.6f}" if status == "tenseal_ok" else "",
        "aggregate_decryption_seconds_total": f"{aggregate_decryption_seconds:.6f}" if status == "tenseal_ok" else "",
        "aggregate_serialization_seconds_total": f"{aggregate_serialization_seconds:.6f}" if status == "tenseal_ok" else "",
        "aggregate_site_sum_max_abs_error": (
            f"{aggregate_site_sum_max_abs_error:.12g}"
            if isinstance(aggregate_site_sum_max_abs_error, float)
            else ""
        ),
        "aggregate_site_mean_max_abs_error": (
            f"{aggregate_site_mean_max_abs_error:.12g}"
            if isinstance(aggregate_site_mean_max_abs_error, float)
            else ""
        ),
        "aggregate_public_context_bytes_total": aggregate_public_context_bytes if status == "tenseal_ok" else "",
        "aggregate_input_ciphertext_bytes_total": aggregate_input_ciphertext_bytes if status == "tenseal_ok" else "",
        "aggregate_result_ciphertext_bytes_total": aggregate_result_ciphertext_bytes if status == "tenseal_ok" else "",
        "aggregate_eval_speedup_vs_baseline": ratio(
            baseline_evaluation_seconds,
            aggregate_evaluation_seconds,
        ) if status == "tenseal_ok" else "",
        "aggregate_input_ciphertext_reduction": ratio(
            baseline_input_ciphertext_bytes,
            aggregate_input_ciphertext_bytes,
        ) if status == "tenseal_ok" else "",
        "aggregate_result_ciphertext_reduction": ratio(
            baseline_result_ciphertext_bytes,
            aggregate_result_ciphertext_bytes,
        ) if status == "tenseal_ok" else "",
        "slots_per_ciphertext": slots_per_ciphertext if status == "tenseal_ok" else "",
        "poly_modulus_degree": poly_modulus_degree if status == "tenseal_ok" else "",
        "scale": str(scale) if status == "tenseal_ok" else "",
        "coeff_mod_bit_sizes": coeff_mod_bit_sizes if status == "tenseal_ok" else "",
        "security_model": (
            "real chr22 PGS-matched dosages partitioned across sites; "
            "site decrypts its own packed aggregate PRS sum"
        ),
    }


def main() -> int:
    args = parse_args()
    args.output.parent.mkdir(parents=True, exist_ok=True)

    start = perf_counter()
    scoring = read_pgs_scoring_file(args.pgs, chrom=args.chrom)
    read_pgs_seconds = perf_counter() - start

    rows = []
    for max_snps in args.max_snps:
        row = run_row(args, scoring, max_snps)
        row["read_pgs_seconds"] = f"{read_pgs_seconds:.6f}"
        rows.append(row)
        print(
            "Real CAD federated packed aggregate:",
            f"max_snps={max_snps}",
            f"matched={row['matched_snps']}",
            f"block={row['block_size']}",
            f"baseline_eval={row['baseline_evaluation_seconds_total']}",
            f"aggregate_eval={row['aggregate_evaluation_seconds_total']}",
            f"speedup={row['aggregate_eval_speedup_vs_baseline']}",
            f"status={row['status']}",
        )

    with args.output.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    print(f"Wrote {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
