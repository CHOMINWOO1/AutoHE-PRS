from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
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
from scripts.run_real_pgs_federated_packed_aggregate import (  # noqa: E402
    DEFAULT_CACHE,
    DEFAULT_PGS,
    DEFAULT_VCF,
    choose_block_size,
)


DEFAULT_OUTPUT = ROOT / "results" / "real_pgs_windowed_federated_packed_aggregate.csv"


@dataclass
class AggregateMetrics:
    site_sums: list[float]
    context_seconds: float = 0.0
    encryption_seconds: float = 0.0
    evaluation_seconds: float = 0.0
    decryption_seconds: float = 0.0
    serialization_seconds: float = 0.0
    public_context_bytes: int = 0
    input_ciphertext_bytes: int = 0
    result_ciphertext_bytes: int = 0
    blocks_total: int = 0
    calls_total: int = 0
    slots_per_ciphertext: int | str = ""
    poly_modulus_degree: int | str = ""
    scale: float | str = ""
    coeff_mod_bit_sizes: str = ""


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
        help="Matched SNP subset sizes for windowed real-data packed aggregate scaling.",
    )
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    return parser.parse_args()


def choose_window_snps(
    preferred_block_size: int,
    poly_modulus_degree: int = DEFAULT_POLY_MODULUS_DEGREE,
) -> int:
    if preferred_block_size <= 0:
        raise ValueError("preferred_block_size must be positive")
    slots = poly_modulus_degree // 2
    return max(1, slots // preferred_block_size)


def window_ranges(n_snps: int, window_snps: int) -> list[tuple[int, int]]:
    if n_snps <= 0:
        raise ValueError("n_snps must be positive")
    if window_snps <= 0:
        raise ValueError("window_snps must be positive")
    return [
        (start, min(start + window_snps, n_snps))
        for start in range(0, n_snps, window_snps)
    ]


def slice_genotype_window(
    genotypes: list[list[int]],
    start: int,
    end: int,
) -> list[list[int]]:
    return [row[start:end] for row in genotypes]


def baseline_site_metrics(
    site_genotypes: list[list[list[int]]],
    weights: list[float],
) -> AggregateMetrics:
    metrics = AggregateMetrics(site_sums=[])
    for site_rows in site_genotypes:
        baseline = encrypted_prs_benchmark(site_rows, weights)
        metrics.site_sums.append(sum(baseline.scores))
        metrics.context_seconds += baseline.context_seconds
        metrics.encryption_seconds += baseline.encryption_seconds
        metrics.evaluation_seconds += baseline.evaluation_seconds
        metrics.decryption_seconds += baseline.decryption_seconds
        metrics.serialization_seconds += baseline.serialization_seconds
        metrics.public_context_bytes += baseline.public_context_bytes
        metrics.input_ciphertext_bytes += baseline.input_ciphertext_bytes
        metrics.result_ciphertext_bytes += baseline.result_ciphertext_bytes
        metrics.calls_total += 1
        metrics.poly_modulus_degree = baseline.poly_modulus_degree
        metrics.scale = baseline.scale
        metrics.coeff_mod_bit_sizes = "-".join(str(value) for value in baseline.coeff_mod_bit_sizes)
    return metrics


def direct_aggregate_metrics(
    site_genotypes: list[list[list[int]]],
    weights: list[float],
    block_size: int,
) -> AggregateMetrics:
    metrics = AggregateMetrics(site_sums=[])
    for site_rows in site_genotypes:
        aggregate = packed_prs_aggregate_benchmark(site_rows, weights, block_size=block_size)
        metrics.site_sums.append(sum(aggregate.block_sums))
        metrics.context_seconds += aggregate.context_seconds
        metrics.encryption_seconds += aggregate.encryption_seconds
        metrics.evaluation_seconds += aggregate.evaluation_seconds
        metrics.decryption_seconds += aggregate.decryption_seconds
        metrics.serialization_seconds += aggregate.serialization_seconds
        metrics.public_context_bytes += aggregate.public_context_bytes
        metrics.input_ciphertext_bytes += aggregate.input_ciphertext_bytes
        metrics.result_ciphertext_bytes += aggregate.result_ciphertext_bytes
        metrics.blocks_total += aggregate.blocks
        metrics.calls_total += 1
        metrics.slots_per_ciphertext = aggregate.slots_per_ciphertext
        metrics.poly_modulus_degree = aggregate.poly_modulus_degree
        metrics.scale = aggregate.scale
        metrics.coeff_mod_bit_sizes = "-".join(str(value) for value in aggregate.coeff_mod_bit_sizes)
    return metrics


def windowed_aggregate_metrics(
    site_genotypes: list[list[list[int]]],
    weights: list[float],
    block_size: int,
    windows: list[tuple[int, int]],
) -> AggregateMetrics:
    metrics = AggregateMetrics(site_sums=[0.0 for _ in site_genotypes])
    for site_index, site_rows in enumerate(site_genotypes):
        for start, end in windows:
            aggregate = packed_prs_aggregate_benchmark(
                slice_genotype_window(site_rows, start, end),
                weights[start:end],
                block_size=block_size,
            )
            metrics.site_sums[site_index] += sum(aggregate.block_sums)
            metrics.context_seconds += aggregate.context_seconds
            metrics.encryption_seconds += aggregate.encryption_seconds
            metrics.evaluation_seconds += aggregate.evaluation_seconds
            metrics.decryption_seconds += aggregate.decryption_seconds
            metrics.serialization_seconds += aggregate.serialization_seconds
            metrics.public_context_bytes += aggregate.public_context_bytes
            metrics.input_ciphertext_bytes += aggregate.input_ciphertext_bytes
            metrics.result_ciphertext_bytes += aggregate.result_ciphertext_bytes
            metrics.blocks_total += aggregate.blocks
            metrics.calls_total += 1
            metrics.slots_per_ciphertext = aggregate.slots_per_ciphertext
            metrics.poly_modulus_degree = aggregate.poly_modulus_degree
            metrics.scale = aggregate.scale
            metrics.coeff_mod_bit_sizes = "-".join(str(value) for value in aggregate.coeff_mod_bit_sizes)
    return metrics


def max_abs_error(reference: list[float], observed: list[float]) -> float:
    return max(abs(plain - secure) for plain, secure in zip(reference, observed))


def max_mean_abs_error(reference: list[float], observed: list[float], sizes: list[int]) -> float:
    return max(
        abs(plain - secure) / size
        for plain, secure, size in zip(reference, observed, sizes)
    )


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
    direct_block_size = choose_block_size(dataset.n_snps, args.preferred_block_size)
    window_snps = choose_window_snps(args.preferred_block_size)
    windows = window_ranges(dataset.n_snps, window_snps)
    max_slots_used = max((end - start) * args.preferred_block_size for start, end in windows)

    site_genotypes = partition_rows(dataset.genotypes, args.sites)
    site_sizes = [len(site_rows) for site_rows in site_genotypes]

    start = perf_counter()
    plain_scores = prs_scores(dataset.genotypes, dataset.weights)
    plain_seconds = perf_counter() - start
    plain_site_sums = site_sums_from_scores(plain_scores, site_sizes)

    try:
        baseline = baseline_site_metrics(site_genotypes, dataset.weights)
        direct = direct_aggregate_metrics(site_genotypes, dataset.weights, direct_block_size)
        windowed = windowed_aggregate_metrics(
            site_genotypes,
            dataset.weights,
            args.preferred_block_size,
            windows,
        )
        status = "tenseal_ok"
    except TenSEALUnavailable as exc:
        status = f"tenseal_unavailable: {exc}"
        baseline = AggregateMetrics(site_sums=[])
        direct = AggregateMetrics(site_sums=[])
        windowed = AggregateMetrics(site_sums=[])

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
        "direct_block_size": direct_block_size,
        "windowed_block_size": args.preferred_block_size,
        "window_snps": window_snps,
        "window_count_per_site": len(windows),
        "window_count_total": len(windows) * len(site_genotypes),
        "max_slots_used_per_window": max_slots_used,
        "requested_weights_on_chrom": match.requested_weights_on_chrom,
        "skipped_allele_mismatch": match.skipped_allele_mismatch,
        "status": status,
        "match_mode": match_mode,
        "match_vcf_seconds": f"{match_seconds:.6f}",
        "plain_prs_seconds": f"{plain_seconds:.6f}",
        "plain_site_sums": format_float_list(plain_site_sums) if status == "tenseal_ok" else "",
        "baseline_site_sums": format_float_list(baseline.site_sums) if status == "tenseal_ok" else "",
        "direct_aggregate_site_sums": format_float_list(direct.site_sums) if status == "tenseal_ok" else "",
        "windowed_aggregate_site_sums": format_float_list(windowed.site_sums) if status == "tenseal_ok" else "",
        "baseline_evaluation_seconds_total": f"{baseline.evaluation_seconds:.6f}" if status == "tenseal_ok" else "",
        "direct_aggregate_evaluation_seconds_total": f"{direct.evaluation_seconds:.6f}" if status == "tenseal_ok" else "",
        "windowed_aggregate_evaluation_seconds_total": f"{windowed.evaluation_seconds:.6f}" if status == "tenseal_ok" else "",
        "baseline_context_seconds_total": f"{baseline.context_seconds:.6f}" if status == "tenseal_ok" else "",
        "direct_aggregate_context_seconds_total": f"{direct.context_seconds:.6f}" if status == "tenseal_ok" else "",
        "windowed_aggregate_context_seconds_total": f"{windowed.context_seconds:.6f}" if status == "tenseal_ok" else "",
        "baseline_encryption_seconds_total": f"{baseline.encryption_seconds:.6f}" if status == "tenseal_ok" else "",
        "direct_aggregate_encryption_seconds_total": f"{direct.encryption_seconds:.6f}" if status == "tenseal_ok" else "",
        "windowed_aggregate_encryption_seconds_total": f"{windowed.encryption_seconds:.6f}" if status == "tenseal_ok" else "",
        "baseline_decryption_seconds_total": f"{baseline.decryption_seconds:.6f}" if status == "tenseal_ok" else "",
        "direct_aggregate_decryption_seconds_total": f"{direct.decryption_seconds:.6f}" if status == "tenseal_ok" else "",
        "windowed_aggregate_decryption_seconds_total": f"{windowed.decryption_seconds:.6f}" if status == "tenseal_ok" else "",
        "baseline_serialization_seconds_total": f"{baseline.serialization_seconds:.6f}" if status == "tenseal_ok" else "",
        "direct_aggregate_serialization_seconds_total": f"{direct.serialization_seconds:.6f}" if status == "tenseal_ok" else "",
        "windowed_aggregate_serialization_seconds_total": f"{windowed.serialization_seconds:.6f}" if status == "tenseal_ok" else "",
        "baseline_site_sum_max_abs_error": f"{max_abs_error(plain_site_sums, baseline.site_sums):.12g}" if status == "tenseal_ok" else "",
        "direct_aggregate_site_sum_max_abs_error": f"{max_abs_error(plain_site_sums, direct.site_sums):.12g}" if status == "tenseal_ok" else "",
        "windowed_aggregate_site_sum_max_abs_error": f"{max_abs_error(plain_site_sums, windowed.site_sums):.12g}" if status == "tenseal_ok" else "",
        "windowed_aggregate_site_mean_max_abs_error": f"{max_mean_abs_error(plain_site_sums, windowed.site_sums, site_sizes):.12g}" if status == "tenseal_ok" else "",
        "baseline_input_ciphertext_bytes_total": baseline.input_ciphertext_bytes if status == "tenseal_ok" else "",
        "direct_aggregate_input_ciphertext_bytes_total": direct.input_ciphertext_bytes if status == "tenseal_ok" else "",
        "windowed_aggregate_input_ciphertext_bytes_total": windowed.input_ciphertext_bytes if status == "tenseal_ok" else "",
        "baseline_result_ciphertext_bytes_total": baseline.result_ciphertext_bytes if status == "tenseal_ok" else "",
        "direct_aggregate_result_ciphertext_bytes_total": direct.result_ciphertext_bytes if status == "tenseal_ok" else "",
        "windowed_aggregate_result_ciphertext_bytes_total": windowed.result_ciphertext_bytes if status == "tenseal_ok" else "",
        "direct_aggregate_blocks_total": direct.blocks_total if status == "tenseal_ok" else "",
        "windowed_aggregate_blocks_total": windowed.blocks_total if status == "tenseal_ok" else "",
        "windowed_aggregate_calls_total": windowed.calls_total if status == "tenseal_ok" else "",
        "windowed_eval_speedup_vs_baseline": ratio(
            baseline.evaluation_seconds,
            windowed.evaluation_seconds,
        ) if status == "tenseal_ok" else "",
        "windowed_eval_speedup_vs_direct_aggregate": ratio(
            direct.evaluation_seconds,
            windowed.evaluation_seconds,
        ) if status == "tenseal_ok" else "",
        "windowed_input_reduction_vs_baseline": ratio(
            baseline.input_ciphertext_bytes,
            windowed.input_ciphertext_bytes,
        ) if status == "tenseal_ok" else "",
        "windowed_result_reduction_vs_baseline": ratio(
            baseline.result_ciphertext_bytes,
            windowed.result_ciphertext_bytes,
        ) if status == "tenseal_ok" else "",
        "slots_per_ciphertext": windowed.slots_per_ciphertext if status == "tenseal_ok" else "",
        "poly_modulus_degree": windowed.poly_modulus_degree if status == "tenseal_ok" else "",
        "scale": str(windowed.scale) if status == "tenseal_ok" else "",
        "coeff_mod_bit_sizes": windowed.coeff_mod_bit_sizes if status == "tenseal_ok" else "",
        "security_model": (
            "real chr22 PGS-matched dosages partitioned across sites and SNP windows; "
            "site decrypts its own sum of packed aggregate PRS windows"
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
            "Real CAD windowed federated packed aggregate:",
            f"max_snps={max_snps}",
            f"matched={row['matched_snps']}",
            f"direct_block={row['direct_block_size']}",
            f"windows={row['window_count_per_site']}",
            f"windowed_eval={row['windowed_aggregate_evaluation_seconds_total']}",
            f"speedup_vs_baseline={row['windowed_eval_speedup_vs_baseline']}",
            f"speedup_vs_direct={row['windowed_eval_speedup_vs_direct_aggregate']}",
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
