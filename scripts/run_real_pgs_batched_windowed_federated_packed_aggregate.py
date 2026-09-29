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
from genome_he_pilot.packed_ckks import packed_windowed_prs_aggregate_benchmark  # noqa: E402
from genome_he_pilot.pgs import read_pgs_scoring_file  # noqa: E402
from genome_he_pilot.pgs_vcf import load_pgs_matched_vcf_dataset  # noqa: E402
from genome_he_pilot.tenseal_backend import TenSEALUnavailable  # noqa: E402
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
from scripts.run_real_pgs_windowed_federated_packed_aggregate import (  # noqa: E402
    AggregateMetrics,
    baseline_site_metrics,
    choose_window_snps,
    direct_aggregate_metrics,
    max_abs_error,
    max_mean_abs_error,
    window_ranges,
    windowed_aggregate_metrics,
)


DEFAULT_OUTPUT = ROOT / "results" / "real_pgs_batched_windowed_federated_packed_aggregate.csv"


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
        help="Matched SNP subset sizes for batched windowed real-data packed aggregate scaling.",
    )
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    return parser.parse_args()


def metric_compute_seconds(metrics: AggregateMetrics) -> float:
    return (
        metrics.context_seconds
        + metrics.encryption_seconds
        + metrics.evaluation_seconds
        + metrics.decryption_seconds
    )


def batched_compute_seconds(result) -> float:
    return (
        result.context_seconds
        + result.encryption_seconds
        + result.evaluation_seconds
        + result.decryption_seconds
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
        independent = windowed_aggregate_metrics(
            site_genotypes,
            dataset.weights,
            args.preferred_block_size,
            windows,
        )
        batched = packed_windowed_prs_aggregate_benchmark(
            site_genotypes,
            dataset.weights,
            block_size=args.preferred_block_size,
            window_snps=window_snps,
        )
        status = "tenseal_ok"
    except TenSEALUnavailable as exc:
        baseline = AggregateMetrics(site_sums=[])
        direct = AggregateMetrics(site_sums=[])
        independent = AggregateMetrics(site_sums=[])
        batched = None
        status = f"tenseal_unavailable: {exc}"

    if status == "tenseal_ok" and batched is not None:
        independent_compute = metric_compute_seconds(independent)
        batched_compute = batched_compute_seconds(batched)
        coeff_mod_bit_sizes = "-".join(str(value) for value in batched.coeff_mod_bit_sizes)
        batched_site_sums = batched.site_sums
    else:
        independent_compute = 0.0
        batched_compute = 0.0
        coeff_mod_bit_sizes = ""
        batched_site_sums = []

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
        "batched_windowed_block_size": args.preferred_block_size,
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
        "independent_windowed_site_sums": format_float_list(independent.site_sums) if status == "tenseal_ok" else "",
        "batched_windowed_site_sums": format_float_list(batched_site_sums) if status == "tenseal_ok" else "",
        "baseline_evaluation_seconds_total": f"{baseline.evaluation_seconds:.6f}" if status == "tenseal_ok" else "",
        "direct_aggregate_evaluation_seconds_total": f"{direct.evaluation_seconds:.6f}" if status == "tenseal_ok" else "",
        "independent_windowed_evaluation_seconds_total": f"{independent.evaluation_seconds:.6f}" if status == "tenseal_ok" else "",
        "batched_windowed_evaluation_seconds_total": f"{batched.evaluation_seconds:.6f}" if batched is not None else "",
        "baseline_context_seconds_total": f"{baseline.context_seconds:.6f}" if status == "tenseal_ok" else "",
        "direct_aggregate_context_seconds_total": f"{direct.context_seconds:.6f}" if status == "tenseal_ok" else "",
        "independent_windowed_context_seconds_total": f"{independent.context_seconds:.6f}" if status == "tenseal_ok" else "",
        "batched_windowed_context_seconds_total": f"{batched.context_seconds:.6f}" if batched is not None else "",
        "baseline_encryption_seconds_total": f"{baseline.encryption_seconds:.6f}" if status == "tenseal_ok" else "",
        "direct_aggregate_encryption_seconds_total": f"{direct.encryption_seconds:.6f}" if status == "tenseal_ok" else "",
        "independent_windowed_encryption_seconds_total": f"{independent.encryption_seconds:.6f}" if status == "tenseal_ok" else "",
        "batched_windowed_encryption_seconds_total": f"{batched.encryption_seconds:.6f}" if batched is not None else "",
        "baseline_decryption_seconds_total": f"{baseline.decryption_seconds:.6f}" if status == "tenseal_ok" else "",
        "direct_aggregate_decryption_seconds_total": f"{direct.decryption_seconds:.6f}" if status == "tenseal_ok" else "",
        "independent_windowed_decryption_seconds_total": f"{independent.decryption_seconds:.6f}" if status == "tenseal_ok" else "",
        "batched_windowed_decryption_seconds_total": f"{batched.decryption_seconds:.6f}" if batched is not None else "",
        "independent_windowed_compute_seconds_total": f"{independent_compute:.6f}" if status == "tenseal_ok" else "",
        "batched_windowed_compute_seconds_total": f"{batched_compute:.6f}" if status == "tenseal_ok" else "",
        "independent_windowed_serialization_seconds_total": f"{independent.serialization_seconds:.6f}" if status == "tenseal_ok" else "",
        "batched_windowed_serialization_seconds_total": f"{batched.serialization_seconds:.6f}" if batched is not None else "",
        "baseline_site_sum_max_abs_error": f"{max_abs_error(plain_site_sums, baseline.site_sums):.12g}" if status == "tenseal_ok" else "",
        "direct_aggregate_site_sum_max_abs_error": f"{max_abs_error(plain_site_sums, direct.site_sums):.12g}" if status == "tenseal_ok" else "",
        "independent_windowed_site_sum_max_abs_error": f"{max_abs_error(plain_site_sums, independent.site_sums):.12g}" if status == "tenseal_ok" else "",
        "batched_windowed_site_sum_max_abs_error": f"{max_abs_error(plain_site_sums, batched_site_sums):.12g}" if status == "tenseal_ok" else "",
        "batched_windowed_site_mean_max_abs_error": f"{max_mean_abs_error(plain_site_sums, batched_site_sums, site_sizes):.12g}" if status == "tenseal_ok" else "",
        "independent_windowed_public_context_bytes_total": independent.public_context_bytes if status == "tenseal_ok" else "",
        "batched_windowed_public_context_bytes_total": batched.public_context_bytes if batched is not None else "",
        "independent_windowed_input_ciphertext_bytes_total": independent.input_ciphertext_bytes if status == "tenseal_ok" else "",
        "batched_windowed_input_ciphertext_bytes_total": batched.input_ciphertext_bytes if batched is not None else "",
        "independent_windowed_result_ciphertext_bytes_total": independent.result_ciphertext_bytes if status == "tenseal_ok" else "",
        "batched_windowed_result_ciphertext_bytes_total": batched.result_ciphertext_bytes if batched is not None else "",
        "independent_windowed_blocks_total": independent.blocks_total if status == "tenseal_ok" else "",
        "batched_windowed_blocks_total": batched.blocks if batched is not None else "",
        "independent_windowed_calls_total": independent.calls_total if status == "tenseal_ok" else "",
        "batched_windowed_calls_total": batched.calls if batched is not None else "",
        "batched_windowed_ciphertexts_total": batched.ciphertexts if batched is not None else "",
        "batched_eval_speedup_vs_baseline": ratio(
            baseline.evaluation_seconds,
            batched.evaluation_seconds,
        ) if batched is not None else "",
        "batched_eval_speedup_vs_direct_aggregate": ratio(
            direct.evaluation_seconds,
            batched.evaluation_seconds,
        ) if batched is not None else "",
        "batched_eval_speedup_vs_independent_windowed": ratio(
            independent.evaluation_seconds,
            batched.evaluation_seconds,
        ) if batched is not None else "",
        "batched_compute_speedup_vs_independent_windowed": ratio(
            independent_compute,
            batched_compute,
        ) if batched is not None else "",
        "batched_context_reduction_vs_independent_windowed": ratio(
            independent.context_seconds,
            batched.context_seconds,
        ) if batched is not None else "",
        "batched_public_context_reduction_vs_independent_windowed": ratio(
            independent.public_context_bytes,
            batched.public_context_bytes,
        ) if batched is not None else "",
        "slots_per_ciphertext": batched.slots_per_ciphertext if batched is not None else "",
        "poly_modulus_degree": batched.poly_modulus_degree if batched is not None else "",
        "scale": str(batched.scale) if batched is not None else "",
        "coeff_mod_bit_sizes": coeff_mod_bit_sizes,
        "security_model": (
            "real chr22 PGS-matched dosages partitioned across sites and SNP windows; "
            "one reusable public CKKS context batches all site-window aggregate ciphertexts"
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
            "Real CAD batched windowed federated packed aggregate:",
            f"max_snps={max_snps}",
            f"matched={row['matched_snps']}",
            f"windows={row['window_count_per_site']}",
            f"batched_eval={row['batched_windowed_evaluation_seconds_total']}",
            f"eval_speedup_vs_independent={row['batched_eval_speedup_vs_independent_windowed']}",
            f"compute_speedup_vs_independent={row['batched_compute_speedup_vs_independent_windowed']}",
            f"context_reduction={row['batched_context_reduction_vs_independent_windowed']}",
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
