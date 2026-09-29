from __future__ import annotations

import csv
from pathlib import Path
import sys
from time import perf_counter

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from genome_he_pilot import make_synthetic_dataset, prs_scores  # noqa: E402
from genome_he_pilot.packed_ckks import packed_prs_aggregate_benchmark  # noqa: E402
from genome_he_pilot.tenseal_backend import (  # noqa: E402
    TenSEALUnavailable,
    encrypted_prs_benchmark,
)


EXPERIMENTS = [
    {
        "name": "SCI_FED_PACK_AGG_1_smoke",
        "sites": 2,
        "samples_per_site": 8,
        "snps": 256,
        "block_size": 4,
    },
    {
        "name": "SCI_FED_PACK_AGG_2_kci_plus",
        "sites": 4,
        "samples_per_site": 8,
        "snps": 512,
        "block_size": 4,
    },
    {
        "name": "SCI_FED_PACK_AGG_3_scale_probe",
        "sites": 4,
        "samples_per_site": 16,
        "snps": 1024,
        "block_size": 4,
    },
]


def partition_rows(rows: list[list[int]], sites: int) -> list[list[list[int]]]:
    if sites <= 0:
        raise ValueError("sites must be positive")
    chunk_size = (len(rows) + sites - 1) // sites
    return [rows[index * chunk_size:(index + 1) * chunk_size] for index in range(sites)]


def site_sums_from_scores(scores: list[float], site_sizes: list[int]) -> list[float]:
    sums: list[float] = []
    offset = 0
    for size in site_sizes:
        sums.append(sum(scores[offset:offset + size]))
        offset += size
    if offset != len(scores):
        raise ValueError("site sizes must account for all scores")
    return sums


def format_float_list(values: list[float]) -> str:
    return ";".join(f"{value:.12g}" for value in values)


def ratio(numerator: float | int, denominator: float | int) -> str:
    if denominator == 0:
        return ""
    return f"{float(numerator) / float(denominator):.6f}"


def run_one(
    name: str,
    sites: int,
    samples_per_site: int,
    snps: int,
    block_size: int,
    seed: int,
) -> dict[str, str | int]:
    total_samples = sites * samples_per_site

    start = perf_counter()
    dataset = make_synthetic_dataset(
        n_samples=total_samples,
        n_snps=snps,
        seed=seed,
    )
    data_seconds = perf_counter() - start

    start = perf_counter()
    plain_scores = prs_scores(dataset.genotypes, dataset.weights)
    plain_seconds = perf_counter() - start

    site_genotypes = partition_rows(dataset.genotypes, sites)
    site_sizes = [len(site_rows) for site_rows in site_genotypes]
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
        status = "ok"
    except TenSEALUnavailable as exc:
        status = f"tenseal_unavailable: {exc}"
        baseline_site_sum_max_abs_error = ""
        aggregate_site_sum_max_abs_error = ""
        aggregate_site_mean_max_abs_error = ""

    return {
        "experiment": name,
        "sites": sites,
        "samples_per_site": samples_per_site,
        "total_samples": total_samples,
        "snps": snps,
        "block_size": block_size,
        "seed": seed,
        "status": status,
        "data_seconds": f"{data_seconds:.6f}",
        "plain_prs_seconds": f"{plain_seconds:.6f}",
        "site_sample_counts": ";".join(str(size) for size in site_sizes),
        "plain_site_sums": format_float_list(plain_site_sums) if status == "ok" else "",
        "baseline_site_sums": format_float_list(baseline_site_sums) if status == "ok" else "",
        "aggregate_site_sums": format_float_list(aggregate_site_sums) if status == "ok" else "",
        "baseline_context_seconds_total": f"{baseline_context_seconds:.6f}" if status == "ok" else "",
        "baseline_encryption_seconds_total": f"{baseline_encryption_seconds:.6f}" if status == "ok" else "",
        "baseline_evaluation_seconds_total": f"{baseline_evaluation_seconds:.6f}" if status == "ok" else "",
        "baseline_decryption_seconds_total": f"{baseline_decryption_seconds:.6f}" if status == "ok" else "",
        "baseline_serialization_seconds_total": f"{baseline_serialization_seconds:.6f}" if status == "ok" else "",
        "baseline_site_sum_max_abs_error": (
            f"{baseline_site_sum_max_abs_error:.12g}"
            if isinstance(baseline_site_sum_max_abs_error, float)
            else ""
        ),
        "baseline_public_context_bytes_total": baseline_public_context_bytes if status == "ok" else "",
        "baseline_input_ciphertext_bytes_total": baseline_input_ciphertext_bytes if status == "ok" else "",
        "baseline_result_ciphertext_bytes_total": baseline_result_ciphertext_bytes if status == "ok" else "",
        "aggregate_blocks_total": aggregate_blocks_total if status == "ok" else "",
        "aggregate_context_seconds_total": f"{aggregate_context_seconds:.6f}" if status == "ok" else "",
        "aggregate_encryption_seconds_total": f"{aggregate_encryption_seconds:.6f}" if status == "ok" else "",
        "aggregate_evaluation_seconds_total": f"{aggregate_evaluation_seconds:.6f}" if status == "ok" else "",
        "aggregate_decryption_seconds_total": f"{aggregate_decryption_seconds:.6f}" if status == "ok" else "",
        "aggregate_serialization_seconds_total": f"{aggregate_serialization_seconds:.6f}" if status == "ok" else "",
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
        "aggregate_public_context_bytes_total": aggregate_public_context_bytes if status == "ok" else "",
        "aggregate_input_ciphertext_bytes_total": aggregate_input_ciphertext_bytes if status == "ok" else "",
        "aggregate_result_ciphertext_bytes_total": aggregate_result_ciphertext_bytes if status == "ok" else "",
        "aggregate_eval_speedup_vs_baseline": ratio(
            baseline_evaluation_seconds,
            aggregate_evaluation_seconds,
        ) if status == "ok" else "",
        "aggregate_input_ciphertext_reduction": ratio(
            baseline_input_ciphertext_bytes,
            aggregate_input_ciphertext_bytes,
        ) if status == "ok" else "",
        "aggregate_result_ciphertext_reduction": ratio(
            baseline_result_ciphertext_bytes,
            aggregate_result_ciphertext_bytes,
        ) if status == "ok" else "",
        "aggregate_slots_used_per_block": block_size * snps,
        "slots_per_ciphertext": slots_per_ciphertext if status == "ok" else "",
        "poly_modulus_degree": poly_modulus_degree if status == "ok" else "",
        "scale": str(scale) if status == "ok" else "",
        "coeff_mod_bit_sizes": coeff_mod_bit_sizes if status == "ok" else "",
        "security_model": (
            "per-site CKKS context; site sends encrypted packed genotype blocks; "
            "evaluator computes encrypted block PRS sums; site decrypts site aggregate"
        ),
    }


def main() -> int:
    output = ROOT / "results" / "sci_federated_packed_aggregate.csv"
    output.parent.mkdir(parents=True, exist_ok=True)

    rows = []
    for index, experiment in enumerate(EXPERIMENTS):
        print(
            "Running",
            experiment["name"],
            f"sites={experiment['sites']}",
            f"samples/site={experiment['samples_per_site']}",
            f"snps={experiment['snps']}",
            f"block={experiment['block_size']}",
            flush=True,
        )
        rows.append(run_one(seed=20260820 + index, **experiment))

    with output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    print(f"Wrote {output}")
    for row in rows:
        print(
            row["experiment"],
            f"status={row['status']}",
            f"sites={row['sites']}",
            f"snps={row['snps']}",
            f"baseline_eval={row['baseline_evaluation_seconds_total']}",
            f"aggregate_eval={row['aggregate_evaluation_seconds_total']}",
            f"speedup={row['aggregate_eval_speedup_vs_baseline']}",
            f"aggregate_error={row['aggregate_site_sum_max_abs_error']}",
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
