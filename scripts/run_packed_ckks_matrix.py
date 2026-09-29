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
from genome_he_pilot.packed_ckks import (  # noqa: E402
    packed_prs_aggregate_benchmark,
    packed_prs_benchmark,
)
from genome_he_pilot.tenseal_backend import (  # noqa: E402
    TenSEALUnavailable,
    encrypted_prs_benchmark,
)


EXPERIMENTS = [
    {"name": "PACK1_smoke", "samples": 8, "snps": 128, "block_size": 2},
    {"name": "PACK2_small", "samples": 16, "snps": 256, "block_size": 4},
    {"name": "PACK3_slot_probe", "samples": 32, "snps": 512, "block_size": 4},
]


def run_one(
    name: str,
    samples: int,
    snps: int,
    block_size: int,
    seed: int,
) -> dict[str, str | int]:
    start = perf_counter()
    dataset = make_synthetic_dataset(n_samples=samples, n_snps=snps, seed=seed)
    data_seconds = perf_counter() - start

    start = perf_counter()
    plain_scores = prs_scores(dataset.genotypes, dataset.weights)
    plain_seconds = perf_counter() - start

    try:
        baseline = encrypted_prs_benchmark(dataset.genotypes, dataset.weights)
        packed = packed_prs_benchmark(
            dataset.genotypes,
            dataset.weights,
            block_size=block_size,
        )
        aggregate = packed_prs_aggregate_benchmark(
            dataset.genotypes,
            dataset.weights,
            block_size=block_size,
        )
        max_abs_error = max(
            abs(plain - secure)
            for plain, secure in zip(plain_scores, packed.scores)
        )
        plain_block_sums = [
            sum(plain_scores[index:index + block_size])
            for index in range(0, len(plain_scores), block_size)
        ]
        aggregate_max_abs_error = max(
            abs(plain - secure)
            for plain, secure in zip(plain_block_sums, aggregate.block_sums)
        )
        status = "ok"
    except TenSEALUnavailable as exc:
        status = f"tenseal_unavailable: {exc}"
        baseline = None
        packed = None
        aggregate = None
        max_abs_error = ""
        aggregate_max_abs_error = ""

    return {
        "experiment": name,
        "samples": samples,
        "snps": snps,
        "block_size": block_size,
        "seed": seed,
        "status": status,
        "data_seconds": f"{data_seconds:.6f}",
        "plain_prs_seconds": f"{plain_seconds:.6f}",
        "baseline_evaluation_seconds": f"{baseline.evaluation_seconds:.6f}" if baseline else "",
        "baseline_input_ciphertext_bytes": baseline.input_ciphertext_bytes if baseline else "",
        "baseline_result_ciphertext_bytes": baseline.result_ciphertext_bytes if baseline else "",
        "packed_blocks": packed.blocks if packed else "",
        "packed_context_seconds": f"{packed.context_seconds:.6f}" if packed else "",
        "packed_encryption_seconds": f"{packed.encryption_seconds:.6f}" if packed else "",
        "packed_evaluation_seconds": f"{packed.evaluation_seconds:.6f}" if packed else "",
        "packed_decryption_seconds": f"{packed.decryption_seconds:.6f}" if packed else "",
        "packed_serialization_seconds": f"{packed.serialization_seconds:.6f}" if packed else "",
        "packed_max_abs_error": f"{max_abs_error:.12g}" if isinstance(max_abs_error, float) else "",
        "packed_public_context_bytes": packed.public_context_bytes if packed else "",
        "packed_input_ciphertext_bytes": packed.input_ciphertext_bytes if packed else "",
        "packed_result_ciphertext_bytes": packed.result_ciphertext_bytes if packed else "",
        "packed_slots_used": block_size * snps,
        "aggregate_blocks": aggregate.blocks if aggregate else "",
        "aggregate_context_seconds": f"{aggregate.context_seconds:.6f}" if aggregate else "",
        "aggregate_encryption_seconds": f"{aggregate.encryption_seconds:.6f}" if aggregate else "",
        "aggregate_evaluation_seconds": f"{aggregate.evaluation_seconds:.6f}" if aggregate else "",
        "aggregate_decryption_seconds": f"{aggregate.decryption_seconds:.6f}" if aggregate else "",
        "aggregate_serialization_seconds": f"{aggregate.serialization_seconds:.6f}" if aggregate else "",
        "aggregate_max_abs_error": f"{aggregate_max_abs_error:.12g}" if isinstance(aggregate_max_abs_error, float) else "",
        "aggregate_public_context_bytes": aggregate.public_context_bytes if aggregate else "",
        "aggregate_input_ciphertext_bytes": aggregate.input_ciphertext_bytes if aggregate else "",
        "aggregate_result_ciphertext_bytes": aggregate.result_ciphertext_bytes if aggregate else "",
        "slots_per_ciphertext": packed.slots_per_ciphertext if packed else "",
        "poly_modulus_degree": packed.poly_modulus_degree if packed else "",
        "scale": str(packed.scale) if packed else "",
        "coeff_mod_bit_sizes": "-".join(str(value) for value in packed.coeff_mod_bit_sizes) if packed else "",
    }


def main() -> int:
    output = ROOT / "results" / "packed_ckks_prs_matrix.csv"
    output.parent.mkdir(parents=True, exist_ok=True)

    rows = []
    for index, experiment in enumerate(EXPERIMENTS):
        print(
            "Running",
            experiment["name"],
            f"samples={experiment['samples']}",
            f"snps={experiment['snps']}",
            f"block={experiment['block_size']}",
            flush=True,
        )
        rows.append(run_one(seed=20260810 + index, **experiment))

    with output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    print(f"Wrote {output}")
    for row in rows:
        print(
            row["experiment"],
            f"status={row['status']}",
            f"baseline_eval={row['baseline_evaluation_seconds']}",
            f"packed_eval={row['packed_evaluation_seconds']}",
            f"aggregate_eval={row['aggregate_evaluation_seconds']}",
            f"packed_error={row['packed_max_abs_error']}",
            f"aggregate_error={row['aggregate_max_abs_error']}",
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
