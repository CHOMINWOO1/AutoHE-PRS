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
from genome_he_pilot.memory import (  # noqa: E402
    process_memory,
    stop_process_memory_measurement,
)
from genome_he_pilot.tenseal_backend import (  # noqa: E402
    TenSEALUnavailable,
    encrypted_prs_benchmark,
)


EXPERIMENTS = [
    {"name": "HE1_validation", "samples": 8, "snps": 128},
    {"name": "HE2_small", "samples": 16, "snps": 256},
    {"name": "HE3_kci_table", "samples": 32, "snps": 512},
    {"name": "HE4_medium", "samples": 64, "snps": 1024},
    {"name": "HE5_large_pilot", "samples": 96, "snps": 2048},
]


def run_one(name: str, samples: int, snps: int, seed: int) -> dict[str, str | int]:
    memory_before = process_memory()

    start = perf_counter()
    dataset = make_synthetic_dataset(n_samples=samples, n_snps=snps, seed=seed)
    data_seconds = perf_counter() - start

    start = perf_counter()
    plain_scores = prs_scores(dataset.genotypes, dataset.weights)
    plain_prs_seconds = perf_counter() - start

    try:
        encrypted = encrypted_prs_benchmark(dataset.genotypes, dataset.weights)
        max_abs_error = max(
            abs(plain - secure)
            for plain, secure in zip(plain_scores, encrypted.scores)
        )
        status = "ok"
        context_seconds = f"{encrypted.context_seconds:.6f}"
        encryption_seconds = f"{encrypted.encryption_seconds:.6f}"
        evaluation_seconds = f"{encrypted.evaluation_seconds:.6f}"
        decryption_seconds = f"{encrypted.decryption_seconds:.6f}"
        serialization_seconds = f"{encrypted.serialization_seconds:.6f}"
        public_context_bytes = encrypted.public_context_bytes
        input_ciphertext_bytes = encrypted.input_ciphertext_bytes
        result_ciphertext_bytes = encrypted.result_ciphertext_bytes
        mean_input_ciphertext_bytes = encrypted.input_ciphertext_bytes // samples
        mean_result_ciphertext_bytes = encrypted.result_ciphertext_bytes // samples
        poly_modulus_degree = encrypted.poly_modulus_degree
        scale = str(encrypted.scale)
        coeff_mod_bit_sizes = "-".join(str(value) for value in encrypted.coeff_mod_bit_sizes)
    except TenSEALUnavailable as exc:
        status = f"tenseal_unavailable: {exc}"
        max_abs_error = ""
        context_seconds = ""
        encryption_seconds = ""
        evaluation_seconds = ""
        decryption_seconds = ""
        serialization_seconds = ""
        public_context_bytes = ""
        input_ciphertext_bytes = ""
        result_ciphertext_bytes = ""
        mean_input_ciphertext_bytes = ""
        mean_result_ciphertext_bytes = ""
        poly_modulus_degree = ""
        scale = ""
        coeff_mod_bit_sizes = ""

    memory = stop_process_memory_measurement(memory_before)

    return {
        "experiment": name,
        "samples": samples,
        "snps": snps,
        "seed": seed,
        "status": status,
        "data_seconds": f"{data_seconds:.6f}",
        "plain_prs_seconds": f"{plain_prs_seconds:.6f}",
        "context_seconds": context_seconds,
        "encryption_seconds": encryption_seconds,
        "evaluation_seconds": evaluation_seconds,
        "decryption_seconds": decryption_seconds,
        "serialization_seconds": serialization_seconds,
        "max_abs_error": f"{max_abs_error:.12g}" if isinstance(max_abs_error, float) else "",
        "public_context_bytes": public_context_bytes,
        "input_ciphertext_bytes": input_ciphertext_bytes,
        "mean_input_ciphertext_bytes": mean_input_ciphertext_bytes,
        "result_ciphertext_bytes": result_ciphertext_bytes,
        "mean_result_ciphertext_bytes": mean_result_ciphertext_bytes,
        "poly_modulus_degree": poly_modulus_degree,
        "scale": scale,
        "coeff_mod_bit_sizes": coeff_mod_bit_sizes,
        "rss_before_bytes": memory.rss_before_bytes,
        "rss_after_bytes": memory.rss_after_bytes,
        "rss_delta_bytes": memory.rss_delta_bytes,
        "peak_rss_bytes": memory.peak_rss_bytes,
        "python_peak_bytes": memory.python_peak_bytes,
    }


def main() -> int:
    output = ROOT / "results" / "encrypted_experiment_matrix.csv"
    output.parent.mkdir(parents=True, exist_ok=True)

    rows = []
    for index, experiment in enumerate(EXPERIMENTS):
        print(
            "Running",
            experiment["name"],
            f"samples={experiment['samples']}",
            f"snps={experiment['snps']}",
            flush=True,
        )
        rows.append(run_one(seed=20260711 + index, **experiment))

    with output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    print(f"Wrote {output}")
    for row in rows:
        print(
            row["experiment"],
            f"status={row['status']}",
            f"plain={row['plain_prs_seconds']}",
            f"enc={row['encryption_seconds']}",
            f"eval={row['evaluation_seconds']}",
            f"dec={row['decryption_seconds']}",
            f"error={row['max_abs_error']}",
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
