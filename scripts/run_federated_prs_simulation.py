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
from genome_he_pilot.tenseal_backend import (  # noqa: E402
    TenSEALUnavailable,
    encrypted_prs_benchmark,
)


EXPERIMENTS = [
    {"name": "SCI_FED_1_smoke", "sites": 2, "samples_per_site": 8, "snps": 256},
    {"name": "SCI_FED_2_kci_plus", "sites": 4, "samples_per_site": 8, "snps": 512},
    {"name": "SCI_FED_3_scale_probe", "sites": 4, "samples_per_site": 16, "snps": 1024},
]


def partition_rows(rows: list[list[int]], sites: int) -> list[list[list[int]]]:
    if sites <= 0:
        raise ValueError("sites must be positive")
    chunk_size = (len(rows) + sites - 1) // sites
    return [rows[index * chunk_size:(index + 1) * chunk_size] for index in range(sites)]


def run_one(
    name: str,
    sites: int,
    samples_per_site: int,
    snps: int,
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
    encrypted_scores: list[float] = []
    context_seconds = 0.0
    encryption_seconds = 0.0
    evaluation_seconds = 0.0
    decryption_seconds = 0.0
    serialization_seconds = 0.0
    public_context_bytes = 0
    input_ciphertext_bytes = 0
    result_ciphertext_bytes = 0
    poly_modulus_degree = ""
    scale = ""
    coeff_mod_bit_sizes = ""

    try:
        for site_rows in site_genotypes:
            encrypted = encrypted_prs_benchmark(site_rows, dataset.weights)
            encrypted_scores.extend(encrypted.scores)
            context_seconds += encrypted.context_seconds
            encryption_seconds += encrypted.encryption_seconds
            evaluation_seconds += encrypted.evaluation_seconds
            decryption_seconds += encrypted.decryption_seconds
            serialization_seconds += encrypted.serialization_seconds
            public_context_bytes += encrypted.public_context_bytes
            input_ciphertext_bytes += encrypted.input_ciphertext_bytes
            result_ciphertext_bytes += encrypted.result_ciphertext_bytes
            poly_modulus_degree = str(encrypted.poly_modulus_degree)
            scale = str(encrypted.scale)
            coeff_mod_bit_sizes = "-".join(str(value) for value in encrypted.coeff_mod_bit_sizes)

        max_abs_error = max(
            abs(plain - secure)
            for plain, secure in zip(plain_scores, encrypted_scores)
        )
        status = "ok"
    except TenSEALUnavailable as exc:
        status = f"tenseal_unavailable: {exc}"
        max_abs_error = ""

    return {
        "experiment": name,
        "sites": sites,
        "samples_per_site": samples_per_site,
        "total_samples": total_samples,
        "snps": snps,
        "seed": seed,
        "status": status,
        "data_seconds": f"{data_seconds:.6f}",
        "plain_prs_seconds": f"{plain_seconds:.6f}",
        "context_seconds_total": f"{context_seconds:.6f}" if status == "ok" else "",
        "encryption_seconds_total": f"{encryption_seconds:.6f}" if status == "ok" else "",
        "evaluation_seconds_total": f"{evaluation_seconds:.6f}" if status == "ok" else "",
        "decryption_seconds_total": f"{decryption_seconds:.6f}" if status == "ok" else "",
        "serialization_seconds_total": f"{serialization_seconds:.6f}" if status == "ok" else "",
        "max_abs_error": f"{max_abs_error:.12g}" if isinstance(max_abs_error, float) else "",
        "public_context_bytes_total": public_context_bytes if status == "ok" else "",
        "input_ciphertext_bytes_total": input_ciphertext_bytes if status == "ok" else "",
        "result_ciphertext_bytes_total": result_ciphertext_bytes if status == "ok" else "",
        "poly_modulus_degree": poly_modulus_degree if status == "ok" else "",
        "scale": scale if status == "ok" else "",
        "coeff_mod_bit_sizes": coeff_mod_bit_sizes if status == "ok" else "",
        "security_model": "per-site CKKS context; evaluator sees encrypted genotypes only",
    }


def main() -> int:
    output = ROOT / "results" / "sci_federated_prs_simulation.csv"
    output.parent.mkdir(parents=True, exist_ok=True)

    rows = []
    for index, experiment in enumerate(EXPERIMENTS):
        print(
            "Running",
            experiment["name"],
            f"sites={experiment['sites']}",
            f"samples/site={experiment['samples_per_site']}",
            f"snps={experiment['snps']}",
            flush=True,
        )
        rows.append(run_one(seed=20260731 + index, **experiment))

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
            f"eval_total={row['evaluation_seconds_total']}",
            f"error={row['max_abs_error']}",
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
