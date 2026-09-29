from __future__ import annotations

import csv
from pathlib import Path
import sys
from time import perf_counter

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from genome_he_pilot import (  # noqa: E402
    allele_counts,
    load_vcf_dosage_subset,
    prs_scores,
    summarize_scores,
)
from genome_he_pilot.memory import (  # noqa: E402
    process_memory,
    stop_process_memory_measurement,
)
from genome_he_pilot.tenseal_backend import (  # noqa: E402
    TenSEALUnavailable,
    encrypted_prs_benchmark,
)


VCF = (
    ROOT
    / "data"
    / "1000genomes"
    / "phase3"
    / "ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"
)

EXPERIMENTS = [
    {"name": "G1_plain_validation", "samples": 32, "snps": 512, "backend": "plaintext"},
    {"name": "G2_plain_public", "samples": 128, "snps": 2048, "backend": "plaintext"},
    {"name": "G3_tenseal_small", "samples": 16, "snps": 256, "backend": "tenseal"},
    {"name": "G4_tenseal_kci", "samples": 32, "snps": 512, "backend": "tenseal"},
]


def run_one(
    name: str,
    samples: int,
    snps: int,
    backend: str,
    seed: int,
    min_minor_allele_count: int,
) -> dict[str, str | int]:
    memory_before = process_memory()

    start = perf_counter()
    dataset = load_vcf_dosage_subset(
        VCF,
        n_samples=samples,
        n_snps=snps,
        seed=seed,
        min_minor_allele_count=min_minor_allele_count,
    )
    load_seconds = perf_counter() - start

    start = perf_counter()
    plain_scores = prs_scores(dataset.genotypes, dataset.weights)
    plain_prs_seconds = perf_counter() - start
    summary = summarize_scores(plain_scores)

    start = perf_counter()
    counts = allele_counts(dataset.genotypes)
    count_seconds = perf_counter() - start

    status = "plaintext_only"
    context_seconds = ""
    encryption_seconds = ""
    evaluation_seconds = ""
    decryption_seconds = ""
    max_abs_error = ""
    public_context_bytes = ""
    input_ciphertext_bytes = ""
    result_ciphertext_bytes = ""
    poly_modulus_degree = ""
    scale = ""
    coeff_mod_bit_sizes = ""

    if backend == "tenseal":
        try:
            encrypted = encrypted_prs_benchmark(dataset.genotypes, dataset.weights)
            max_abs_error = max(
                abs(plain - secure)
                for plain, secure in zip(plain_scores, encrypted.scores)
            )
            status = "tenseal_ok"
            context_seconds = f"{encrypted.context_seconds:.6f}"
            encryption_seconds = f"{encrypted.encryption_seconds:.6f}"
            evaluation_seconds = f"{encrypted.evaluation_seconds:.6f}"
            decryption_seconds = f"{encrypted.decryption_seconds:.6f}"
            public_context_bytes = encrypted.public_context_bytes
            input_ciphertext_bytes = encrypted.input_ciphertext_bytes
            result_ciphertext_bytes = encrypted.result_ciphertext_bytes
            poly_modulus_degree = encrypted.poly_modulus_degree
            scale = str(encrypted.scale)
            coeff_mod_bit_sizes = "-".join(str(value) for value in encrypted.coeff_mod_bit_sizes)
        except TenSEALUnavailable as exc:
            status = f"tenseal_unavailable: {exc}"

    first_variant = dataset.variants[0]
    last_variant = dataset.variants[-1]
    memory = stop_process_memory_measurement(memory_before)
    return {
        "experiment": name,
        "samples": dataset.n_samples,
        "snps": dataset.n_snps,
        "backend": backend,
        "status": status,
        "seed": seed,
        "min_minor_allele_count": min_minor_allele_count,
        "load_seconds": f"{load_seconds:.6f}",
        "plain_prs_seconds": f"{plain_prs_seconds:.6f}",
        "count_seconds": f"{count_seconds:.6f}",
        "context_seconds": context_seconds,
        "encryption_seconds": encryption_seconds,
        "evaluation_seconds": evaluation_seconds,
        "decryption_seconds": decryption_seconds,
        "max_abs_error": f"{max_abs_error:.12g}" if isinstance(max_abs_error, float) else "",
        "public_context_bytes": public_context_bytes,
        "input_ciphertext_bytes": input_ciphertext_bytes,
        "result_ciphertext_bytes": result_ciphertext_bytes,
        "poly_modulus_degree": poly_modulus_degree,
        "scale": scale,
        "coeff_mod_bit_sizes": coeff_mod_bit_sizes,
        "prs_mean": f"{summary.mean:.8f}",
        "prs_sd": f"{summary.standard_deviation:.8f}",
        "prs_min": f"{summary.minimum:.8f}",
        "prs_max": f"{summary.maximum:.8f}",
        "first_snp_allele_count": counts[0] if counts else "",
        "mean_alt_allele_frequency": f"{sum(dataset.allele_frequencies) / len(dataset.allele_frequencies):.8f}",
        "first_variant": f"{first_variant.chrom}:{first_variant.position}:{first_variant.ref}>{first_variant.alt}",
        "last_variant": f"{last_variant.chrom}:{last_variant.position}:{last_variant.ref}>{last_variant.alt}",
        "rss_before_bytes": memory.rss_before_bytes,
        "rss_after_bytes": memory.rss_after_bytes,
        "rss_delta_bytes": memory.rss_delta_bytes,
        "peak_rss_bytes": memory.peak_rss_bytes,
        "python_peak_bytes": memory.python_peak_bytes,
    }


def main() -> int:
    output = ROOT / "results" / "1000g_chr22_public_matrix.csv"
    output.parent.mkdir(parents=True, exist_ok=True)

    rows = []
    for index, experiment in enumerate(EXPERIMENTS):
        print(
            "Running",
            experiment["name"],
            f"samples={experiment['samples']}",
            f"snps={experiment['snps']}",
            f"backend={experiment['backend']}",
            flush=True,
        )
        rows.append(
            run_one(
                seed=20260731 + index,
                min_minor_allele_count=1,
                **experiment,
            )
        )

    with output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    print(f"Wrote {output}")
    for row in rows:
        print(
            row["experiment"],
            f"status={row['status']}",
            f"load={row['load_seconds']}",
            f"plain={row['plain_prs_seconds']}",
            f"eval={row['evaluation_seconds']}",
            f"error={row['max_abs_error']}",
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
