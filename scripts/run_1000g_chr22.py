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

from genome_he_pilot import (  # noqa: E402
    allele_counts,
    case_control_allele_counts,
    load_vcf_dosage_subset,
    prs_scores,
    select_panel_samples,
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


DEFAULT_VCF = (
    ROOT
    / "data"
    / "1000genomes"
    / "phase3"
    / "ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"
)
DEFAULT_PANEL = (
    ROOT
    / "data"
    / "1000genomes"
    / "phase3"
    / "integrated_call_samples_v3.20130502.ALL.panel"
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--samples", type=int, default=32)
    parser.add_argument("--snps", type=int, default=512)
    parser.add_argument("--seed", type=int, default=20260721)
    parser.add_argument(
        "--min-mac",
        type=int,
        default=1,
        help="Minimum minor allele count within selected samples.",
    )
    parser.add_argument("--vcf", type=Path, default=DEFAULT_VCF)
    parser.add_argument("--panel", type=Path, default=DEFAULT_PANEL)
    parser.add_argument("--population")
    parser.add_argument("--super-population")
    parser.add_argument(
        "--backend",
        choices=["plaintext", "tenseal"],
        default="plaintext",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=ROOT / "results" / "1000g_chr22_summary.csv",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    memory_before = process_memory()

    requested_sample_ids = None
    if args.population or args.super_population:
        requested_sample_ids = select_panel_samples(
            args.panel,
            limit=args.samples,
            population=args.population,
            super_population=args.super_population,
        )

    start = perf_counter()
    dataset = load_vcf_dosage_subset(
        args.vcf,
        n_samples=args.samples,
        n_snps=args.snps,
        sample_ids=requested_sample_ids,
        seed=args.seed,
        min_minor_allele_count=args.min_mac,
    )
    load_seconds = perf_counter() - start

    start = perf_counter()
    plain_scores = prs_scores(dataset.genotypes, dataset.weights)
    plain_prs_seconds = perf_counter() - start
    summary = summarize_scores(plain_scores)

    start = perf_counter()
    counts = allele_counts(dataset.genotypes)
    control_counts, case_counts = case_control_allele_counts(
        dataset.genotypes,
        dataset.phenotypes,
    )
    count_seconds = perf_counter() - start

    backend_status = "plaintext_only"
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

    if args.backend == "tenseal":
        try:
            encrypted = encrypted_prs_benchmark(dataset.genotypes, dataset.weights)
            max_abs_error = max(
                abs(plain - secure)
                for plain, secure in zip(plain_scores, encrypted.scores)
            )
            backend_status = "tenseal_ok"
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
            backend_status = f"tenseal_unavailable: {exc}"

    memory = stop_process_memory_measurement(memory_before)

    first_variant = dataset.variants[0]
    last_variant = dataset.variants[-1]
    row = {
        "source": str(dataset.source_path),
        "samples": dataset.n_samples,
        "snps": dataset.n_snps,
        "backend": args.backend,
        "backend_status": backend_status,
        "population": args.population or "",
        "super_population": args.super_population or "",
        "min_minor_allele_count": args.min_mac,
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
        "first_snp_control_count": control_counts[0] if control_counts else "",
        "first_snp_case_count": case_counts[0] if case_counts else "",
        "first_variant": f"{first_variant.chrom}:{first_variant.position}:{first_variant.ref}>{first_variant.alt}",
        "last_variant": f"{last_variant.chrom}:{last_variant.position}:{last_variant.ref}>{last_variant.alt}",
        "first_sample": dataset.sample_ids[0],
        "last_sample": dataset.sample_ids[-1],
        "rss_before_bytes": memory.rss_before_bytes,
        "rss_after_bytes": memory.rss_after_bytes,
        "rss_delta_bytes": memory.rss_delta_bytes,
        "peak_rss_bytes": memory.peak_rss_bytes,
        "python_peak_bytes": memory.python_peak_bytes,
    }

    write_header = not args.output.exists()
    with args.output.open("a", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(row))
        if write_header:
            writer.writeheader()
        writer.writerow(row)

    print(f"Wrote {args.output}")
    print(
        "1000G chr22:",
        f"samples={dataset.n_samples}",
        f"snps={dataset.n_snps}",
        f"load_seconds={load_seconds:.6f}",
        f"prs_seconds={plain_prs_seconds:.6f}",
        f"count_seconds={count_seconds:.6f}",
    )
    print(
        "Variants:",
        row["first_variant"],
        "to",
        row["last_variant"],
    )
    if backend_status != "plaintext_only":
        print(f"Backend status: {backend_status}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
