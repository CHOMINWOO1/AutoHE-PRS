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

from genome_he_pilot import prs_scores, summarize_scores  # noqa: E402
from genome_he_pilot.pgs import read_pgs_scoring_file  # noqa: E402
from genome_he_pilot.pgs_vcf import load_pgs_matched_vcf_dataset  # noqa: E402
from genome_he_pilot.vcf_cache import load_pgs_matched_cached_dataset  # noqa: E402
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
DEFAULT_PGS = ROOT / "data" / "pgs_catalog" / "PGS000348" / "PGS000348.txt.gz"
DEFAULT_CACHE = ROOT / "data" / "cache" / "1000g_chr22_samples32.sqlite"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--pgs", type=Path, default=DEFAULT_PGS)
    parser.add_argument("--vcf", type=Path, default=DEFAULT_VCF)
    parser.add_argument("--chrom", default="22")
    parser.add_argument("--samples", type=int, default=32)
    parser.add_argument("--cache", type=Path, default=DEFAULT_CACHE)
    parser.add_argument("--no-cache", action="store_true")
    parser.add_argument(
        "--max-snps",
        type=int,
        help="Limit matched SNPs for pilot encrypted runs.",
    )
    parser.add_argument("--backend", choices=["plaintext", "tenseal"], default="plaintext")
    parser.add_argument(
        "--output",
        type=Path,
        default=ROOT / "results" / "real_pgs_chr22_summary.csv",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    args.output.parent.mkdir(parents=True, exist_ok=True)

    start = perf_counter()
    scoring = read_pgs_scoring_file(args.pgs, chrom=args.chrom)
    read_pgs_seconds = perf_counter() - start

    start = perf_counter()
    if args.cache.exists() and not args.no_cache:
        match = load_pgs_matched_cached_dataset(
            cache_path=args.cache,
            scoring=scoring,
            chrom=args.chrom,
            max_snps=args.max_snps,
        )
        match_mode = "sqlite_cache"
    else:
        match = load_pgs_matched_vcf_dataset(
            vcf_path=args.vcf,
            scoring=scoring,
            chrom=args.chrom,
            n_samples=args.samples,
            max_snps=args.max_snps,
        )
        match_mode = "sequential_vcf_scan"
    match_seconds = perf_counter() - start
    dataset = match.dataset

    start = perf_counter()
    plain_scores = prs_scores(dataset.genotypes, dataset.weights)
    plain_prs_seconds = perf_counter() - start
    summary = summarize_scores(plain_scores)

    status = "plaintext_only"
    evaluation_seconds = ""
    encryption_seconds = ""
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
            status = "tenseal_ok"
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

    row = {
        "pgs_id": scoring.metadata.get("pgs_id", ""),
        "trait_reported": scoring.metadata.get("trait_reported", ""),
        "trait_mapped": scoring.metadata.get("trait_mapped", ""),
        "genome_build": scoring.metadata.get("genome_build", ""),
        "chrom": args.chrom,
        "samples": dataset.n_samples,
        "matched_snps": dataset.n_snps,
        "max_snps": args.max_snps or "",
        "requested_weights_on_chrom": match.requested_weights_on_chrom,
        "skipped_allele_mismatch": match.skipped_allele_mismatch,
        "backend": args.backend,
        "status": status,
        "match_mode": match_mode,
        "read_pgs_seconds": f"{read_pgs_seconds:.6f}",
        "match_vcf_seconds": f"{match_seconds:.6f}",
        "plain_prs_seconds": f"{plain_prs_seconds:.6f}",
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
        "first_variant": (
            f"{dataset.variants[0].chrom}:{dataset.variants[0].position}:"
            f"{dataset.variants[0].ref}>{dataset.variants[0].alt}"
        ),
    }

    with args.output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(row))
        writer.writeheader()
        writer.writerow(row)

    print(f"Wrote {args.output}")
    print(
        "Real PGS chr22:",
        f"pgs={row['pgs_id']}",
        f"trait={row['trait_reported']}",
        f"samples={dataset.n_samples}",
        f"matched_snps={dataset.n_snps}",
        f"match_mode={match_mode}",
        f"status={status}",
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
