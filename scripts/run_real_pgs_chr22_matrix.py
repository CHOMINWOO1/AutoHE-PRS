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
from genome_he_pilot.tenseal_backend import (  # noqa: E402
    TenSEALUnavailable,
    encrypted_prs_benchmark,
)
from genome_he_pilot.vcf_cache import load_pgs_matched_cached_dataset  # noqa: E402


DEFAULT_VCF = (
    ROOT
    / "data"
    / "1000genomes"
    / "phase3"
    / "ALL.chr22.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"
)
DEFAULT_PGS = ROOT / "data" / "pgs_catalog" / "PGS004941" / "PGS004941.txt.gz"
DEFAULT_CACHE = ROOT / "data" / "cache" / "1000g_chr22_samples32.sqlite"
DEFAULT_OUTPUT = ROOT / "results" / "real_pgs_chr22_pgs004941_ckks_matrix.csv"


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
        nargs="+",
        default=[512, 1024, 2048, 4096],
        help="Matched SNP subset sizes for CKKS pilot scaling.",
    )
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    return parser.parse_args()


def run_row(args: argparse.Namespace, scoring, max_snps: int) -> dict[str, str]:
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

    start = perf_counter()
    plain_scores = prs_scores(dataset.genotypes, dataset.weights)
    plain_prs_seconds = perf_counter() - start
    summary = summarize_scores(plain_scores)

    status = "tenseal_ok"
    encryption_seconds = ""
    evaluation_seconds = ""
    decryption_seconds = ""
    max_abs_error = ""

    try:
        encrypted = encrypted_prs_benchmark(dataset.genotypes, dataset.weights)
        max_abs_error_value = max(
            abs(plain - secure)
            for plain, secure in zip(plain_scores, encrypted.scores)
        )
        max_abs_error = f"{max_abs_error_value:.12g}"
        encryption_seconds = f"{encrypted.encryption_seconds:.6f}"
        evaluation_seconds = f"{encrypted.evaluation_seconds:.6f}"
        decryption_seconds = f"{encrypted.decryption_seconds:.6f}"
        input_ciphertext_mean_bytes = f"{encrypted.input_ciphertext_bytes / dataset.n_samples:.1f}"
        result_ciphertext_mean_bytes = f"{encrypted.result_ciphertext_bytes / dataset.n_samples:.1f}"
        public_context_bytes = str(encrypted.public_context_bytes)
        poly_modulus_degree = str(encrypted.poly_modulus_degree)
        scale = str(encrypted.scale)
        coeff_mod_bit_sizes = "-".join(str(value) for value in encrypted.coeff_mod_bit_sizes)
    except TenSEALUnavailable as exc:
        status = f"tenseal_unavailable: {exc}"
        input_ciphertext_mean_bytes = ""
        result_ciphertext_mean_bytes = ""
        public_context_bytes = ""
        poly_modulus_degree = ""
        scale = ""
        coeff_mod_bit_sizes = ""

    return {
        "pgs_id": scoring.metadata.get("pgs_id", ""),
        "trait_reported": scoring.metadata.get("trait_reported", ""),
        "genome_build": scoring.metadata.get("genome_build", ""),
        "chrom": args.chrom,
        "samples": str(dataset.n_samples),
        "matched_snps": str(dataset.n_snps),
        "max_snps": str(max_snps),
        "requested_weights_on_chrom": str(match.requested_weights_on_chrom),
        "skipped_allele_mismatch": str(match.skipped_allele_mismatch),
        "backend": "tenseal",
        "status": status,
        "match_mode": match_mode,
        "match_vcf_seconds": f"{match_seconds:.6f}",
        "plain_prs_seconds": f"{plain_prs_seconds:.6f}",
        "encryption_seconds": encryption_seconds,
        "evaluation_seconds": evaluation_seconds,
        "decryption_seconds": decryption_seconds,
        "max_abs_error": max_abs_error,
        "prs_mean": f"{summary.mean:.8f}",
        "prs_sd": f"{summary.standard_deviation:.8f}",
        "prs_min": f"{summary.minimum:.8f}",
        "prs_max": f"{summary.maximum:.8f}",
        "input_ciphertext_mean_bytes": input_ciphertext_mean_bytes,
        "result_ciphertext_mean_bytes": result_ciphertext_mean_bytes,
        "public_context_bytes": public_context_bytes,
        "poly_modulus_degree": poly_modulus_degree,
        "scale": scale,
        "coeff_mod_bit_sizes": coeff_mod_bit_sizes,
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
            "CAD PGS CKKS matrix:",
            f"max_snps={max_snps}",
            f"matched={row['matched_snps']}",
            f"eval={row['evaluation_seconds']}",
            f"error={row['max_abs_error']}",
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
