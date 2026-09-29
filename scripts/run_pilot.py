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
    make_synthetic_dataset,
    prs_scores,
    summarize_scores,
)
from genome_he_pilot.tenseal_backend import (  # noqa: E402
    TenSEALUnavailable,
    encrypted_prs_scores,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--samples", type=int, default=128)
    parser.add_argument("--snps", type=int, default=2048)
    parser.add_argument("--seed", type=int, default=20260701)
    parser.add_argument(
        "--backend",
        choices=["plaintext", "tenseal"],
        default="plaintext",
        help="Run plaintext baseline or optional TenSEAL encrypted PRS.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=ROOT / "results" / "pilot_summary.csv",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    args.output.parent.mkdir(parents=True, exist_ok=True)

    start = perf_counter()
    dataset = make_synthetic_dataset(
        n_samples=args.samples,
        n_snps=args.snps,
        seed=args.seed,
    )
    data_seconds = perf_counter() - start

    start = perf_counter()
    plain_scores = prs_scores(dataset.genotypes, dataset.weights)
    plain_seconds = perf_counter() - start
    summary = summarize_scores(plain_scores)

    start = perf_counter()
    counts = allele_counts(dataset.genotypes)
    control_counts, case_counts = case_control_allele_counts(
        dataset.genotypes, dataset.phenotypes
    )
    count_seconds = perf_counter() - start

    encrypted_seconds = ""
    max_abs_error = ""
    backend_status = "plaintext_only"

    if args.backend == "tenseal":
        try:
            start = perf_counter()
            encrypted_scores = encrypted_prs_scores(dataset.genotypes, dataset.weights)
            encrypted_seconds = perf_counter() - start
            max_abs_error = max(
                abs(a - b) for a, b in zip(plain_scores, encrypted_scores)
            )
            backend_status = "tenseal_ok"
        except TenSEALUnavailable as exc:
            backend_status = f"tenseal_unavailable: {exc}"

    row = {
        "samples": dataset.n_samples,
        "snps": dataset.n_snps,
        "backend": args.backend,
        "backend_status": backend_status,
        "data_seconds": f"{data_seconds:.6f}",
        "plain_prs_seconds": f"{plain_seconds:.6f}",
        "count_seconds": f"{count_seconds:.6f}",
        "encrypted_prs_seconds": encrypted_seconds,
        "max_abs_error": max_abs_error,
        "prs_mean": f"{summary.mean:.8f}",
        "prs_sd": f"{summary.standard_deviation:.8f}",
        "prs_min": f"{summary.minimum:.8f}",
        "prs_max": f"{summary.maximum:.8f}",
        "first_snp_allele_count": counts[0] if counts else "",
        "first_snp_control_count": control_counts[0] if control_counts else "",
        "first_snp_case_count": case_counts[0] if case_counts else "",
    }

    write_header = not args.output.exists()
    with args.output.open("a", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(row))
        if write_header:
            writer.writeheader()
        writer.writerow(row)

    print(f"Wrote {args.output}")
    print(
        "Baseline:",
        f"samples={dataset.n_samples}",
        f"snps={dataset.n_snps}",
        f"prs_seconds={plain_seconds:.6f}",
        f"count_seconds={count_seconds:.6f}",
    )
    if backend_status != "plaintext_only":
        print(f"Backend status: {backend_status}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

