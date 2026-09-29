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
    case_control_allele_counts,
    make_synthetic_dataset,
    prs_scores,
    summarize_scores,
)


EXPERIMENTS = [
    {"name": "E1_validation", "samples": 32, "snps": 128},
    {"name": "E2_small_pilot", "samples": 128, "snps": 2048},
    {"name": "E3_kci_main", "samples": 512, "snps": 8192},
    {"name": "E4_stress", "samples": 2048, "snps": 32768},
]


def run_one(name: str, samples: int, snps: int, seed: int) -> dict[str, str | int]:
    start = perf_counter()
    dataset = make_synthetic_dataset(n_samples=samples, n_snps=snps, seed=seed)
    data_seconds = perf_counter() - start

    start = perf_counter()
    scores = prs_scores(dataset.genotypes, dataset.weights)
    prs_seconds = perf_counter() - start
    summary = summarize_scores(scores)

    start = perf_counter()
    counts = allele_counts(dataset.genotypes)
    control_counts, case_counts = case_control_allele_counts(
        dataset.genotypes, dataset.phenotypes
    )
    count_seconds = perf_counter() - start

    return {
        "experiment": name,
        "samples": samples,
        "snps": snps,
        "seed": seed,
        "data_seconds": f"{data_seconds:.6f}",
        "plain_prs_seconds": f"{prs_seconds:.6f}",
        "count_seconds": f"{count_seconds:.6f}",
        "prs_mean": f"{summary.mean:.8f}",
        "prs_sd": f"{summary.standard_deviation:.8f}",
        "prs_min": f"{summary.minimum:.8f}",
        "prs_max": f"{summary.maximum:.8f}",
        "first_snp_allele_count": counts[0] if counts else "",
        "first_snp_control_count": control_counts[0] if control_counts else "",
        "first_snp_case_count": case_counts[0] if case_counts else "",
    }


def main() -> int:
    output = ROOT / "results" / "kci_experiment_matrix.csv"
    output.parent.mkdir(parents=True, exist_ok=True)

    rows = []
    for index, experiment in enumerate(EXPERIMENTS):
        rows.append(run_one(seed=20260701 + index, **experiment))

    with output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    print(f"Wrote {output}")
    for row in rows:
        print(
            row["experiment"],
            f"samples={row['samples']}",
            f"snps={row['snps']}",
            f"prs_seconds={row['plain_prs_seconds']}",
            f"count_seconds={row['count_seconds']}",
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

