"""Reproducible synthetic PRS workload matrix."""

from __future__ import annotations

from dataclasses import dataclass
import math
import random


@dataclass(frozen=True)
class SyntheticWorkload:
    dosages: list[list[float | None]]
    weights: list[float]
    genotype_encoding: str
    missing_rate: float
    seed: int


def make_synthetic_workload(
    *,
    samples: int,
    variants: int,
    genotype_encoding: str = "hard_call",
    weight_standard_deviation: float = 0.05,
    weight_decimal_places: int = 7,
    sparsity: float = 0.0,
    missing_rate: float = 0.0,
    seed: int = 20260728,
) -> SyntheticWorkload:
    if samples <= 0 or variants <= 0:
        raise ValueError("samples and variants must be positive")
    if genotype_encoding not in {"hard_call", "imputed"}:
        raise ValueError("genotype encoding must be hard_call or imputed")
    if not 0.0 <= missing_rate < 1.0 or not 0.0 <= sparsity < 1.0:
        raise ValueError("missing rate and sparsity must be in [0,1)")
    rng = random.Random(seed)
    allele_frequencies = [rng.uniform(0.01, 0.5) for _ in range(variants)]
    dosages: list[list[float | None]] = []
    for _ in range(samples):
        row: list[float | None] = []
        for frequency in allele_frequencies:
            hard = float(
                int(rng.random() < frequency) + int(rng.random() < frequency)
            )
            value = (
                hard
                if genotype_encoding == "hard_call"
                else min(2.0, max(0.0, hard + rng.gauss(0.0, 0.08)))
            )
            row.append(None if rng.random() < missing_rate else value)
        dosages.append(row)
    weights = [
        0.0 if rng.random() < sparsity else round(
            rng.gauss(0.0, weight_standard_deviation), weight_decimal_places
        )
        for _ in range(variants)
    ]
    return SyntheticWorkload(dosages, weights, genotype_encoding, missing_rate, seed)


def apply_missing_policy(
    workload: SyntheticWorkload, policy: str
) -> tuple[list[list[float]], list[float]]:
    if policy not in {"error", "exclude_variant", "zero", "mean_dosage"}:
        raise ValueError(f"unsupported missing policy: {policy}")
    missing_columns = {
        column
        for column in range(len(workload.weights))
        if any(row[column] is None for row in workload.dosages)
    }
    if missing_columns and policy == "error":
        raise ValueError("synthetic workload contains missing dosage")
    if policy == "exclude_variant":
        keep = [index for index in range(len(workload.weights)) if index not in missing_columns]
        return (
            [[float(row[index]) for index in keep] for row in workload.dosages],
            [workload.weights[index] for index in keep],
        )
    means = []
    for column in range(len(workload.weights)):
        observed = [
            float(row[column]) for row in workload.dosages if row[column] is not None
        ]
        means.append(math.fsum(observed) / len(observed) if observed else 0.0)
    return (
        [
            [
                float(value) if value is not None
                else (0.0 if policy == "zero" else means[column])
                for column, value in enumerate(row)
            ]
            for row in workload.dosages
        ],
        list(workload.weights),
    )
