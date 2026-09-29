"""Plaintext reference computation and numerical comparison."""

from __future__ import annotations

from dataclasses import dataclass
import math


@dataclass(frozen=True)
class ErrorMetrics:
    max_absolute_error: float
    mean_absolute_error: float
    max_relative_error: float
    mean_relative_error: float
    finite: bool


def prs(dosages: list[list[float]], weights: list[float]) -> list[float]:
    if any(len(row) != len(weights) for row in dosages):
        raise ValueError("every dosage row must match the weight count")
    return [math.fsum(dosage * weight for dosage, weight in zip(row, weights)) for row in dosages]


def prs_matrix(dosages: list[list[float]], score_weights: list[list[float]]) -> list[list[float]]:
    if not score_weights:
        return [[] for _ in dosages]
    variants = len(dosages[0]) if dosages else len(score_weights)
    if len(score_weights) != variants:
        raise ValueError("weight matrix rows must equal dosage columns")
    scores = len(score_weights[0])
    if any(len(row) != scores for row in score_weights):
        raise ValueError("weight matrix must be rectangular")
    if any(len(row) != variants for row in dosages):
        raise ValueError("dosage matrix must be rectangular")
    return [
        [
            math.fsum(dosages[sample][variant] * score_weights[variant][score]
                      for variant in range(variants))
            for score in range(scores)
        ]
        for sample in range(len(dosages))
    ]


def compare(reference: list[float], observed: list[float], epsilon: float = 1e-15) -> ErrorMetrics:
    if len(reference) != len(observed):
        raise ValueError("reference and observed lengths differ")
    absolute = [abs(left - right) for left, right in zip(reference, observed)]
    relative = [error / max(abs(left), epsilon) for left, error in zip(reference, absolute)]
    return ErrorMetrics(
        max(absolute, default=0.0),
        math.fsum(absolute) / len(absolute) if absolute else 0.0,
        max(relative, default=0.0),
        math.fsum(relative) / len(relative) if relative else 0.0,
        all(math.isfinite(value) for value in observed),
    )
