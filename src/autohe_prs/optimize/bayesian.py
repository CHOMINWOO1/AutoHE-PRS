"""Small deterministic Gaussian-process refinement for discrete plan spaces."""

from __future__ import annotations

from dataclasses import asdict, dataclass
import math
from typing import Callable

from autohe_prs.config import FHEPlan


def _vector(plan: FHEPlan) -> tuple[float, ...]:
    scheme = {"CKKS": 0.0, "BFV": 1.0, "BGV": 2.0}[plan.scheme.upper()]
    return (
        scheme,
        math.log2(plan.ring_dimension),
        math.log2(plan.window_size),
        float(plan.scale_bits or 0) / 60.0,
        math.log10(plan.fixed_point_scaling_factor or 1),
    )


def _kernel(left: tuple[float, ...], right: tuple[float, ...], length: float = 1.5) -> float:
    distance = sum((a - b) ** 2 for a, b in zip(left, right))
    return math.exp(-distance / (2 * length * length))


def _solve(matrix: list[list[float]], values: list[float]) -> list[float]:
    size = len(values)
    work = [row[:] + [value] for row, value in zip(matrix, values)]
    for col in range(size):
        pivot = max(range(col, size), key=lambda row: abs(work[row][col]))
        work[col], work[pivot] = work[pivot], work[col]
        divisor = work[col][col] or 1e-12
        work[col] = [value / divisor for value in work[col]]
        for row in range(size):
            if row == col:
                continue
            factor = work[row][col]
            work[row] = [a - factor * b for a, b in zip(work[row], work[col])]
    return [work[row][-1] for row in range(size)]


@dataclass(frozen=True)
class SurrogatePrediction:
    plan: FHEPlan
    mean: float
    standard_deviation: float
    expected_improvement: float


def suggest(
    candidates: list[FHEPlan],
    observations: list[tuple[FHEPlan, float]],
    *,
    noise: float = 1e-6,
) -> list[SurrogatePrediction]:
    if not observations:
        return [
            SurrogatePrediction(plan, 0.0, 1.0, 1.0)
            for plan in sorted(candidates, key=repr)
        ]
    observed_vectors = [_vector(plan) for plan, _ in observations]
    observed_values = [value for _, value in observations]
    gram = [
        [_kernel(left, right) + (noise if i == j else 0.0)
         for j, right in enumerate(observed_vectors)]
        for i, left in enumerate(observed_vectors)
    ]
    alpha = _solve(gram, observed_values)
    best = min(observed_values)
    result: list[SurrogatePrediction] = []
    observed_repr = {repr(plan) for plan, _ in observations}
    for plan in candidates:
        if repr(plan) in observed_repr:
            continue
        vector = _vector(plan)
        covariance = [_kernel(vector, known) for known in observed_vectors]
        mean = sum(value * weight for value, weight in zip(covariance, alpha))
        projection = _solve(gram, covariance)
        variance = max(1e-12, 1.0 - sum(a * b for a, b in zip(covariance, projection)))
        sigma = math.sqrt(variance)
        improvement = best - mean
        z = improvement / sigma
        cdf = 0.5 * (1.0 + math.erf(z / math.sqrt(2.0)))
        pdf = math.exp(-0.5 * z * z) / math.sqrt(2.0 * math.pi)
        expected = max(0.0, improvement * cdf + sigma * pdf)
        result.append(SurrogatePrediction(plan, mean, sigma, expected))
    return sorted(result, key=lambda item: (-item.expected_improvement, item.mean, repr(item.plan)))


def bayesian_refine(
    candidates: list[FHEPlan],
    initial_observations: list[tuple[FHEPlan, float]],
    evaluator: Callable[[FHEPlan], float | None],
    *,
    iterations: int = 5,
) -> list[tuple[FHEPlan, float]]:
    observations = list(initial_observations)
    for _ in range(iterations):
        suggestions = suggest(candidates, observations)
        if not suggestions:
            break
        plan = suggestions[0].plan
        measured = evaluator(plan)
        if measured is not None and math.isfinite(measured):
            observations.append((plan, measured))
        else:
            # Failed observations receive a finite penalty and remain auditable
            # through the caller's benchmark record.
            worst = max((value for _, value in observations), default=1.0)
            observations.append((plan, worst * 10.0))
    return observations
