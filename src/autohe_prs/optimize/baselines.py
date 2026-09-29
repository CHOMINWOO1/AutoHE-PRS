"""Fixed, minimal, grid, random, and uninformed Bayesian search baselines."""

from __future__ import annotations

from dataclasses import dataclass
import math
import random
from typing import Callable

from autohe_prs.config import FHEPlan
from .bayesian import suggest


@dataclass(frozen=True)
class SearchObservation:
    iteration: int
    plan: FHEPlan
    objective_value: float | None
    success: bool


@dataclass(frozen=True)
class BaselineSearchResult:
    method: str
    observations: tuple[SearchObservation, ...]
    best_value: float | None
    best_plan: FHEPlan | None
    regret: float | None


def _finish(
    method: str,
    observations: list[SearchObservation],
    known_best: float | None,
) -> BaselineSearchResult:
    successful = [
        item for item in observations
        if item.success and item.objective_value is not None and math.isfinite(item.objective_value)
    ]
    best = min(successful, key=lambda item: float(item.objective_value)) if successful else None
    value = float(best.objective_value) if best else None
    regret = value - known_best if value is not None and known_best is not None else None
    return BaselineSearchResult(method, tuple(observations), value, best.plan if best else None, regret)


def run_baseline_search(
    candidates: list[FHEPlan],
    evaluator: Callable[[FHEPlan], float | None],
    *,
    method: str,
    budget: int,
    seed: int = 20260728,
    fixed_plan: FHEPlan | None = None,
    known_best: float | None = None,
) -> BaselineSearchResult:
    if budget <= 0:
        raise ValueError("search budget must be positive")
    if not candidates and fixed_plan is None:
        raise ValueError("search requires candidates")
    method = method.lower()
    if method == "fixed":
        order = [fixed_plan or candidates[0]]
    elif method == "minimal":
        order = [min(
            candidates,
            key=lambda plan: (
                plan.ring_dimension,
                plan.multiplicative_depth,
                plan.window_size,
                plan.scheme,
            ),
        )]
    elif method == "grid":
        order = sorted(candidates, key=repr)
    elif method == "random":
        order = list(candidates)
        random.Random(seed).shuffle(order)
    elif method == "bayesian":
        order = []
    else:
        raise ValueError(f"unsupported baseline search method: {method}")

    observations: list[SearchObservation] = []
    if method != "bayesian":
        for iteration, plan in enumerate(order[:budget], start=1):
            value = evaluator(plan)
            observations.append(
                SearchObservation(iteration, plan, value, value is not None and math.isfinite(value))
            )
        return _finish(method, observations, known_best)

    # Uninformed warm start: deterministic pseudo-random first observation.
    remaining = list(candidates)
    first = random.Random(seed).choice(remaining)
    for iteration in range(1, min(budget, len(candidates)) + 1):
        if iteration == 1:
            plan = first
        else:
            successful = [
                (item.plan, float(item.objective_value))
                for item in observations
                if item.success and item.objective_value is not None
            ]
            if successful:
                proposals = suggest(candidates, successful)
                used = {repr(item.plan) for item in observations}
                plan = next(item.plan for item in proposals if repr(item.plan) not in used)
            else:
                used = {repr(item.plan) for item in observations}
                plan = next(item for item in sorted(candidates, key=repr) if repr(item) not in used)
        value = evaluator(plan)
        observations.append(
            SearchObservation(iteration, plan, value, value is not None and math.isfinite(value))
        )
    return _finish(method, observations, known_best)
