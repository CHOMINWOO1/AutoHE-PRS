"""End-to-end learned ranking plus top-K measured validation."""

from __future__ import annotations

from dataclasses import asdict, dataclass
import math
from typing import Any, Callable

from autohe_prs.config import Constraints, FHEPlan
from autohe_prs.fhe.runner import ExecutionResult, run_openfhe_plan
from autohe_prs.models.train import CostModelBundle
from .pareto import pareto_frontier
from .rank import RankedCandidate, rank_candidates
from .bayesian import bayesian_refine


@dataclass(frozen=True)
class TuningResult:
    selection_status: str
    selected_plan: FHEPlan | None
    ranking: tuple[RankedCandidate, ...]
    measured_records: tuple[dict[str, Any], ...]
    pareto_records: tuple[dict[str, Any], ...]
    selected_record: dict[str, Any] | None


MEASURED_OBJECTIVES = {
    "evaluation_latency": "evaluation_seconds",
    "total_latency": "total_seconds",
    "memory": "peak_memory_mb",
    "ciphertext_size": "ciphertext_bytes",
    "key_size": "key_bytes",
}


def _annotate_measurement(
    result: ExecutionResult | dict[str, Any],
    candidate: RankedCandidate | None,
    search_stage: str,
) -> dict[str, Any]:
    record = dict(result) if isinstance(result, dict) else asdict(result)
    predictions = candidate.predictions if candidate else {}
    prediction_scopes = candidate.prediction_target_scopes if candidate else {}
    record["search_stage"] = search_stage
    record["predicted_targets"] = dict(predictions)
    record["prediction_target_scopes"] = dict(prediction_scopes)
    if candidate is not None:
        record["predicted_rank"] = candidate.rank
    if record.get("status") == "ok":
        record["prediction_residuals"] = {
            key: float(record[key]) - float(predicted)
            for key, predicted in predictions.items()
            if key in record
            and record.get(key) is not None
            and not key.endswith("_probability")
            and (
                key != "key_bytes"
                or not prediction_scopes.get(key)
                or prediction_scopes[key] == record.get("key_bytes_scope")
            )
        }
        record["prediction_comparison_status"] = "paired"
    else:
        record["prediction_residuals"] = {}
        record["prediction_comparison_status"] = "unavailable_failed_measurement"
    return record


def _attach_measured_objective(
    records: list[dict[str, Any]], objective: str
) -> None:
    if objective != "balanced_multi_objective":
        key = MEASURED_OBJECTIVES.get(objective)
        if key is None:
            raise ValueError(f"unsupported objective: {objective}")
        for record in records:
            value = record.get(key)
            record["measured_objective_score"] = (
                float(value) if value is not None else None
            )
        return
    targets = (
        "evaluation_seconds",
        "peak_memory_mb",
        "ciphertext_bytes",
        "key_bytes",
        "max_absolute_error",
    )
    available = [
        record for record in records
        if all(record.get(target) is not None for target in targets)
    ]
    ranges: dict[str, tuple[float, float]] = {}
    for target in targets:
        values = [math.log1p(max(0.0, float(record[target]))) for record in available]
        if values:
            ranges[target] = (min(values), max(values))
    for record in records:
        if record not in available:
            record["measured_objective_score"] = None
            continue
        normalized = []
        for target, (low, high) in ranges.items():
            value = math.log1p(max(0.0, float(record[target])))
            normalized.append((value - low) / (high - low) if high > low else 0.0)
        record["measured_objective_score"] = sum(normalized)


def _raw_objective(record: dict[str, Any], objective: str) -> float | None:
    if objective == "balanced_multi_objective":
        targets = (
            "evaluation_seconds",
            "peak_memory_mb",
            "ciphertext_bytes",
            "key_bytes",
            "max_absolute_error",
        )
        if any(record.get(target) is None for target in targets):
            return None
        return sum(
            math.log1p(max(0.0, float(record[target])))
            for target in targets
        )
    key = MEASURED_OBJECTIVES.get(objective)
    if key is None:
        raise ValueError(f"unsupported objective: {objective}")
    value = record.get(key)
    return float(value) if value is not None else None


def tune(
    candidates: list[FHEPlan],
    model: CostModelBundle,
    workload_features: dict[str, Any],
    constraints: Constraints,
    dosages: list[list[float]],
    weights: list[float],
    *,
    objective: str = "evaluation_latency",
    top_k: int = 10,
    evaluator: Callable[[FHEPlan], ExecutionResult | dict[str, Any]] | None = None,
    use_bayesian_refinement: bool = False,
    bayesian_iterations: int = 5,
) -> TuningResult:
    if top_k <= 0:
        raise ValueError("top_k must be positive")
    ranking = rank_candidates(candidates, model, workload_features, constraints, objective=objective)
    feasible = [candidate for candidate in ranking if candidate.feasible][:top_k]
    evaluator = evaluator or (
        lambda plan: run_openfhe_plan(plan, dosages, weights, constraints=constraints)
    )
    measured: list[dict[str, Any]] = []
    candidates_by_plan = {repr(item.plan): item for item in ranking}
    for candidate in feasible:
        result = evaluator(candidate.plan)
        record = _annotate_measurement(result, candidate, "top_k_actual_validation")
        measured.append(record)
    if use_bayesian_refinement:
        observations = [
            (FHEPlan(**record["plan"]), float(value))
            for record in measured
            if record.get("status") == "ok"
            and (value := _raw_objective(record, objective)) is not None
        ]
        if observations:
            def refinement_evaluator(plan: FHEPlan) -> float | None:
                result = evaluator(plan)
                record = _annotate_measurement(
                    result,
                    candidates_by_plan.get(repr(plan)),
                    "bayesian_refinement",
                )
                measured.append(record)
                value = _raw_objective(record, objective)
                return (
                    float(value)
                    if record.get("status") == "ok" and value is not None
                    else None
                )

            bayesian_refine(
                [candidate.plan for candidate in ranking if candidate.feasible],
                observations,
                refinement_evaluator,
                iterations=bayesian_iterations,
            )
    _attach_measured_objective(measured, objective)
    successful = [
        record for record in measured
        if record.get("status") == "ok"
        and record.get("measured_objective_score") is not None
    ]
    if successful:
        successful.sort(
            key=lambda row: float(row["measured_objective_score"])
        )
        selected_record = successful[0]
        selected_raw = selected_record["plan"]
        selected = FHEPlan(**selected_raw)
        status = "measured_and_validated"
    else:
        selected_record = None
        selected = feasible[0].plan if feasible else None
        status = "predicted_only_no_successful_actual_validation" if selected else "no_feasible_candidate"
    frontier = pareto_frontier(successful)
    return TuningResult(
        status,
        selected,
        tuple(ranking),
        tuple(measured),
        tuple(frontier),
        selected_record,
    )
