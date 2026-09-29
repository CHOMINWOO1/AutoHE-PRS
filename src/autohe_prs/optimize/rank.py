"""Constraint filtering and multi-objective ranking."""

from __future__ import annotations

from dataclasses import asdict, dataclass
import math
from typing import Any

from autohe_prs.config import Constraints, FHEPlan
from autohe_prs.fhe.validator import validate_plan
from autohe_prs.fhe.packing import operation_counts
from autohe_prs.models.train import CostModelBundle


@dataclass(frozen=True)
class RankedCandidate:
    rank: int
    score: float
    plan: FHEPlan
    predictions: dict[str, float]
    prediction_target_scopes: dict[str, str]
    feasible: bool
    rejection_reasons: tuple[str, ...]


OBJECTIVE_TARGET = {
    "evaluation_latency": "evaluation_seconds",
    "total_latency": "total_seconds",
    "memory": "peak_memory_mb",
    "ciphertext_size": "ciphertext_bytes",
    "key_size": "key_bytes",
}


def rank_candidates(
    candidates: list[FHEPlan],
    model: CostModelBundle,
    workload_features: dict[str, Any],
    constraints: Constraints,
    *,
    objective: str = "evaluation_latency",
    min_success_probability: float = 0.5,
) -> list[RankedCandidate]:
    evaluated: list[tuple[float, FHEPlan, dict[str, float], tuple[str, ...]]] = []
    prediction_target_scopes = dict(getattr(model, "target_scopes", {}))
    if (
        "key_bytes" in getattr(model, "regressors", {})
        and "key_bytes" not in prediction_target_scopes
    ):
        prediction_target_scopes["key_bytes"] = "legacy_context+public_key_proxy"
    for plan in candidates:
        validation = validate_plan(
            plan,
            samples=int(workload_features["sample_count"]),
            variants=int(workload_features["matched_variant_count"]),
            genotype_min=float(workload_features.get("genotype_min", 0.0)),
            genotype_max=float(workload_features.get("genotype_max", 2.0)),
            sum_abs_weights=float(workload_features["sum_absolute_weights"]),
            constraints=constraints,
            imputed_dosage=bool(workload_features.get("imputed_dosage", False)),
        )
        reasons = [issue.failure_type for issue in validation.issues]
        if reasons:
            evaluated.append((math.inf, plan, {}, tuple(reasons)))
            continue
        counts = operation_counts(
            samples=int(workload_features["sample_count"]),
            variants=int(workload_features["matched_variant_count"]),
            variants_per_ciphertext=plan.variants_per_ciphertext,
            aggregate_strategy=plan.chunk_aggregation_strategy,
        )
        row = {
            **workload_features,
            **asdict(plan),
            **counts,
            "scheme": plan.scheme,
            "slots_used": plan.window_size,
            "fixed_point_scale": plan.fixed_point_scaling_factor,
            "ckks_scale_bits": plan.scale_bits,
        }
        predictions = model.predict(row)
        if (constraints.max_memory_mb is not None
                and predictions.get("peak_memory_mb", 0.0) > constraints.max_memory_mb):
            reasons.append("predicted_memory_exceeded")
        if (constraints.max_latency_seconds is not None
                and predictions.get("total_seconds", 0.0) > constraints.max_latency_seconds):
            reasons.append("predicted_latency_exceeded")
        if (constraints.max_absolute_error is not None
                and predictions.get("max_absolute_error", 0.0) > constraints.max_absolute_error):
            reasons.append("predicted_absolute_error_exceeded")
        if (constraints.max_relative_error is not None
                and predictions.get("max_relative_error", 0.0) > constraints.max_relative_error):
            reasons.append("predicted_relative_error_exceeded")
        if (
            "execution_success_probability" in predictions
            and predictions["execution_success_probability"] < min_success_probability
        ):
            reasons.append("predicted_execution_failure")
        if (
            "error_constraint_feasibility_probability" in predictions
            and predictions["error_constraint_feasibility_probability"]
            < min_success_probability
        ):
            reasons.append("predicted_error_constraint_failure")
        target = OBJECTIVE_TARGET.get(objective)
        if objective == "balanced_multi_objective":
            values = [
                math.log1p(predictions.get(name, 0.0))
                for name in ("evaluation_seconds", "peak_memory_mb", "ciphertext_bytes",
                             "key_bytes", "max_absolute_error")
            ]
            score = sum(values)
        elif target is None:
            raise ValueError(f"unsupported objective: {objective}")
        else:
            score = predictions.get(target, math.inf)
        evaluated.append((score, plan, predictions, tuple(reasons)))
    evaluated.sort(key=lambda item: (bool(item[3]), item[0], repr(item[1])))
    return [
        RankedCandidate(
            index + 1,
            score,
            plan,
            predictions,
            prediction_target_scopes,
            not reasons,
            reasons,
        )
        for index, (score, plan, predictions, reasons) in enumerate(evaluated)
    ]
