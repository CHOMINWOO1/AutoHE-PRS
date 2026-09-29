"""Named ablation switches used by reproducible study drivers."""

from __future__ import annotations

from dataclasses import dataclass, replace

from autohe_prs.config import Constraints, FHEPlan


@dataclass(frozen=True)
class AblationSettings:
    name: str
    select_scheme: bool = True
    optimize_packing: bool = True
    optimize_chunk_strategy: bool = True
    optimize_reduction_strategy: bool = True
    use_learned_cost_model: bool = True
    use_bayesian_refinement: bool = True
    use_genomic_features: bool = True
    use_top_k_validation: bool = True
    enforce_error_constraint: bool = True
    enforce_fixed_point_range: bool = True


ABLATIONS = {
    "no_scheme_selection": {"select_scheme": False},
    "no_packing_optimization": {"optimize_packing": False},
    "no_chunk_strategy_optimization": {"optimize_chunk_strategy": False},
    "no_reduction_strategy_optimization": {"optimize_reduction_strategy": False},
    "no_learned_cost_model": {"use_learned_cost_model": False},
    "no_bayesian_refinement": {"use_bayesian_refinement": False},
    "no_genomic_workload_features": {"use_genomic_features": False},
    "no_top_k_actual_validation": {"use_top_k_validation": False},
    "no_error_constraint": {"enforce_error_constraint": False},
    "no_fixed_point_range_validator": {"enforce_fixed_point_range": False},
}


def settings(name: str) -> AblationSettings:
    if name not in ABLATIONS:
        raise ValueError(f"unknown ablation: {name}")
    return replace(AblationSettings(name), **ABLATIONS[name])


def apply_plan_ablation(
    candidates: list[FHEPlan],
    experiment: AblationSettings,
    *,
    fixed_scheme: str = "CKKS",
    fixed_window: int = 4096,
) -> list[FHEPlan]:
    values = candidates
    if not experiment.select_scheme:
        values = [plan for plan in values if plan.scheme.upper() == fixed_scheme.upper()]
    if not experiment.optimize_packing:
        values = [plan for plan in values if plan.window_size == fixed_window]
    if not experiment.optimize_chunk_strategy:
        values = [
            plan for plan in values if plan.chunk_aggregation_strategy == "reduce_then_add"
        ]
    if not experiment.optimize_reduction_strategy:
        values = [plan for plan in values if plan.reduction_tree_strategy == "binary_tree"]
    return values


def apply_constraint_ablation(
    constraints: Constraints, experiment: AblationSettings
) -> Constraints:
    return replace(
        constraints,
        max_absolute_error=(
            constraints.max_absolute_error if experiment.enforce_error_constraint else None
        ),
        max_relative_error=(
            constraints.max_relative_error if experiment.enforce_error_constraint else None
        ),
    )
