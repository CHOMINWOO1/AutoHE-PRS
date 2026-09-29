"""Reproducible generalization split suite for AutoHE-PRS cost models."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

from .auto import train_auto_cost_models


GROUP_SPLITS = {
    "unseen_variant_count": "group_variant_count",
    "unseen_sample_count": "group_sample_count",
    "unseen_window_size": "group_window_size",
    "unseen_parameter_combination": "group_parameter_combination",
    "unseen_scheme_configuration": "group_scheme_configuration",
    "one_pgs_to_another": "group_pgs",
}


@dataclass(frozen=True)
class SplitEvaluation:
    name: str
    status: str
    group_key: str
    train_rows: int
    test_rows: int
    test_groups: tuple[str, ...]
    metrics: dict[str, dict[str, float]]
    reason: str | None = None


def _evaluate(
    rows: list[dict[str, Any]],
    *,
    name: str,
    group_key: str,
    backend: str,
    test_fraction: float,
    seed: int,
    held_out_groups: tuple[str, ...] | None = None,
) -> SplitEvaluation:
    groups = {
        str(row.get(group_key, ""))
        for row in rows
        if str(row.get(group_key, ""))
    }
    if held_out_groups is None and len(groups) < 2:
        return SplitEvaluation(
            name, "unavailable", group_key, 0, 0, (), {},
            "at least two non-empty groups are required",
        )
    try:
        model = train_auto_cost_models(
            rows,
            backend=backend,
            group_key=group_key,
            test_fraction=test_fraction,
            seed=seed,
            held_out_groups=held_out_groups,
        )
    except (ValueError, RuntimeError) as exc:
        return SplitEvaluation(
            name, "unavailable", group_key, 0, 0, (), {}, str(exc)
        )
    return SplitEvaluation(
        name=name,
        status="ok",
        group_key=group_key,
        train_rows=int(model.split.get("train_rows", 0)),
        test_rows=int(model.split.get("test_rows", 0)),
        test_groups=tuple(str(value) for value in model.split.get("test_groups", ())),
        metrics=model.metrics,
    )


def evaluate_generalization_splits(
    rows: list[dict[str, Any]],
    *,
    backend: str = "ridge",
    test_fraction: float = 0.2,
    seed: int = 20260728,
) -> dict[str, SplitEvaluation]:
    """Evaluate every mandated non-random generalization split.

    A split without sufficient source groups is retained as an auditable
    ``unavailable`` result rather than silently falling back to a random row
    split.
    """
    results = {
        name: _evaluate(
            rows,
            name=name,
            group_key=group_key,
            backend=backend,
            test_fraction=test_fraction,
            seed=seed,
        )
        for name, group_key in GROUP_SPLITS.items()
    }
    domains = {
        str(row.get("group_source_domain", ""))
        for row in rows
        if str(row.get("group_source_domain", "")) in {"synthetic", "real"}
    }
    if domains == {"synthetic", "real"}:
        results["synthetic_to_real"] = _evaluate(
            rows,
            name="synthetic_to_real",
            group_key="group_source_domain",
            backend=backend,
            test_fraction=test_fraction,
            seed=seed,
            held_out_groups=("real",),
        )
    else:
        results["synthetic_to_real"] = SplitEvaluation(
            "synthetic_to_real",
            "unavailable",
            "group_source_domain",
            0,
            0,
            (),
            {},
            "both synthetic and real source domains are required",
        )
    return results
