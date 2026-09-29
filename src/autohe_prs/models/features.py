"""Stable numerical/categorical features for cost modeling."""

from __future__ import annotations

import math
from typing import Any


NUMERIC_FEATURES = (
    "matched_variant_count",
    "sample_count",
    "score_count",
    "total_pgs_variants",
    "ring_dimension",
    "multiplicative_depth",
    "slots_used",
    "window_size",
    "chunk_count",
    "ciphertext_count",
    "multiplication_count",
    "addition_count",
    "rotation_count",
    "key_count",
    "samples_per_ciphertext",
    "scores_per_ciphertext",
    "plaintext_modulus",
    "fixed_point_scale",
    "ckks_scale_bits",
    "scaling_modulus_size",
    "first_modulus_size",
    "genotype_min",
    "genotype_max",
    "missing_rate",
    "max_absolute_weight",
    "sum_absolute_weights",
    "weight_mean",
    "weight_standard_deviation",
    "weight_decimal_precision",
    "expected_score_bound",
    "score_sparsity",
    "score_overlap_ratio",
    "physical_core_count",
    "logical_core_count",
    "ram_bytes",
)
SCHEMES = ("CKKS", "BFV", "BGV")


def feature_names() -> tuple[str, ...]:
    return NUMERIC_FEATURES + tuple(f"scheme_{scheme}" for scheme in SCHEMES)


def feature_vector(
    row: dict[str, Any], names: tuple[str, ...] | None = None
) -> list[float]:
    values: list[float] = []
    for name in names or feature_names():
        if name.startswith("scheme_"):
            values.append(float(str(row.get("scheme", "")).upper() == name[7:]))
            continue
        raw = row.get(name)
        try:
            value = float(raw) if raw not in {None, ""} else 0.0
        except (TypeError, ValueError):
            value = 0.0
        if not math.isfinite(value):
            value = 0.0
        values.append(math.log1p(max(0.0, value)) if value >= 0 else -math.log1p(-value))
    return values
