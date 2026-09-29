"""Canonical feature rows from benchmark output."""

from __future__ import annotations

from typing import Any, Iterable

from autohe_prs.fhe.packing import operation_counts


def _first_present(*values: Any) -> Any:
    return next((value for value in values if value not in {None, ""}), None)


def _number(value: Any) -> float:
    return float(value) if value not in {None, ""} else 0.0


def canonicalize(record: dict[str, Any]) -> dict[str, Any]:
    plan = record.get("plan", {})
    if isinstance(plan, str):
        import json
        plan = json.loads(plan)
    samples = int(record.get("samples") or record.get("sample_count") or 0)
    variants = int(
        record.get("variants")
        or record.get("matched_snps")
        or record.get("snps")
        or record.get("matched_variant_count")
        or 0
    )
    window = int(
        plan.get("window_size")
        or record.get("window_size")
        or record.get("window_snps")
        or max(1, variants)
    )
    operations = operation_counts(
        samples=samples,
        variants=variants,
        variants_per_ciphertext=int(plan.get("variants_per_ciphertext") or window),
        aggregate_strategy=plan.get("chunk_aggregation_strategy", "reduce_then_add"),
    )
    backend = str(record.get("backend") or "")
    measurement_kind = record.get("measurement_kind")
    if not measurement_kind:
        if "openfhe" in backend.lower():
            measurement_kind = "measured_historical"
        elif record.get("status") is not None:
            measurement_kind = "measured_unspecified"
        else:
            measurement_kind = "unknown"
    total_seconds = _first_present(
        record.get("total_seconds"),
        record.get("total_he_compute_seconds"),
        record.get("total_ckks_compute_seconds"),
    )
    ciphertext_bytes = _first_present(record.get("ciphertext_bytes"))
    if ciphertext_bytes is None:
        input_bytes = _first_present(
            record.get("input_ciphertext_total_bytes"),
            record.get("input_ciphertext_bytes"),
        )
        result_bytes = _first_present(
            record.get("result_ciphertext_total_bytes"),
            record.get("result_ciphertext_bytes"),
        )
        if input_bytes is not None or result_bytes is not None:
            ciphertext_bytes = _number(input_bytes) + _number(result_bytes)
    key_bytes = _first_present(
        record.get("total_key_material_bytes"),
        record.get("key_bytes"),
        record.get("public_context_bytes"),
    )
    if record.get("key_bytes_scope"):
        key_bytes_scope = str(record["key_bytes_scope"])
    elif record.get("total_key_material_bytes") not in {None, ""}:
        key_bytes_scope = "public+secret+multiplication+rotation"
    elif record.get("key_bytes") not in {None, ""}:
        key_bytes_scope = "source_defined"
    elif record.get("public_context_bytes") not in {None, ""}:
        key_bytes_scope = "legacy_context+public_key_proxy"
    else:
        key_bytes_scope = "unavailable"
    max_absolute_error = _first_present(
        record.get("max_absolute_error"), record.get("max_abs_error")
    )
    error_feasible = record.get("error_constraint_feasible")
    if error_feasible in {None, ""}:
        if record.get("failure_type") == "error_constraint":
            error_feasible = 0
        elif record.get("status") == "ok" and max_absolute_error not in {None, ""}:
            # Successful runner records have already passed their attached
            # deterministic absolute/relative error constraints.
            error_feasible = 1
    peak_memory_mb = record.get("peak_memory_mb")
    if peak_memory_mb in {None, ""} and record.get("peak_rss_bytes") not in {None, ""}:
        peak_memory_mb = float(record["peak_rss_bytes"]) / (1024 * 1024)
    dataset_group = str(
        record.get("workload_group")
        or record.get("dataset_group")
        or record.get("pgs_id")
        or "unspecified"
    )
    source_text = " ".join(
        str(record.get(key, ""))
        for key in ("workload_type", "source", "dataset_group", "workload_group")
    ).lower()
    if "synthetic" in source_text:
        source_domain = "synthetic"
    elif any(token in source_text for token in ("1000g", "1000 genomes", "real", "pgs")):
        source_domain = "real"
    else:
        source_domain = "unspecified"
    scheme = str(record.get("scheme") or plan.get("scheme") or "").upper()
    ring_dimension = int(
        plan.get("ring_dimension")
        or record.get("ring_dimension")
        or record.get("ring_dim")
        or record.get("poly_modulus_degree")
        or 0
    )
    predictions = record.get("predicted_targets")
    prediction_scopes = record.get("prediction_target_scopes")
    prediction_scopes = (
        prediction_scopes if isinstance(prediction_scopes, dict) else {}
    )
    prediction_residuals: dict[str, float] = {}
    prediction_comparison_status = record.get("prediction_comparison_status")
    if isinstance(predictions, dict):
        if record.get("status") == "ok":
            for target, predicted in predictions.items():
                measured = record.get(target)
                if (
                    target.endswith("_probability")
                    or measured is None
                    or (
                        target == "key_bytes"
                        and prediction_scopes.get(target)
                        and prediction_scopes[target] != key_bytes_scope
                    )
                ):
                    continue
                try:
                    prediction_residuals[target] = float(measured) - float(predicted)
                except (TypeError, ValueError):
                    continue
            prediction_comparison_status = "paired"
        else:
            prediction_comparison_status = "unavailable_failed_measurement"
    return {
        **record,
        "prediction_residuals": (
            prediction_residuals
            if isinstance(predictions, dict)
            else record.get("prediction_residuals", {})
        ),
        "prediction_comparison_status": prediction_comparison_status,
        "measurement_kind": measurement_kind,
        "workload_group": dataset_group,
        "group_variant_count": f"variants:{variants}",
        "group_sample_count": f"samples:{samples}",
        "group_window_size": f"window:{window}",
        "group_parameter_combination": (
            f"{scheme}:{ring_dimension}:{window}:"
            f"{plan.get('scale_bits') or record.get('ckks_scaling_mod_size') or ''}:"
            f"{plan.get('fixed_point_scaling_factor') or record.get('fixed_point_scale') or ''}"
        ),
        "group_scheme_configuration": f"{scheme}:{ring_dimension}",
        # Do not relabel a generic dataset/source group as a PGS identifier.
        # Cross-PGS evaluation is available only when records carry real PGS IDs.
        "group_pgs": str(record.get("pgs_id") or ""),
        "group_source_domain": source_domain,
        "sample_count": samples,
        "matched_variant_count": variants,
        "total_pgs_variants": int(
            record.get("total_pgs_variants")
            or record.get("requested_weights_on_chrom")
            or variants
        ),
        "score_count": int(plan.get("scores_per_ciphertext") or 1),
        "scheme": scheme,
        "ring_dimension": ring_dimension,
        "security_bits": int(plan.get("security_bits") or 128),
        "multiplicative_depth": int(plan.get("multiplicative_depth") or 1),
        "first_modulus_size": plan.get("first_modulus_size"),
        "scaling_modulus_size": (
            plan.get("scaling_modulus_size")
            or record.get("ckks_scaling_mod_size")
        ),
        "plaintext_modulus": plan.get("plaintext_modulus") or record.get("plain_modulus"),
        "fixed_point_scale": (
            plan.get("fixed_point_scaling_factor") or record.get("fixed_point_scale")
        ),
        "ckks_scale_bits": (
            plan.get("scale_bits") or record.get("ckks_scaling_mod_size")
        ),
        "slots_used": window,
        "window_size": window,
        "chunk_count": int(plan.get("chunk_count") or 1),
        "samples_per_ciphertext": int(plan.get("samples_per_ciphertext") or 1),
        "scores_per_ciphertext": int(plan.get("scores_per_ciphertext") or 1),
        "key_count": (
            len(plan.get("rotation_indices") or ())
            + int(record.get("multiplication_key_bytes") not in {None, ""})
            + 2
        ),
        "packing_strategy": "feature" if plan.get("feature_packing", True) else "none",
        "reduction_strategy": plan.get("reduction_tree_strategy", "binary_tree"),
        "execution_success": int(record.get("status") == "ok"),
        "error_constraint_feasible": error_feasible,
        "total_seconds": total_seconds,
        "peak_memory_mb": peak_memory_mb,
        "ciphertext_bytes": ciphertext_bytes,
        "communication_bytes": record.get("communication_bytes") or ciphertext_bytes,
        "key_bytes": key_bytes,
        "key_bytes_scope": key_bytes_scope,
        "max_absolute_error": max_absolute_error,
        **operations,
    }


def build_dataset(records: Iterable[dict[str, Any]]) -> list[dict[str, Any]]:
    rows = [canonicalize(record) for record in records]
    return [
        row for row in rows
        if (
            str(row.get("measurement_kind", "")).startswith("measured")
            or row.get("measurement_kind") == "validation_only"
        )
        and row.get("status") is not None
    ]
