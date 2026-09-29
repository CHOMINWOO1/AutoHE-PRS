"""Paper-ready Markdown tables from benchmark dictionaries."""

from __future__ import annotations

from pathlib import Path
from typing import Any


def write_scheme_table(records: list[dict[str, Any]], path: str | Path) -> Path:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    header = (
        "| Scheme | Provenance | Variants | Samples | Status | Evaluation seconds | "
        "Peak MB | Ciphertext bytes | Key bytes | Key scope | Max abs error |\n"
        "| --- | --- | ---: | ---: | --- | ---: | ---: | ---: | ---: | --- | ---: |"
    )
    lines = [header]
    for row in sorted(
        records,
        key=lambda value: (
            int(value.get("variants") or value.get("matched_variant_count") or 0),
            str(value.get("scheme", "")),
        ),
    ):
        lines.append(
            "| {scheme} | {provenance} | {variants} | {samples} | {status} | {evaluation} | "
            "{memory} | {ciphertext} | {key} | {key_scope} | {error} |".format(
                scheme=row.get("scheme", ""),
                provenance=row.get("measurement_kind", "unknown"),
                variants=row.get("variants", row.get("matched_variant_count", "")),
                samples=row.get("samples", row.get("sample_count", "")),
                status=row.get("status", ""),
                evaluation=row.get("evaluation_seconds", ""),
                memory=row.get("peak_memory_mb", ""),
                ciphertext=row.get("ciphertext_bytes", ""),
                key=row.get("key_bytes", ""),
                key_scope=row.get("key_bytes_scope", ""),
                error=row.get("max_absolute_error", ""),
            )
        )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def write_prediction_comparison_table(
    records: list[dict[str, Any]], path: str | Path
) -> Path:
    """Write paired predicted-vs-measured values; never synthesize a pair."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    header = (
        "| Scheme | Rank | Stage | Target | Predicted | Measured | Residual |\n"
        "| --- | ---: | --- | --- | ---: | ---: | ---: |"
    )
    lines = [header]
    for row in records:
        if (
            row.get("status") != "ok"
            or row.get("prediction_comparison_status") != "paired"
        ):
            continue
        predicted = row.get("predicted_targets")
        if not isinstance(predicted, dict):
            continue
        residuals = row.get("prediction_residuals")
        residuals = residuals if isinstance(residuals, dict) else {}
        prediction_scopes = row.get("prediction_target_scopes")
        prediction_scopes = (
            prediction_scopes if isinstance(prediction_scopes, dict) else {}
        )
        for target, predicted_value in sorted(predicted.items()):
            if target.endswith("_probability") or row.get(target) is None:
                continue
            if (
                target == "key_bytes"
                and prediction_scopes.get(target)
                and prediction_scopes[target] != row.get("key_bytes_scope")
            ):
                continue
            lines.append(
                "| {scheme} | {rank} | {stage} | {target} | {predicted} | "
                "{measured} | {residual} |".format(
                    scheme=row.get("scheme", ""),
                    rank=row.get("predicted_rank", ""),
                    stage=row.get("search_stage", ""),
                    target=target,
                    predicted=predicted_value,
                    measured=row[target],
                    residual=residuals.get(target, ""),
                )
            )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path
