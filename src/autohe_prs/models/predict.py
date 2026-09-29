from __future__ import annotations

from pathlib import Path
from typing import Any

from .auto import load_cost_model


def predict_rows(model_path: str | Path, rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    model = load_cost_model(model_path)
    return [
        {
            **row,
            **{f"predicted_{key}": value for key, value in model.predict(row).items()},
            "measurement_kind": "predicted",
        }
        for row in rows
    ]
