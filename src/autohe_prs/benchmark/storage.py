"""Canonical JSONL/CSV/optional Parquet benchmark storage."""

from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Any, Iterable


def _json_value(value: Any) -> Any:
    if isinstance(value, (dict, list, tuple)):
        return json.dumps(value, sort_keys=True, separators=(",", ":"))
    return value


def write_records(records: Iterable[dict[str, Any]], path: str | Path) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    rows = list(records)
    suffix = path.suffix.lower()
    if suffix == ".json":
        payload: Any = rows[0] if len(rows) == 1 else rows
        path.write_text(
            json.dumps(payload, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        return
    if suffix == ".jsonl":
        with path.open("w", encoding="utf-8") as handle:
            for row in rows:
                handle.write(json.dumps(row, sort_keys=True) + "\n")
        return
    if suffix == ".csv":
        fields = sorted({key for row in rows for key in row})
        with path.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields)
            writer.writeheader()
            writer.writerows({key: _json_value(value) for key, value in row.items()} for row in rows)
        return
    if suffix == ".parquet":
        try:
            import pandas as pd  # type: ignore[import-not-found]
        except ImportError as exc:
            raise RuntimeError("Parquet output requires the 'data' optional dependencies") from exc
        normalized = [
            {key: _json_value(value) for key, value in row.items()}
            for row in rows
        ]
        frame = pd.DataFrame(normalized)
        for column in frame.columns:
            if frame[column].dtype != "object":
                continue
            present = frame[column].dropna()
            present = present[present.astype(str) != ""]
            if present.empty:
                continue
            try:
                numeric = pd.to_numeric(present, errors="raise")
            except (TypeError, ValueError):
                frame[column] = frame[column].map(
                    lambda value: None if value is None else str(value)
                )
            else:
                frame.loc[present.index, column] = numeric
                frame[column] = pd.to_numeric(frame[column], errors="coerce")
        frame.to_parquet(path, index=False)
        return
    raise ValueError("benchmark output must use .json, .jsonl, .csv, or .parquet")


def read_records(path: str | Path) -> list[dict[str, Any]]:
    path = Path(path)
    suffix = path.suffix.lower()
    if suffix == ".json":
        payload = json.loads(path.read_text(encoding="utf-8"))
        if isinstance(payload, dict):
            return [payload]
        if isinstance(payload, list) and all(isinstance(row, dict) for row in payload):
            return payload
        raise ValueError("JSON benchmark input must contain one object or a list of objects")
    if suffix == ".jsonl":
        return [json.loads(line) for line in path.read_text(encoding="utf-8").splitlines() if line]
    if suffix == ".csv":
        with path.open("r", newline="", encoding="utf-8-sig") as handle:
            return list(csv.DictReader(handle))
    if suffix == ".parquet":
        try:
            import pandas as pd  # type: ignore[import-not-found]
        except ImportError as exc:
            raise RuntimeError("Parquet input requires the 'data' optional dependencies") from exc
        frame = pd.read_parquet(path)
        frame = frame.astype(object).where(frame.notna(), None)
        return frame.to_dict(orient="records")
    raise ValueError("benchmark input must use .json, .jsonl, .csv, or .parquet")
