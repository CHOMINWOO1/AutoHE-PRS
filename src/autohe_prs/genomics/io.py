"""Portable CSV I/O for harmonized dosage and weights."""

from __future__ import annotations

import csv
from pathlib import Path


def read_harmonized(
    dosage_path: str | Path, weight_path: str | Path
) -> tuple[list[str], list[list[float]], list[float]]:
    dosage_path, weight_path = Path(dosage_path), Path(weight_path)
    if dosage_path.suffix.lower() == ".parquet":
        try:
            import pandas as pd  # type: ignore[import-not-found]
        except ImportError as exc:
            raise RuntimeError("Parquet harmonized input requires the 'data' dependencies") from exc
        frame = pd.read_parquet(dosage_path)
        if not len(frame.columns) or str(frame.columns[0]) != "sample_id":
            raise ValueError("dosage Parquet must begin with sample_id")
        sample_ids = [str(value) for value in frame.iloc[:, 0].tolist()]
        dosages = [[float(value) for value in row] for row in frame.iloc[:, 1:].values.tolist()]
    else:
        with dosage_path.open("r", newline="", encoding="utf-8-sig") as handle:
            reader = csv.reader(handle)
            header = next(reader, None)
            if not header or header[0] != "sample_id":
                raise ValueError("dosage CSV must begin with a sample_id column")
            sample_ids = []
            dosages = []
            for row in reader:
                sample_ids.append(row[0])
                dosages.append([float(value) for value in row[1:]])
    if weight_path.suffix.lower() == ".parquet":
        try:
            import pandas as pd  # type: ignore[import-not-found]
        except ImportError as exc:
            raise RuntimeError("Parquet harmonized input requires the 'data' dependencies") from exc
        frame = pd.read_parquet(weight_path)
        if "effect_weight" not in frame.columns:
            raise ValueError("weights Parquet requires effect_weight")
        weights = [float(value) for value in frame["effect_weight"].tolist()]
    else:
        with weight_path.open("r", newline="", encoding="utf-8-sig") as handle:
            reader = csv.DictReader(handle)
            if not reader.fieldnames or "effect_weight" not in reader.fieldnames:
                raise ValueError("weights CSV requires effect_weight")
            weights = [float(row["effect_weight"]) for row in reader]
    if any(len(row) != len(weights) for row in dosages):
        raise ValueError("dosage and weight variant counts differ")
    return sample_ids, dosages, weights
