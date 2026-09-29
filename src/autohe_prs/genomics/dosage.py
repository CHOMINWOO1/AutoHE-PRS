"""Preprocessed dosage-matrix input with allele-bearing column identifiers."""

from __future__ import annotations

import csv
from pathlib import Path
from typing import Any

from .vcf_parser import VCFData, VCFVariant


def _variant_header(value: str) -> tuple[str, int, str, str, str]:
    """Parse `chrom:position:REF:ALT[:rsID]`."""

    parts = value.split(":")
    if len(parts) not in {4, 5}:
        raise ValueError(
            f"dosage column must be chrom:position:REF:ALT[:rsID], got {value!r}"
        )
    chrom = parts[0][3:] if parts[0].lower().startswith("chr") else parts[0]
    try:
        position = int(parts[1])
    except ValueError as exc:
        raise ValueError(f"invalid dosage column position: {value!r}") from exc
    ref, alt = parts[2].upper(), parts[3].upper()
    return chrom, position, ref, alt, parts[4] if len(parts) == 5 else "."


def _from_rows(
    path: Path,
    fieldnames: list[str],
    rows: list[dict[str, Any]],
    genome_build: str | None,
) -> VCFData:
    if not fieldnames or fieldnames[0] != "sample_id":
        raise ValueError("dosage matrix must begin with sample_id")
    headers = [_variant_header(value) for value in fieldnames[1:]]
    sample_ids: list[str] = []
    columns: list[list[float | None]] = [[] for _ in headers]
    for row in rows:
        sample_ids.append(str(row["sample_id"]))
        for index, name in enumerate(fieldnames[1:]):
            raw = row.get(name)
            if raw in {None, "", ".", "NA", "NaN", "nan"}:
                columns[index].append(None)
                continue
            value = float(raw)
            columns[index].append(value if 0.0 <= value <= 2.0 else None)
    variants = tuple(
        VCFVariant(chrom, position, rsid, ref, alt, tuple(values), "matrix")
        for (chrom, position, ref, alt, rsid), values in zip(headers, columns)
    )
    return VCFData(
        path,
        tuple(sample_ids),
        variants,
        genome_build,
        len(variants),
        0,
    )


def read_dosage_matrix(
    path: str | Path, *, genome_build: str | None = None
) -> VCFData:
    path = Path(path)
    if path.suffix.lower() == ".parquet":
        try:
            import pandas as pd  # type: ignore[import-not-found]
        except ImportError as exc:
            raise RuntimeError("Parquet dosage input requires the 'data' dependencies") from exc
        frame = pd.read_parquet(path)
        return _from_rows(
            path,
            [str(column) for column in frame.columns],
            frame.to_dict(orient="records"),
            genome_build,
        )
    delimiter = "\t" if path.suffix.lower() in {".tsv", ".txt"} else ","
    with path.open("r", newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle, delimiter=delimiter)
        fieldnames = reader.fieldnames or []
        rows = list(reader)
    return _from_rows(path, fieldnames, rows, genome_build)
