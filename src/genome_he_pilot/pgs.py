"""PGS Catalog scoring-file parsing utilities."""

from __future__ import annotations

from dataclasses import dataclass
import gzip
from pathlib import Path


@dataclass(frozen=True)
class PGSWeight:
    pgs_id: str
    chrom: str
    position: int
    rsid: str
    effect_allele: str
    other_allele: str
    effect_weight: float


@dataclass(frozen=True)
class PGSScoringFile:
    path: Path
    metadata: dict[str, str]
    weights: list[PGSWeight]


def read_pgs_scoring_file(path: Path, chrom: str | None = None) -> PGSScoringFile:
    """Read a PGS Catalog scoring file in formatted or harmonized text format."""

    metadata: dict[str, str] = {}
    weights: list[PGSWeight] = []
    header: list[str] | None = None
    target_chrom = chrom.removeprefix("chr") if chrom else None

    with _open_text(path) as handle:
        for raw_line in handle:
            line = raw_line.rstrip("\n")
            if not line:
                continue
            if line.startswith("#"):
                _read_metadata_line(line, metadata)
                continue
            fields = line.split("\t")
            if header is None:
                header = fields
                continue
            row = dict(zip(header, fields))
            weight = _parse_weight_row(row, metadata.get("pgs_id", ""))
            if weight is not None and (target_chrom is None or weight.chrom == target_chrom):
                weights.append(weight)

    if header is None:
        raise ValueError(f"PGS scoring header was not found: {path}")
    if not weights:
        suffix = f" on chromosome {target_chrom}" if target_chrom else ""
        raise ValueError(f"no usable PGS weights found{suffix}: {path}")

    return PGSScoringFile(path=path, metadata=metadata, weights=weights)


def _open_text(path: Path):
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open("r", encoding="utf-8")


def _read_metadata_line(line: str, metadata: dict[str, str]) -> None:
    trimmed = line.lstrip("#")
    if "=" not in trimmed:
        return
    key, value = trimmed.split("=", 1)
    metadata[key.strip()] = value.strip()


def _parse_weight_row(row: dict[str, str], pgs_id: str) -> PGSWeight | None:
    chrom = _clean_chrom(row.get("chr_name", ""))
    position = row.get("chr_position", "")
    effect_allele = row.get("effect_allele", "")
    effect_weight = row.get("effect_weight", "")
    if not chrom or not position or not effect_allele or not effect_weight:
        return None

    try:
        parsed_position = int(position)
        parsed_weight = float(effect_weight)
    except ValueError:
        return None

    return PGSWeight(
        pgs_id=pgs_id,
        chrom=chrom,
        position=parsed_position,
        rsid=row.get("rsID", ""),
        effect_allele=effect_allele.upper(),
        other_allele=row.get("other_allele", row.get("reference_allele", "")).upper(),
        effect_weight=parsed_weight,
    )


def _clean_chrom(chrom: str) -> str:
    cleaned = chrom.strip()
    if cleaned.lower().startswith("chr"):
        cleaned = cleaned[3:]
    return cleaned
