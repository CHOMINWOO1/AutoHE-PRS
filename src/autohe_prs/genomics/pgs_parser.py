"""PGS Catalog formatted/harmonized scoring-file parser."""

from __future__ import annotations

from dataclasses import dataclass
import gzip
import math
from pathlib import Path
from typing import TextIO


@dataclass(frozen=True)
class PGSVariant:
    pgs_id: str
    genome_build: str | None
    chrom: str
    position: int
    rsid: str
    effect_allele: str
    other_allele: str
    effect_weight: float
    weight_type: str
    source_row: int

    @property
    def key(self) -> tuple[str, int]:
        return self.chrom, self.position


@dataclass(frozen=True)
class PGSExclusion:
    source_row: int
    reason: str
    detail: str


@dataclass(frozen=True)
class PGSFile:
    path: Path
    metadata: dict[str, str]
    variants: tuple[PGSVariant, ...]
    exclusions: tuple[PGSExclusion, ...]


def _open(path: Path) -> TextIO:
    if path.suffix.lower() == ".gz":
        return gzip.open(path, "rt", encoding="utf-8-sig")
    return path.open("r", encoding="utf-8-sig")


def _chrom(value: str) -> str:
    value = value.strip()
    return value[3:] if value.lower().startswith("chr") else value


def read_pgs(path: str | Path, *, chrom: str | None = None) -> PGSFile:
    path = Path(path)
    metadata: dict[str, str] = {}
    variants: list[PGSVariant] = []
    exclusions: list[PGSExclusion] = []
    header: list[str] | None = None
    target_chrom = _chrom(chrom) if chrom else None
    data_rows = 0
    outside_target = 0
    with _open(path) as handle:
        for source_row, raw in enumerate(handle, start=1):
            line = raw.rstrip("\r\n")
            if not line:
                continue
            if line.startswith("#"):
                body = line.lstrip("#")
                if "=" in body:
                    key, value = body.split("=", 1)
                    metadata[key.strip()] = value.strip()
                continue
            if header is None:
                header = line.split("\t")
                continue
            row = dict(zip(header, line.split("\t")))
            data_rows += 1
            pgs_id = metadata.get("pgs_id", path.name.split(".", 1)[0])
            build = metadata.get("HmPOS_build") or metadata.get("genome_build")
            chrom = _chrom(row.get("chr_name", row.get("chromosome", "")) or "")
            if target_chrom is not None and chrom != target_chrom:
                outside_target += 1
                continue
            position = row.get("chr_position", row.get("position", "")) or ""
            effect = (row.get("effect_allele", "") or "").upper()
            other = (row.get("other_allele", row.get("reference_allele", "")) or "").upper()
            weight = row.get("effect_weight", "") or ""
            weight_type = (
                row.get("weight_type", metadata.get("weight_type", "additive")) or "additive"
            ).lower()
            model_description = " ".join(
                str(row.get(name, "") or "").lower()
                for name in (
                    "effect_type",
                    "variant_description",
                    "locus_name",
                    "inclusion_criteria",
                )
            )
            unsupported_flags = [
                name
                for name in (
                    "is_interaction",
                    "is_dominant",
                    "is_recessive",
                    "is_haplotype",
                )
                if str(row.get(name, "") or "").strip().lower()
                in {"1", "true", "yes", "y"}
            ]
            if (
                unsupported_flags
                or any(
                    token in f"{weight_type} {model_description}"
                    for token in ("dominant", "recessive", "interaction", "haplotype")
                )
            ):
                exclusions.append(PGSExclusion(source_row, "unsupported_weight_type", weight_type))
                continue
            if not chrom or not position or not effect or not other or not weight:
                exclusions.append(PGSExclusion(source_row, "incomplete_variant", "required field missing"))
                continue
            try:
                parsed_position = int(position)
                parsed_weight = float(weight)
            except ValueError:
                exclusions.append(PGSExclusion(source_row, "invalid_numeric_value", f"{position}/{weight}"))
                continue
            if (
                len(effect) != 1
                or len(other) != 1
                or effect not in {"A", "C", "G", "T"}
                or other not in {"A", "C", "G", "T"}
                or effect == other
            ):
                exclusions.append(PGSExclusion(source_row, "multiallelic_or_haplotype", f"{effect}/{other}"))
                continue
            if parsed_position <= 0 or not math.isfinite(parsed_weight):
                exclusions.append(
                    PGSExclusion(
                        source_row,
                        "invalid_numeric_value",
                        f"{position}/{weight}",
                    )
                )
                continue
            variants.append(
                PGSVariant(
                    pgs_id=pgs_id,
                    genome_build=build,
                    chrom=chrom,
                    position=parsed_position,
                    rsid=row.get("rsID", row.get("rsid", "")) or "",
                    effect_allele=effect,
                    other_allele=other,
                    effect_weight=parsed_weight,
                    weight_type=weight_type,
                    source_row=source_row,
                )
            )
    if header is None:
        raise ValueError(f"PGS header is missing: {path}")
    metadata["autohe_total_data_rows"] = str(data_rows)
    metadata["autohe_outside_target_chromosome"] = str(outside_target)
    if target_chrom is not None:
        metadata["autohe_target_chromosome"] = target_chrom
    if outside_target:
        exclusions.append(
            PGSExclusion(
                0,
                "outside_target_chromosome",
                f"{outside_target} rows excluded by --chrom {target_chrom}",
            )
        )
    return PGSFile(path, metadata, tuple(variants), tuple(exclusions))
