"""Streaming VCF parser for GT hard calls and DS imputed dosages."""

from __future__ import annotations

from dataclasses import dataclass
import gzip
from pathlib import Path
from typing import TextIO


@dataclass(frozen=True)
class VCFVariant:
    chrom: str
    position: int
    rsid: str
    ref: str
    alt: str
    alt_dosages: tuple[float | None, ...]
    dosage_source: str

    @property
    def key(self) -> tuple[str, int]:
        return self.chrom, self.position


@dataclass(frozen=True)
class VCFData:
    path: Path
    sample_ids: tuple[str, ...]
    variants: tuple[VCFVariant, ...]
    genome_build: str | None = None
    total_variant_records: int = 0
    filtered_out_records: int = 0


def _open(path: Path) -> TextIO:
    if path.suffix.lower() == ".gz":
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open("r", encoding="utf-8")


def _chrom(value: str) -> str:
    value = value.strip()
    return value[3:] if value.lower().startswith("chr") else value


def gt_to_alt_dosage(value: str) -> float | None:
    if not value or "." in value:
        return None
    alleles = value.replace("|", "/").split("/")
    if len(alleles) != 2 or any(allele not in {"0", "1"} for allele in alleles):
        return None
    return float(sum(int(allele) for allele in alleles))


def _sample_dosage(sample: str, fields: list[str], prefer_ds: bool) -> tuple[float | None, str]:
    values = sample.split(":")
    lookup = dict(zip(fields, values))
    if prefer_ds and "DS" in lookup and lookup["DS"] not in {"", "."}:
        try:
            dosage = float(lookup["DS"])
        except ValueError:
            return None, "DS"
        return (dosage if 0.0 <= dosage <= 2.0 else None), "DS"
    return gt_to_alt_dosage(lookup.get("GT", "")), "GT"


def read_vcf(
    path: str | Path,
    *,
    sample_ids: list[str] | None = None,
    sample_limit: int | None = None,
    variant_keys: set[tuple[str, int]] | None = None,
    prefer_ds: bool = True,
    genome_build: str | None = None,
) -> VCFData:
    path = Path(path)
    selected_ids: tuple[str, ...] | None = None
    selected_indices: list[int] = []
    variants: list[VCFVariant] = []
    total_variant_records = 0
    filtered_out_records = 0
    with _open(path) as handle:
        for raw in handle:
            if raw.startswith("##"):
                continue
            if raw.startswith("#CHROM"):
                header_samples = raw.rstrip("\r\n").split("\t")[9:]
                if sample_ids is None:
                    if sample_limit is not None:
                        header_samples = header_samples[:sample_limit]
                    selected_ids = tuple(header_samples)
                    selected_indices = list(range(len(header_samples)))
                else:
                    index = {sample: i for i, sample in enumerate(header_samples)}
                    missing = [sample for sample in sample_ids if sample not in index]
                    if missing:
                        raise ValueError(f"samples absent from VCF: {', '.join(missing)}")
                    selected_ids = tuple(sample_ids)
                    selected_indices = [index[sample] for sample in sample_ids]
                continue
            if raw.startswith("#"):
                continue
            if selected_ids is None:
                raise ValueError("VCF #CHROM header is missing")
            parts = raw.rstrip("\r\n").split("\t")
            if len(parts) < 9:
                continue
            total_variant_records += 1
            chrom, position, rsid, ref, alt = parts[:5]
            normalized_key = (_chrom(chrom), int(position))
            if variant_keys is not None and normalized_key not in variant_keys:
                filtered_out_records += 1
                continue
            if "," in alt:
                # Retain the unsupported record so harmonization can report it.
                variants.append(
                    VCFVariant(normalized_key[0], normalized_key[1], rsid, ref.upper(), alt.upper(),
                               tuple(None for _ in selected_ids), "unsupported_multiallelic")
                )
                continue
            format_fields = parts[8].split(":")
            dosage_values: list[float | None] = []
            sources: set[str] = set()
            sample_values = parts[9:]
            for sample_index in selected_indices:
                if sample_index >= len(sample_values):
                    dosage, source = None, "missing"
                else:
                    dosage, source = _sample_dosage(
                        sample_values[sample_index], format_fields, prefer_ds
                    )
                dosage_values.append(dosage)
                sources.add(source)
            dosage_source = next(iter(sources)) if len(sources) == 1 else "mixed"
            variants.append(
                VCFVariant(
                    normalized_key[0],
                    normalized_key[1],
                    rsid,
                    ref.upper(),
                    alt.upper(),
                    tuple(dosage_values),
                    dosage_source,
                )
            )
    if selected_ids is None:
        raise ValueError("VCF #CHROM header is missing")
    return VCFData(
        path,
        selected_ids,
        tuple(variants),
        genome_build,
        total_variant_records,
        filtered_out_records,
    )
