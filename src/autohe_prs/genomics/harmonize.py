"""Deterministic allele harmonization with explicit exclusion provenance."""

from __future__ import annotations

from collections import Counter
from dataclasses import asdict, dataclass
import csv
import json
from pathlib import Path
from typing import Iterable

from .pgs_parser import PGSFile, PGSVariant
from .vcf_parser import VCFData, VCFVariant


AMBIGUOUS = {frozenset(("A", "T")), frozenset(("C", "G"))}
COMPLEMENT = {"A": "T", "T": "A", "C": "G", "G": "C"}
MISSING_POLICIES = {"error", "exclude_variant", "zero", "mean_dosage"}


@dataclass(frozen=True)
class HarmonizedVariant:
    chrom: str
    position: int
    rsid: str
    ref: str
    alt: str
    effect_allele: str
    orientation: str
    effect_weight: float
    effect_dosages: tuple[float, ...]
    dosage_source: str
    strand_flipped: bool = False


@dataclass(frozen=True)
class ExcludedVariant:
    chrom: str
    position: int
    rsid: str
    reason: str
    detail: str = ""


@dataclass(frozen=True)
class HarmonizationReport:
    total_pgs_variants: int
    matched_variants: int
    alt_oriented_variants: int
    ref_oriented_variants: int
    strand_ambiguous_variants: int
    strand_flipped_variants: int
    multiallelic_exclusions: int
    missing_variants: int
    duplicate_variants: int
    unsupported_variants: int
    final_scored_variants: int
    vcf_only_variants: int
    exclusion_reason_counts: dict[str, int]
    missing_policy: str


@dataclass(frozen=True)
class HarmonizationResult:
    sample_ids: tuple[str, ...]
    variants: tuple[HarmonizedVariant, ...]
    excluded: tuple[ExcludedVariant, ...]
    report: HarmonizationReport

    @property
    def dosage_matrix(self) -> list[list[float]]:
        return [
            [variant.effect_dosages[sample] for variant in self.variants]
            for sample in range(len(self.sample_ids))
        ]

    @property
    def weights(self) -> list[float]:
        return [variant.effect_weight for variant in self.variants]


def _excluded(weight: PGSVariant, reason: str, detail: str = "") -> ExcludedVariant:
    return ExcludedVariant(weight.chrom, weight.position, weight.rsid, reason, detail)


def _resolve_missing(
    values: tuple[float | None, ...], policy: str, weight: PGSVariant
) -> tuple[float, ...] | None:
    if all(value is not None for value in values):
        return tuple(float(value) for value in values if value is not None)
    if policy == "error":
        raise ValueError(
            f"missing genotype at {weight.chrom}:{weight.position}; "
            "select an explicit missing policy"
        )
    if policy == "exclude_variant":
        return None
    observed = [float(value) for value in values if value is not None]
    fill = 0.0 if policy == "zero" else (sum(observed) / len(observed) if observed else 0.0)
    return tuple(fill if value is None else float(value) for value in values)


def harmonize(
    vcf: VCFData,
    pgs: PGSFile,
    *,
    missing_policy: str = "error",
    allow_strand_ambiguous: bool = False,
) -> HarmonizationResult:
    if missing_policy not in MISSING_POLICIES:
        raise ValueError(f"unsupported missing policy: {missing_policy}")
    pgs_build = next((variant.genome_build for variant in pgs.variants if variant.genome_build), None)
    if vcf.genome_build and pgs_build and vcf.genome_build.lower() != pgs_build.lower():
        raise ValueError(f"genome build mismatch: VCF={vcf.genome_build}, PGS={pgs_build}")

    index: dict[tuple[str, int], list[VCFVariant]] = {}
    for variant in vcf.variants:
        index.setdefault(variant.key, []).append(variant)
    duplicate_pgs = Counter(variant.key for variant in pgs.variants)
    matched: list[HarmonizedVariant] = []
    excluded: list[ExcludedVariant] = [
        ExcludedVariant("", item.source_row, "", item.reason, item.detail)
        for item in pgs.exclusions
    ]
    used_vcf_keys: set[tuple[str, int]] = set()

    for weight in pgs.variants:
        if duplicate_pgs[weight.key] > 1:
            excluded.append(_excluded(weight, "duplicate_pgs_variant"))
            continue
        candidates = index.get(weight.key, [])
        if not candidates:
            excluded.append(_excluded(weight, "missing_from_vcf"))
            continue
        if len(candidates) > 1:
            excluded.append(_excluded(weight, "duplicate_vcf_variant"))
            continue
        variant = candidates[0]
        used_vcf_keys.add(variant.key)
        if "," in variant.alt or variant.dosage_source == "unsupported_multiallelic":
            excluded.append(_excluded(weight, "multiallelic_variant", variant.alt))
            continue
        allele_pair = frozenset((variant.ref, variant.alt))
        if allele_pair in AMBIGUOUS and not allow_strand_ambiguous:
            excluded.append(_excluded(weight, "strand_ambiguous", f"{variant.ref}/{variant.alt}"))
            continue
        effect_for_match = weight.effect_allele
        strand_flipped = False
        if weight.other_allele and {weight.effect_allele, weight.other_allele} != {
            variant.ref,
            variant.alt,
        }:
            complemented = {
                COMPLEMENT.get(weight.effect_allele, ""),
                COMPLEMENT.get(weight.other_allele, ""),
            }
            if complemented == {variant.ref, variant.alt}:
                effect_for_match = COMPLEMENT[weight.effect_allele]
                strand_flipped = True
            else:
                excluded.append(_excluded(weight, "allele_mismatch", f"{variant.ref}/{variant.alt}"))
                continue
        alt_values = _resolve_missing(variant.alt_dosages, missing_policy, weight)
        if alt_values is None:
            excluded.append(_excluded(weight, "missing_genotype"))
            continue
        if effect_for_match == variant.alt:
            effect_values = alt_values
            orientation = "ALT"
        elif effect_for_match == variant.ref:
            effect_values = tuple(2.0 - value for value in alt_values)
            orientation = "REF"
        else:
            excluded.append(_excluded(weight, "allele_mismatch", weight.effect_allele))
            continue
        matched.append(
            HarmonizedVariant(
                weight.chrom,
                weight.position,
                variant.rsid if variant.rsid != "." else weight.rsid,
                variant.ref,
                variant.alt,
                weight.effect_allele,
                orientation,
                weight.effect_weight,
                effect_values,
                variant.dosage_source,
                strand_flipped,
            )
        )

    reasons = Counter(item.reason for item in excluded)
    outside_target = int(pgs.metadata.get("autohe_outside_target_chromosome", "0"))
    if outside_target:
        # One representative exclusion row is retained while the report carries
        # the exact aggregate count, avoiding millions of duplicate objects.
        reasons["outside_target_chromosome"] += outside_target - 1
    total_pgs_rows = int(
        pgs.metadata.get(
            "autohe_total_data_rows",
            str(len(pgs.variants) + len(pgs.exclusions)),
        )
    )
    report = HarmonizationReport(
        total_pgs_variants=total_pgs_rows,
        matched_variants=len(matched),
        alt_oriented_variants=sum(item.orientation == "ALT" for item in matched),
        ref_oriented_variants=sum(item.orientation == "REF" for item in matched),
        strand_ambiguous_variants=reasons["strand_ambiguous"],
        strand_flipped_variants=sum(item.strand_flipped for item in matched),
        multiallelic_exclusions=reasons["multiallelic_variant"]
        + reasons["multiallelic_or_haplotype"],
        missing_variants=reasons["missing_from_vcf"] + reasons["missing_genotype"],
        duplicate_variants=reasons["duplicate_pgs_variant"] + reasons["duplicate_vcf_variant"],
        unsupported_variants=sum(
            count
            for reason, count in reasons.items()
            if reason.startswith("unsupported") or reason == "incomplete_variant"
        ),
        final_scored_variants=len(matched),
        vcf_only_variants=(
            vcf.filtered_out_records
            + len({variant.key for variant in vcf.variants} - used_vcf_keys)
        ),
        exclusion_reason_counts=dict(sorted(reasons.items())),
        missing_policy=missing_policy,
    )
    return HarmonizationResult(vcf.sample_ids, tuple(matched), tuple(excluded), report)


def write_harmonized(
    result: HarmonizationResult,
    output_dir: str | Path,
    *,
    output_format: str = "csv",
) -> None:
    if output_format not in {"csv", "parquet", "both"}:
        raise ValueError("harmonized output format must be csv, parquet, or both")
    output = Path(output_dir)
    output.mkdir(parents=True, exist_ok=True)
    weight_fields = [
        "chrom", "position", "rsid", "ref", "alt", "effect_allele",
        "orientation", "effect_weight", "dosage_source", "strand_flipped",
    ]
    weight_rows = []
    for variant in result.variants:
        row = asdict(variant)
        row.pop("effect_dosages")
        weight_rows.append(row)
    if output_format in {"csv", "both"}:
        with (output / "dosage.csv").open("w", newline="", encoding="utf-8") as handle:
            writer = csv.writer(handle)
            writer.writerow([
                "sample_id",
                *[
                    f"{v.chrom}:{v.position}:{v.ref}:{v.alt}:{v.rsid or '.'}"
                    for v in result.variants
                ],
            ])
            for sample_id, values in zip(result.sample_ids, result.dosage_matrix):
                writer.writerow([sample_id, *values])
        with (output / "weights.csv").open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=weight_fields)
            writer.writeheader()
            writer.writerows(weight_rows)
    if output_format in {"parquet", "both"}:
        try:
            import pandas as pd  # type: ignore[import-not-found]
        except ImportError as exc:
            raise RuntimeError("Parquet harmonized output requires the 'data' dependencies") from exc
        dosage_rows = [
            {"sample_id": sample_id, **{
                (
                    f"{variant.chrom}:{variant.position}:{variant.ref}:"
                    f"{variant.alt}:{variant.rsid or '.'}"
                ): value
                for variant, value in zip(result.variants, values)
            }}
            for sample_id, values in zip(result.sample_ids, result.dosage_matrix)
        ]
        pd.DataFrame(dosage_rows).to_parquet(output / "dosage.parquet", index=False)
        pd.DataFrame(weight_rows, columns=weight_fields).to_parquet(
            output / "weights.parquet", index=False
        )
    (output / "harmonization_report.json").write_text(
        json.dumps(asdict(result.report), indent=2, sort_keys=True), encoding="utf-8"
    )
    with (output / "exclusions.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(ExcludedVariant.__dataclass_fields__))
        writer.writeheader()
        writer.writerows(asdict(item) for item in result.excluded)
