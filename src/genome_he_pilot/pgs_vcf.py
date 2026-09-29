"""Match PGS Catalog weights to VCF genotype dosage."""

from __future__ import annotations

from dataclasses import dataclass
import gzip
from pathlib import Path

from .pgs import PGSScoringFile, PGSWeight
from .vcf import (
    VCFGenomicDataset,
    VariantRecord,
    _is_biallelic_snp,
    _parse_gt_dosage,
    _select_vcf_samples,
)


@dataclass(frozen=True)
class PGSVCFMatchResult:
    dataset: VCFGenomicDataset
    matched_weights: int
    skipped_allele_mismatch: int
    requested_weights_on_chrom: int


def load_pgs_matched_vcf_dataset(
    vcf_path: Path,
    scoring: PGSScoringFile,
    chrom: str = "22",
    n_samples: int = 32,
    sample_ids: list[str] | None = None,
    max_snps: int | None = None,
) -> PGSVCFMatchResult:
    """Load effect-allele dosage for PGS variants found in a VCF."""

    target_chrom = chrom.removeprefix("chr")
    weights_by_position: dict[int, list[PGSWeight]] = {}
    for weight in scoring.weights:
        if weight.chrom != target_chrom:
            continue
        weights_by_position.setdefault(weight.position, []).append(weight)

    if not weights_by_position:
        raise ValueError(f"no PGS weights found on chromosome {target_chrom}")

    selected_sample_ids: list[str] | None = None
    sample_field_indices: list[int] | None = None
    variant_columns: list[list[int]] = []
    selected_weights: list[float] = []
    variants: list[VariantRecord] = []
    skipped_allele_mismatch = 0

    with gzip.open(vcf_path, "rt", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                header = line.rstrip("\n").split("\t")
                vcf_sample_ids = header[9:]
                selected_sample_ids, sample_field_indices = _select_vcf_samples(
                    vcf_sample_ids=vcf_sample_ids,
                    n_samples=n_samples,
                    requested_sample_ids=sample_ids,
                )
                continue
            if line.startswith("#"):
                continue
            if max_snps is not None and len(selected_weights) >= max_snps:
                break
            if selected_sample_ids is None or sample_field_indices is None:
                raise ValueError("VCF header with sample IDs was not found")

            fields = line.rstrip("\n").split("\t")
            if len(fields) < 10:
                continue
            vcf_chrom = fields[0].removeprefix("chr")
            if vcf_chrom != target_chrom:
                continue

            position = int(fields[1])
            if position not in weights_by_position:
                continue

            variant_id, ref, alt = fields[2], fields[3].upper(), fields[4].upper()
            if not _is_biallelic_snp(ref, alt):
                continue

            format_fields = fields[8].split(":")
            try:
                gt_index = format_fields.index("GT")
            except ValueError:
                continue

            alt_dosages: list[int] = []
            missing = False
            sample_fields = fields[9:]
            for sample_index in sample_field_indices:
                dosage = _parse_gt_dosage(sample_fields[sample_index], gt_index)
                if dosage is None:
                    missing = True
                    break
                alt_dosages.append(dosage)
            if missing:
                continue

            for weight in weights_by_position[position]:
                if weight.effect_allele == alt:
                    effect_dosages = alt_dosages
                elif weight.effect_allele == ref:
                    effect_dosages = [2 - dosage for dosage in alt_dosages]
                else:
                    skipped_allele_mismatch += 1
                    continue

                effect_frequency = sum(effect_dosages) / (2 * len(effect_dosages))
                variant_columns.append(effect_dosages)
                selected_weights.append(weight.effect_weight)
                variants.append(
                    VariantRecord(
                        chrom=vcf_chrom,
                        position=position,
                        variant_id=variant_id if variant_id != "." else weight.rsid,
                        ref=ref,
                        alt=alt,
                        allele_frequency=effect_frequency,
                    )
                )
                if max_snps is not None and len(selected_weights) >= max_snps:
                    break

    if selected_sample_ids is None:
        raise ValueError("VCF header with sample IDs was not found")
    if not variant_columns:
        raise ValueError("no PGS variants matched the VCF with compatible alleles")

    genotypes = [
        [variant_columns[j][i] for j in range(len(variant_columns))]
        for i in range(len(selected_sample_ids))
    ]
    raw_scores = [sum(g * w for g, w in zip(row, selected_weights)) for row in genotypes]
    threshold = sorted(raw_scores)[len(raw_scores) // 2]
    phenotypes = [int(score >= threshold) for score in raw_scores]

    return PGSVCFMatchResult(
        dataset=VCFGenomicDataset(
            genotypes=genotypes,
            weights=selected_weights,
            phenotypes=phenotypes,
            allele_frequencies=[variant.allele_frequency for variant in variants],
            sample_ids=selected_sample_ids,
            variants=variants,
            source_path=vcf_path,
        ),
        matched_weights=len(selected_weights),
        skipped_allele_mismatch=skipped_allele_mismatch,
        requested_weights_on_chrom=sum(
            len(weights) for weights in weights_by_position.values()
        ),
    )
