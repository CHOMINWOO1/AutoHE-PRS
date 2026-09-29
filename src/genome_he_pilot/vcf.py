"""VCF dosage extraction for public genomic validation."""

from __future__ import annotations

from dataclasses import dataclass
import csv
import gzip
from pathlib import Path
import random


@dataclass(frozen=True)
class VariantRecord:
    chrom: str
    position: int
    variant_id: str
    ref: str
    alt: str
    allele_frequency: float


@dataclass(frozen=True)
class VCFGenomicDataset:
    """A public-data genotype dosage subset.

    Genotypes are encoded as alternate allele dosage: 0, 1, or 2.
    Rows are samples and columns are SNPs.
    """

    genotypes: list[list[int]]
    weights: list[float]
    phenotypes: list[int]
    allele_frequencies: list[float]
    sample_ids: list[str]
    variants: list[VariantRecord]
    source_path: Path

    @property
    def n_samples(self) -> int:
        return len(self.genotypes)

    @property
    def n_snps(self) -> int:
        return len(self.weights)


def read_sample_panel(panel_path: Path) -> dict[str, dict[str, str]]:
    """Read a 1000 Genomes sample panel keyed by sample ID."""

    with panel_path.open("r", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return {row["sample"]: row for row in reader}


def select_panel_samples(
    panel_path: Path,
    limit: int,
    population: str | None = None,
    super_population: str | None = None,
) -> list[str]:
    """Select sample IDs from a 1000 Genomes panel."""

    selected: list[str] = []
    panel = read_sample_panel(panel_path)
    for sample_id, row in panel.items():
        if population and row.get("pop") != population:
            continue
        if super_population and row.get("super_pop") != super_population:
            continue
        selected.append(sample_id)
        if len(selected) >= limit:
            break
    return selected


def make_sparse_weights(
    n_snps: int,
    informative_fraction: float = 0.05,
    seed: int = 20260721,
) -> list[float]:
    """Create reproducible sparse PRS-like weights for public genotypes."""

    if n_snps <= 0:
        raise ValueError("n_snps must be positive")
    if not 0 < informative_fraction <= 1:
        raise ValueError("informative_fraction must be in (0, 1]")

    rng = random.Random(seed)
    n_informative = max(1, int(n_snps * informative_fraction))
    informative_indices = set(rng.sample(range(n_snps), n_informative))
    return [
        rng.gauss(0.0, 0.08) if index in informative_indices else 0.0
        for index in range(n_snps)
    ]


def load_vcf_dosage_subset(
    vcf_path: Path,
    n_samples: int = 32,
    n_snps: int = 512,
    sample_ids: list[str] | None = None,
    seed: int = 20260721,
    informative_fraction: float = 0.05,
    min_minor_allele_count: int = 1,
) -> VCFGenomicDataset:
    """Load a small biallelic SNP dosage matrix from a gzipped VCF."""

    if n_samples <= 0:
        raise ValueError("n_samples must be positive")
    if n_snps <= 0:
        raise ValueError("n_snps must be positive")
    if min_minor_allele_count < 0:
        raise ValueError("min_minor_allele_count must be non-negative")

    selected_sample_ids: list[str] | None = None
    sample_field_indices: list[int] | None = None
    variant_columns: list[list[int]] = []
    variants: list[VariantRecord] = []

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
            if selected_sample_ids is None or sample_field_indices is None:
                raise ValueError("VCF header with sample IDs was not found")
            if len(variants) >= n_snps:
                break

            fields = line.rstrip("\n").split("\t")
            if len(fields) < 10:
                continue
            chrom, pos, variant_id, ref, alt = fields[0], fields[1], fields[2], fields[3], fields[4]
            if not _is_biallelic_snp(ref, alt):
                continue

            format_fields = fields[8].split(":")
            try:
                gt_index = format_fields.index("GT")
            except ValueError:
                continue

            dosages: list[int] = []
            missing = False
            sample_fields = fields[9:]
            for sample_index in sample_field_indices:
                dosage = _parse_gt_dosage(sample_fields[sample_index], gt_index)
                if dosage is None:
                    missing = True
                    break
                dosages.append(dosage)
            if missing:
                continue

            alt_count = sum(dosages)
            minor_allele_count = min(alt_count, 2 * len(dosages) - alt_count)
            if minor_allele_count < min_minor_allele_count:
                continue

            allele_frequency = alt_count / (2 * len(dosages))
            variant_columns.append(dosages)
            variants.append(
                VariantRecord(
                    chrom=chrom,
                    position=int(pos),
                    variant_id=variant_id,
                    ref=ref,
                    alt=alt,
                    allele_frequency=allele_frequency,
                )
            )

    if selected_sample_ids is None:
        raise ValueError("VCF header with sample IDs was not found")
    if len(variants) < n_snps:
        raise ValueError(f"requested {n_snps} SNPs but found {len(variants)} usable SNPs")

    genotypes = [
        [variant_columns[j][i] for j in range(len(variant_columns))]
        for i in range(len(selected_sample_ids))
    ]
    weights = make_sparse_weights(
        n_snps=len(variants),
        informative_fraction=informative_fraction,
        seed=seed,
    )
    raw_scores = [sum(g * w for g, w in zip(row, weights)) for row in genotypes]
    threshold = sorted(raw_scores)[len(raw_scores) // 2]
    phenotypes = [int(score >= threshold) for score in raw_scores]

    return VCFGenomicDataset(
        genotypes=genotypes,
        weights=weights,
        phenotypes=phenotypes,
        allele_frequencies=[variant.allele_frequency for variant in variants],
        sample_ids=selected_sample_ids,
        variants=variants,
        source_path=vcf_path,
    )


def _select_vcf_samples(
    vcf_sample_ids: list[str],
    n_samples: int,
    requested_sample_ids: list[str] | None,
) -> tuple[list[str], list[int]]:
    if requested_sample_ids is None:
        selected = vcf_sample_ids[:n_samples]
        return selected, list(range(len(selected)))

    lookup = {sample_id: index for index, sample_id in enumerate(vcf_sample_ids)}
    selected = []
    indices = []
    for sample_id in requested_sample_ids:
        if sample_id not in lookup:
            continue
        selected.append(sample_id)
        indices.append(lookup[sample_id])
        if len(selected) >= n_samples:
            break

    if len(selected) < n_samples:
        raise ValueError(
            f"requested {n_samples} samples but only found {len(selected)} in VCF"
        )
    return selected, indices


def _is_biallelic_snp(ref: str, alt: str) -> bool:
    return len(ref) == 1 and len(alt) == 1 and "," not in alt


def _parse_gt_dosage(sample_field: str, gt_index: int) -> int | None:
    parts = sample_field.split(":")
    if gt_index >= len(parts):
        return None
    genotype = parts[gt_index]
    if "." in genotype:
        return None
    alleles = genotype.replace("|", "/").split("/")
    if not alleles:
        return None

    dosage = 0
    for allele in alleles:
        if allele not in {"0", "1"}:
            return None
        dosage += int(allele)
    return dosage
