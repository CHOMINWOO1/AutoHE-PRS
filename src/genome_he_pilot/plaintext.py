"""Plaintext genomic statistics used as validation baselines."""

from __future__ import annotations

from dataclasses import dataclass
import math


@dataclass(frozen=True)
class ScoreSummary:
    n: int
    mean: float
    standard_deviation: float
    minimum: float
    maximum: float


def prs_scores(genotypes: list[list[int]], weights: list[float]) -> list[float]:
    """Compute polygenic risk scores as genotype-weight dot products."""

    return [sum(g * w for g, w in zip(row, weights)) for row in genotypes]


def allele_counts(genotypes: list[list[int]]) -> list[int]:
    """Compute per-SNP alternate allele counts."""

    if not genotypes:
        return []
    n_snps = len(genotypes[0])
    counts = [0] * n_snps
    for row in genotypes:
        if len(row) != n_snps:
            raise ValueError("all genotype rows must have the same SNP count")
        for j, dosage in enumerate(row):
            counts[j] += dosage
    return counts


def case_control_allele_counts(
    genotypes: list[list[int]], phenotypes: list[int]
) -> tuple[list[int], list[int]]:
    """Compute per-SNP allele counts for controls and cases."""

    if len(genotypes) != len(phenotypes):
        raise ValueError("genotypes and phenotypes must have the same length")
    if not genotypes:
        return [], []

    n_snps = len(genotypes[0])
    control_counts = [0] * n_snps
    case_counts = [0] * n_snps

    for row, phenotype in zip(genotypes, phenotypes):
        target = case_counts if phenotype else control_counts
        for j, dosage in enumerate(row):
            target[j] += dosage

    return control_counts, case_counts


def summarize_scores(scores: list[float]) -> ScoreSummary:
    """Summarize a PRS vector for benchmark tables."""

    if not scores:
        raise ValueError("scores must not be empty")

    mean = sum(scores) / len(scores)
    variance = sum((score - mean) ** 2 for score in scores) / len(scores)
    return ScoreSummary(
        n=len(scores),
        mean=mean,
        standard_deviation=math.sqrt(variance),
        minimum=min(scores),
        maximum=max(scores),
    )

