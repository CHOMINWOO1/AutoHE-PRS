"""Synthetic genotype data for reproducible pilot experiments."""

from __future__ import annotations

from dataclasses import dataclass
import random


@dataclass(frozen=True)
class SyntheticGenomicDataset:
    """A small SNP dosage dataset.

    Genotypes are encoded as allele dosages: 0, 1, or 2.
    Rows are samples and columns are SNPs.
    """

    genotypes: list[list[int]]
    weights: list[float]
    phenotypes: list[int]
    allele_frequencies: list[float]

    @property
    def n_samples(self) -> int:
        return len(self.genotypes)

    @property
    def n_snps(self) -> int:
        return len(self.weights)


def _binomial_two(rng: random.Random, p: float) -> int:
    return int(rng.random() < p) + int(rng.random() < p)


def make_synthetic_dataset(
    n_samples: int = 128,
    n_snps: int = 2048,
    informative_fraction: float = 0.02,
    seed: int = 20260701,
) -> SyntheticGenomicDataset:
    """Create a reproducible synthetic genotype matrix and sparse PRS weights."""

    if n_samples <= 0:
        raise ValueError("n_samples must be positive")
    if n_snps <= 0:
        raise ValueError("n_snps must be positive")
    if not 0 < informative_fraction <= 1:
        raise ValueError("informative_fraction must be in (0, 1]")

    rng = random.Random(seed)
    allele_frequencies = [rng.uniform(0.05, 0.5) for _ in range(n_snps)]
    genotypes = [
        [_binomial_two(rng, allele_frequencies[j]) for j in range(n_snps)]
        for _ in range(n_samples)
    ]

    n_informative = max(1, int(n_snps * informative_fraction))
    informative_indices = set(rng.sample(range(n_snps), n_informative))
    weights = [
        rng.gauss(0.0, 0.08) if j in informative_indices else 0.0
        for j in range(n_snps)
    ]

    raw_scores = [sum(g * w for g, w in zip(row, weights)) for row in genotypes]
    threshold = sorted(raw_scores)[len(raw_scores) // 2]
    phenotypes = [int(score >= threshold) for score in raw_scores]

    return SyntheticGenomicDataset(
        genotypes=genotypes,
        weights=weights,
        phenotypes=phenotypes,
        allele_frequencies=allele_frequencies,
    )

