"""Slot-aware chunked CKKS benchmarks for larger PRS vectors."""

from __future__ import annotations

from dataclasses import dataclass
from time import perf_counter

from .tenseal_backend import (
    DEFAULT_COEFF_MOD_BIT_SIZES,
    DEFAULT_POLY_MODULUS_DEGREE,
    DEFAULT_SCALE,
    TenSEALUnavailable,
)


@dataclass(frozen=True)
class ChunkedPRSResult:
    scores: list[float]
    context_seconds: float
    encryption_seconds: float
    evaluation_seconds: float
    decryption_seconds: float
    serialization_seconds: float
    public_context_bytes: int
    input_ciphertext_bytes: int
    result_ciphertext_bytes: int
    window_snps: int
    windows_per_sample: int
    input_ciphertexts: int
    result_ciphertexts: int
    slots_per_ciphertext: int
    poly_modulus_degree: int
    scale: float
    coeff_mod_bit_sizes: tuple[int, ...]


def chunk_ranges(n_items: int, chunk_size: int) -> list[tuple[int, int]]:
    if n_items <= 0:
        raise ValueError("n_items must be positive")
    if chunk_size <= 0:
        raise ValueError("chunk_size must be positive")
    return [
        (start, min(start + chunk_size, n_items))
        for start in range(0, n_items, chunk_size)
    ]


def _validate_rectangular(genotypes: list[list[int]]) -> int:
    if not genotypes:
        raise ValueError("genotypes must not be empty")
    n_snps = len(genotypes[0])
    if n_snps == 0:
        raise ValueError("genotype rows must not be empty")
    for row in genotypes:
        if len(row) != n_snps:
            raise ValueError("all genotype rows must have the same SNP count")
    return n_snps


def chunked_prs_benchmark(
    genotypes: list[list[int]],
    weights: list[float],
    window_snps: int | None = None,
    poly_modulus_degree: int = DEFAULT_POLY_MODULUS_DEGREE,
    scale: float = DEFAULT_SCALE,
    coeff_mod_bit_sizes: tuple[int, ...] = DEFAULT_COEFF_MOD_BIT_SIZES,
) -> ChunkedPRSResult:
    """Run a per-sample CKKS PRS benchmark beyond a single ciphertext slot budget.

    Each sample is split into SNP windows. The server evaluates an encrypted dot
    product per window and homomorphically adds the encrypted partial scores for
    that sample. This keeps the one-sample privacy model while exposing the
    runtime and ciphertext cost of crossing the 4,096-slot boundary.
    """

    n_snps = _validate_rectangular(genotypes)
    if len(weights) != n_snps:
        raise ValueError("weights length must match genotype SNP count")

    slots_per_ciphertext = poly_modulus_degree // 2
    if window_snps is None:
        window_snps = slots_per_ciphertext
    if window_snps <= 0:
        raise ValueError("window_snps must be positive")
    if window_snps > slots_per_ciphertext:
        raise ValueError("window_snps exceeds available CKKS slots")

    try:
        import tenseal as ts
    except ImportError as exc:
        raise TenSEALUnavailable(
            "TenSEAL is not installed. Install requirements-he.txt in Python 3.10/3.11."
        ) from exc

    windows = chunk_ranges(n_snps, window_snps)

    start = perf_counter()
    context = ts.context(
        ts.SCHEME_TYPE.CKKS,
        poly_modulus_degree=poly_modulus_degree,
        coeff_mod_bit_sizes=list(coeff_mod_bit_sizes),
    )
    context.generate_galois_keys()
    context.global_scale = scale
    context_seconds = perf_counter() - start

    start = perf_counter()
    encrypted_sample_chunks = [
        [
            ts.ckks_vector(context, [float(value) for value in row[start:end]])
            for start, end in windows
        ]
        for row in genotypes
    ]
    encryption_seconds = perf_counter() - start

    weight_chunks = [
        [float(weight) for weight in weights[start:end]]
        for start, end in windows
    ]

    start = perf_counter()
    encrypted_scores = []
    for sample_chunks in encrypted_sample_chunks:
        encrypted_total = None
        for encrypted_chunk, weight_chunk in zip(sample_chunks, weight_chunks):
            encrypted_partial = encrypted_chunk.dot(weight_chunk)
            if encrypted_total is None:
                encrypted_total = encrypted_partial
            else:
                encrypted_total += encrypted_partial
        if encrypted_total is None:
            raise ValueError("no encrypted chunk was produced")
        encrypted_scores.append(encrypted_total)
    evaluation_seconds = perf_counter() - start

    start = perf_counter()
    scores = [float(encrypted_score.decrypt()[0]) for encrypted_score in encrypted_scores]
    decryption_seconds = perf_counter() - start

    start = perf_counter()
    public_context_bytes = len(context.serialize(save_secret_key=False))
    input_ciphertext_bytes = sum(
        len(encrypted_chunk.serialize())
        for sample_chunks in encrypted_sample_chunks
        for encrypted_chunk in sample_chunks
    )
    result_ciphertext_bytes = sum(
        len(encrypted_score.serialize()) for encrypted_score in encrypted_scores
    )
    serialization_seconds = perf_counter() - start

    return ChunkedPRSResult(
        scores=scores,
        context_seconds=context_seconds,
        encryption_seconds=encryption_seconds,
        evaluation_seconds=evaluation_seconds,
        decryption_seconds=decryption_seconds,
        serialization_seconds=serialization_seconds,
        public_context_bytes=public_context_bytes,
        input_ciphertext_bytes=input_ciphertext_bytes,
        result_ciphertext_bytes=result_ciphertext_bytes,
        window_snps=window_snps,
        windows_per_sample=len(windows),
        input_ciphertexts=len(genotypes) * len(windows),
        result_ciphertexts=len(encrypted_scores),
        slots_per_ciphertext=slots_per_ciphertext,
        poly_modulus_degree=poly_modulus_degree,
        scale=scale,
        coeff_mod_bit_sizes=coeff_mod_bit_sizes,
    )
