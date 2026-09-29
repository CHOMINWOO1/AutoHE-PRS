"""TenSEAL BFV fixed-point benchmarks for PRS dot products."""

from __future__ import annotations

from dataclasses import dataclass
from time import perf_counter

from .tenseal_backend import TenSEALUnavailable


DEFAULT_BFV_POLY_MODULUS_DEGREE = 8192
DEFAULT_BFV_PLAIN_MODULUS_BITS = 31
DEFAULT_BFV_FIXED_POINT_SCALE = 10_000_000


@dataclass(frozen=True)
class BFVPRSResult:
    scores: list[float]
    integer_scores: list[int]
    context_seconds: float
    encryption_seconds: float
    evaluation_seconds: float
    decryption_seconds: float
    serialization_seconds: float
    public_context_bytes: int
    input_ciphertext_bytes: int
    result_ciphertext_bytes: int
    poly_modulus_degree: int
    plain_modulus: int
    plain_modulus_bits: int
    fixed_point_scale: int
    integer_weight_abs_sum: int
    integer_dot_abs_bound: int
    safe_modulus_margin: int


def quantize_weights(weights: list[float], fixed_point_scale: int) -> list[int]:
    return [int(round(weight * fixed_point_scale)) for weight in weights]


def integer_dot_bound(genotypes: list[list[int]], integer_weights: list[int]) -> int:
    if not genotypes:
        return 0
    max_abs_dosage = max(abs(value) for row in genotypes for value in row)
    return max_abs_dosage * sum(abs(weight) for weight in integer_weights)


def bfv_fixed_point_prs_benchmark(
    genotypes: list[list[int]],
    weights: list[float],
    *,
    poly_modulus_degree: int = DEFAULT_BFV_POLY_MODULUS_DEGREE,
    plain_modulus_bits: int = DEFAULT_BFV_PLAIN_MODULUS_BITS,
    fixed_point_scale: int = DEFAULT_BFV_FIXED_POINT_SCALE,
) -> BFVPRSResult:
    """Run a BFV fixed-point PRS benchmark.

    BFV operates over integers. To compare it with CKKS on real-valued PRS
    weights, weights are quantized as round(weight * fixed_point_scale). The
    decrypted integer dot product is divided by the same scale.
    """

    try:
        import tenseal as ts
        from tenseal.sealapi import PlainModulus
    except ImportError as exc:
        raise TenSEALUnavailable(
            "TenSEAL is not installed. Install requirements-he.txt in Python 3.10/3.11."
        ) from exc

    if not genotypes:
        raise ValueError("genotypes must not be empty")
    n_snps = len(genotypes[0])
    if len(weights) != n_snps:
        raise ValueError("weights length must match genotype width")
    if any(len(row) != n_snps for row in genotypes):
        raise ValueError("all genotype rows must have the same SNP count")
    if fixed_point_scale <= 0:
        raise ValueError("fixed_point_scale must be positive")

    integer_weights = quantize_weights(weights, fixed_point_scale)
    integer_abs_sum = sum(abs(weight) for weight in integer_weights)
    dot_bound = integer_dot_bound(genotypes, integer_weights)

    start = perf_counter()
    plain_modulus = PlainModulus.Batching(poly_modulus_degree, plain_modulus_bits).value()
    safe_margin = plain_modulus // 2 - dot_bound
    if safe_margin <= 0:
        raise ValueError(
            "BFV fixed-point dot product may wrap modulo the plaintext modulus: "
            f"bound={dot_bound}, plain_modulus={plain_modulus}, scale={fixed_point_scale}"
        )
    context = ts.context(
        ts.SCHEME_TYPE.BFV,
        poly_modulus_degree=poly_modulus_degree,
        plain_modulus=plain_modulus,
    )
    context.generate_galois_keys()
    context_seconds = perf_counter() - start

    start = perf_counter()
    encrypted_rows = [ts.bfv_vector(context, row) for row in genotypes]
    encryption_seconds = perf_counter() - start

    start = perf_counter()
    encrypted_scores = [encrypted_row.dot(integer_weights) for encrypted_row in encrypted_rows]
    evaluation_seconds = perf_counter() - start

    start = perf_counter()
    integer_scores = [int(encrypted_score.decrypt()[0]) for encrypted_score in encrypted_scores]
    scores = [integer_score / fixed_point_scale for integer_score in integer_scores]
    decryption_seconds = perf_counter() - start

    start = perf_counter()
    public_context_bytes = len(context.serialize(save_secret_key=False))
    input_ciphertext_bytes = sum(len(encrypted_row.serialize()) for encrypted_row in encrypted_rows)
    result_ciphertext_bytes = sum(len(encrypted_score.serialize()) for encrypted_score in encrypted_scores)
    serialization_seconds = perf_counter() - start

    return BFVPRSResult(
        scores=scores,
        integer_scores=integer_scores,
        context_seconds=context_seconds,
        encryption_seconds=encryption_seconds,
        evaluation_seconds=evaluation_seconds,
        decryption_seconds=decryption_seconds,
        serialization_seconds=serialization_seconds,
        public_context_bytes=public_context_bytes,
        input_ciphertext_bytes=input_ciphertext_bytes,
        result_ciphertext_bytes=result_ciphertext_bytes,
        poly_modulus_degree=poly_modulus_degree,
        plain_modulus=plain_modulus,
        plain_modulus_bits=plain_modulus_bits,
        fixed_point_scale=fixed_point_scale,
        integer_weight_abs_sum=integer_abs_sum,
        integer_dot_abs_bound=dot_bound,
        safe_modulus_margin=safe_margin,
    )
