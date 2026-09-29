"""OpenFHE BGV fixed-point benchmarks for PRS dot products."""

from __future__ import annotations

from dataclasses import dataclass
from time import perf_counter

from .bfv_backend import (
    DEFAULT_BFV_FIXED_POINT_SCALE,
    integer_dot_bound,
    quantize_weights,
)


DEFAULT_BGV_POLY_MODULUS_DEGREE = 16384
DEFAULT_BGV_PLAIN_MODULUS = 2_147_352_577
DEFAULT_BGV_FIXED_POINT_SCALE = DEFAULT_BFV_FIXED_POINT_SCALE


class OpenFHEUnavailable(RuntimeError):
    """Raised when OpenFHE-Python is not available in the current runtime."""


@dataclass(frozen=True)
class OpenFHEBGVPRSResult:
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
    fixed_point_scale: int
    integer_weight_abs_sum: int
    integer_dot_abs_bound: int
    safe_modulus_margin: int


def _center_mod(value: int, modulus: int) -> int:
    half = modulus // 2
    return value - modulus if value > half else value


def openfhe_bgv_fixed_point_prs_benchmark(
    genotypes: list[list[int]],
    weights: list[float],
    *,
    poly_modulus_degree: int = DEFAULT_BGV_POLY_MODULUS_DEGREE,
    plain_modulus: int = DEFAULT_BGV_PLAIN_MODULUS,
    fixed_point_scale: int = DEFAULT_BGV_FIXED_POINT_SCALE,
) -> OpenFHEBGVPRSResult:
    """Run a BGV fixed-point PRS benchmark with OpenFHE-Python.

    BGV evaluates integer arithmetic. Real-valued PRS weights are quantized as
    round(weight * fixed_point_scale), then decrypted integer scores are divided
    by the same scale.
    """

    try:
        from openfhe import (  # type: ignore[import-not-found]
            BINARY,
            CCParamsBGVRNS,
            GenCryptoContext,
            PKESchemeFeature,
            Serialize,
        )
    except ImportError as exc:
        raise OpenFHEUnavailable(
            "OpenFHE-Python is not importable. Run this backend in a Linux "
            "Python environment with openfhe installed."
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
    if n_snps > poly_modulus_degree:
        raise ValueError("n_snps exceeds configured BGV batch/ring capacity")

    integer_weights = quantize_weights(weights, fixed_point_scale)
    integer_abs_sum = sum(abs(weight) for weight in integer_weights)
    dot_bound = integer_dot_bound(genotypes, integer_weights)
    safe_margin = plain_modulus // 2 - dot_bound
    if safe_margin <= 0:
        raise ValueError(
            "BGV fixed-point dot product may wrap modulo the plaintext modulus: "
            f"bound={dot_bound}, plain_modulus={plain_modulus}, scale={fixed_point_scale}"
        )

    start = perf_counter()
    params = CCParamsBGVRNS()
    params.SetPlaintextModulus(plain_modulus)
    params.SetMultiplicativeDepth(1)
    params.SetRingDim(poly_modulus_degree)
    params.SetBatchSize(n_snps)
    context = GenCryptoContext(params)
    context.Enable(PKESchemeFeature.PKE)
    context.Enable(PKESchemeFeature.KEYSWITCH)
    context.Enable(PKESchemeFeature.LEVELEDSHE)
    context.Enable(PKESchemeFeature.ADVANCEDSHE)
    key_pair = context.KeyGen()
    context.EvalMultKeyGen(key_pair.secretKey)
    context.EvalSumKeyGen(key_pair.secretKey)
    weight_plaintext = context.MakePackedPlaintext(integer_weights)
    context_seconds = perf_counter() - start

    start = perf_counter()
    encrypted_rows = [
        context.Encrypt(key_pair.publicKey, context.MakePackedPlaintext(row))
        for row in genotypes
    ]
    encryption_seconds = perf_counter() - start

    start = perf_counter()
    encrypted_scores = [
        context.EvalInnerProduct(encrypted_row, weight_plaintext, n_snps)
        for encrypted_row in encrypted_rows
    ]
    evaluation_seconds = perf_counter() - start

    start = perf_counter()
    integer_scores: list[int] = []
    for encrypted_score in encrypted_scores:
        decrypted = context.Decrypt(encrypted_score, key_pair.secretKey)
        decrypted.SetLength(1)
        raw_value = int(decrypted.GetPackedValue()[0])
        integer_scores.append(_center_mod(raw_value, plain_modulus))
    scores = [integer_score / fixed_point_scale for integer_score in integer_scores]
    decryption_seconds = perf_counter() - start

    start = perf_counter()
    public_context_bytes = len(Serialize(context, BINARY)) + len(Serialize(key_pair.publicKey, BINARY))
    input_ciphertext_bytes = sum(len(Serialize(encrypted_row, BINARY)) for encrypted_row in encrypted_rows)
    result_ciphertext_bytes = sum(len(Serialize(encrypted_score, BINARY)) for encrypted_score in encrypted_scores)
    serialization_seconds = perf_counter() - start

    return OpenFHEBGVPRSResult(
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
        fixed_point_scale=fixed_point_scale,
        integer_weight_abs_sum=integer_abs_sum,
        integer_dot_abs_bound=dot_bound,
        safe_modulus_margin=safe_margin,
    )
