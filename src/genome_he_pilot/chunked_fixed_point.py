"""Chunked fixed-point HE benchmarks for larger PRS vectors."""

from __future__ import annotations

from dataclasses import dataclass
from time import perf_counter

from .bfv_backend import (
    DEFAULT_BFV_FIXED_POINT_SCALE,
    DEFAULT_BFV_PLAIN_MODULUS_BITS,
    DEFAULT_BFV_POLY_MODULUS_DEGREE,
    integer_dot_bound,
    quantize_weights,
)
from .chunked_ckks import chunk_ranges
from .openfhe_bgv_backend import (
    DEFAULT_BGV_FIXED_POINT_SCALE,
    DEFAULT_BGV_PLAIN_MODULUS,
    DEFAULT_BGV_POLY_MODULUS_DEGREE,
    OpenFHEUnavailable,
    _center_mod,
)
from .tenseal_backend import TenSEALUnavailable


DEFAULT_FIXED_POINT_WINDOW_SNPS = 4096


@dataclass(frozen=True)
class ChunkedFixedPointPRSResult:
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
    window_snps: int
    windows_per_sample: int
    input_ciphertexts: int
    result_ciphertexts: int
    poly_modulus_degree: int
    plain_modulus: int
    fixed_point_scale: int
    integer_weight_abs_sum: int
    integer_dot_abs_bound: int
    safe_modulus_margin: int


def _validate_rectangular(genotypes: list[list[int]], weights: list[float]) -> int:
    if not genotypes:
        raise ValueError("genotypes must not be empty")
    n_snps = len(genotypes[0])
    if n_snps == 0:
        raise ValueError("genotype rows must not be empty")
    if len(weights) != n_snps:
        raise ValueError("weights length must match genotype SNP count")
    if any(len(row) != n_snps for row in genotypes):
        raise ValueError("all genotype rows must have the same SNP count")
    return n_snps


def _fixed_point_metadata(
    genotypes: list[list[int]],
    weights: list[float],
    fixed_point_scale: int,
    plain_modulus: int,
) -> tuple[list[int], int, int, int]:
    if fixed_point_scale <= 0:
        raise ValueError("fixed_point_scale must be positive")
    integer_weights = quantize_weights(weights, fixed_point_scale)
    integer_abs_sum = sum(abs(weight) for weight in integer_weights)
    dot_bound = integer_dot_bound(genotypes, integer_weights)
    safe_margin = plain_modulus // 2 - dot_bound
    if safe_margin <= 0:
        raise ValueError(
            "fixed-point dot product may wrap modulo the plaintext modulus: "
            f"bound={dot_bound}, plain_modulus={plain_modulus}, scale={fixed_point_scale}"
        )
    return integer_weights, integer_abs_sum, dot_bound, safe_margin


def chunked_bfv_fixed_point_prs_benchmark(
    genotypes: list[list[int]],
    weights: list[float],
    *,
    window_snps: int = DEFAULT_FIXED_POINT_WINDOW_SNPS,
    poly_modulus_degree: int = DEFAULT_BFV_POLY_MODULUS_DEGREE,
    plain_modulus_bits: int = DEFAULT_BFV_PLAIN_MODULUS_BITS,
    fixed_point_scale: int = DEFAULT_BFV_FIXED_POINT_SCALE,
) -> ChunkedFixedPointPRSResult:
    """Run a chunked BFV fixed-point PRS benchmark."""

    try:
        import tenseal as ts
        from tenseal.sealapi import PlainModulus
    except ImportError as exc:
        raise TenSEALUnavailable(
            "TenSEAL is not installed. Install requirements-he.txt in Python 3.10/3.11."
        ) from exc

    n_snps = _validate_rectangular(genotypes, weights)
    if window_snps <= 0:
        raise ValueError("window_snps must be positive")
    if window_snps > poly_modulus_degree:
        raise ValueError("window_snps exceeds configured BFV vector capacity")

    plain_modulus = PlainModulus.Batching(poly_modulus_degree, plain_modulus_bits).value()
    integer_weights, integer_abs_sum, dot_bound, safe_margin = _fixed_point_metadata(
        genotypes, weights, fixed_point_scale, plain_modulus
    )
    windows = chunk_ranges(n_snps, window_snps)

    start = perf_counter()
    context = ts.context(
        ts.SCHEME_TYPE.BFV,
        poly_modulus_degree=poly_modulus_degree,
        plain_modulus=plain_modulus,
    )
    context.generate_galois_keys()
    context_seconds = perf_counter() - start

    start = perf_counter()
    encrypted_sample_chunks = [
        [ts.bfv_vector(context, row[start:end]) for start, end in windows]
        for row in genotypes
    ]
    encryption_seconds = perf_counter() - start

    weight_chunks = [integer_weights[start:end] for start, end in windows]

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
    integer_scores = [int(encrypted_score.decrypt()[0]) for encrypted_score in encrypted_scores]
    scores = [integer_score / fixed_point_scale for integer_score in integer_scores]
    decryption_seconds = perf_counter() - start

    start = perf_counter()
    public_context_bytes = len(context.serialize(save_secret_key=False))
    input_ciphertext_bytes = sum(
        len(encrypted_chunk.serialize())
        for sample_chunks in encrypted_sample_chunks
        for encrypted_chunk in sample_chunks
    )
    result_ciphertext_bytes = sum(len(encrypted_score.serialize()) for encrypted_score in encrypted_scores)
    serialization_seconds = perf_counter() - start

    return ChunkedFixedPointPRSResult(
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
        window_snps=window_snps,
        windows_per_sample=len(windows),
        input_ciphertexts=len(genotypes) * len(windows),
        result_ciphertexts=len(encrypted_scores),
        poly_modulus_degree=poly_modulus_degree,
        plain_modulus=plain_modulus,
        fixed_point_scale=fixed_point_scale,
        integer_weight_abs_sum=integer_abs_sum,
        integer_dot_abs_bound=dot_bound,
        safe_modulus_margin=safe_margin,
    )


def chunked_openfhe_bgv_fixed_point_prs_benchmark(
    genotypes: list[list[int]],
    weights: list[float],
    *,
    window_snps: int = DEFAULT_FIXED_POINT_WINDOW_SNPS,
    poly_modulus_degree: int = DEFAULT_BGV_POLY_MODULUS_DEGREE,
    plain_modulus: int = DEFAULT_BGV_PLAIN_MODULUS,
    fixed_point_scale: int = DEFAULT_BGV_FIXED_POINT_SCALE,
) -> ChunkedFixedPointPRSResult:
    """Run a chunked OpenFHE BGV fixed-point PRS benchmark."""

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

    n_snps = _validate_rectangular(genotypes, weights)
    if window_snps <= 0:
        raise ValueError("window_snps must be positive")
    if window_snps > poly_modulus_degree:
        raise ValueError("window_snps exceeds configured BGV batch/ring capacity")

    integer_weights, integer_abs_sum, dot_bound, safe_margin = _fixed_point_metadata(
        genotypes, weights, fixed_point_scale, plain_modulus
    )
    windows = chunk_ranges(n_snps, window_snps)

    start = perf_counter()
    params = CCParamsBGVRNS()
    params.SetPlaintextModulus(plain_modulus)
    params.SetMultiplicativeDepth(1)
    params.SetRingDim(poly_modulus_degree)
    params.SetBatchSize(window_snps)
    context = GenCryptoContext(params)
    context.Enable(PKESchemeFeature.PKE)
    context.Enable(PKESchemeFeature.KEYSWITCH)
    context.Enable(PKESchemeFeature.LEVELEDSHE)
    context.Enable(PKESchemeFeature.ADVANCEDSHE)
    key_pair = context.KeyGen()
    context.EvalMultKeyGen(key_pair.secretKey)
    context.EvalSumKeyGen(key_pair.secretKey)
    weight_plaintexts = [
        context.MakePackedPlaintext(integer_weights[start:end])
        for start, end in windows
    ]
    context_seconds = perf_counter() - start

    start = perf_counter()
    encrypted_sample_chunks = [
        [
            context.Encrypt(key_pair.publicKey, context.MakePackedPlaintext(row[start:end]))
            for start, end in windows
        ]
        for row in genotypes
    ]
    encryption_seconds = perf_counter() - start

    start = perf_counter()
    encrypted_scores = []
    for sample_chunks in encrypted_sample_chunks:
        encrypted_total = None
        for encrypted_chunk, weight_plaintext, (start_index, end_index) in zip(sample_chunks, weight_plaintexts, windows):
            encrypted_partial = context.EvalInnerProduct(
                encrypted_chunk,
                weight_plaintext,
                end_index - start_index,
            )
            if encrypted_total is None:
                encrypted_total = encrypted_partial
            else:
                encrypted_total = context.EvalAdd(encrypted_total, encrypted_partial)
        if encrypted_total is None:
            raise ValueError("no encrypted chunk was produced")
        encrypted_scores.append(encrypted_total)
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
    input_ciphertext_bytes = sum(
        len(Serialize(encrypted_chunk, BINARY))
        for sample_chunks in encrypted_sample_chunks
        for encrypted_chunk in sample_chunks
    )
    result_ciphertext_bytes = sum(len(Serialize(encrypted_score, BINARY)) for encrypted_score in encrypted_scores)
    serialization_seconds = perf_counter() - start

    return ChunkedFixedPointPRSResult(
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
        window_snps=window_snps,
        windows_per_sample=len(windows),
        input_ciphertexts=len(genotypes) * len(windows),
        result_ciphertexts=len(encrypted_scores),
        poly_modulus_degree=poly_modulus_degree,
        plain_modulus=plain_modulus,
        fixed_point_scale=fixed_point_scale,
        integer_weight_abs_sum=integer_abs_sum,
        integer_dot_abs_bound=dot_bound,
        safe_modulus_margin=safe_margin,
    )
