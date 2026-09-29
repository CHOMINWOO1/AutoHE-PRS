"""Optional TenSEAL backend for encrypted PRS experiments."""

from __future__ import annotations

from dataclasses import dataclass
from time import perf_counter


class TenSEALUnavailable(RuntimeError):
    """Raised when TenSEAL is requested but not installed."""


DEFAULT_POLY_MODULUS_DEGREE = 8192
DEFAULT_SCALE = 2**40
DEFAULT_COEFF_MOD_BIT_SIZES = (60, 40, 40, 60)


@dataclass(frozen=True)
class EncryptedPRSResult:
    scores: list[float]
    context_seconds: float
    encryption_seconds: float
    evaluation_seconds: float
    decryption_seconds: float
    serialization_seconds: float
    public_context_bytes: int
    input_ciphertext_bytes: int
    result_ciphertext_bytes: int
    poly_modulus_degree: int
    scale: float
    coeff_mod_bit_sizes: tuple[int, ...]


def encrypted_prs_benchmark(
    genotypes: list[list[int]],
    weights: list[float],
    poly_modulus_degree: int = DEFAULT_POLY_MODULUS_DEGREE,
    scale: float = DEFAULT_SCALE,
    coeff_mod_bit_sizes: tuple[int, ...] = DEFAULT_COEFF_MOD_BIT_SIZES,
) -> EncryptedPRSResult:
    """Run a measured CKKS encrypted PRS benchmark.

    The implementation encrypts one sample vector per ciphertext. This is easy
    to explain for a KCI pilot and leaves packed multi-sample optimization for
    the SCI extension.
    """

    try:
        import tenseal as ts
    except ImportError as exc:
        raise TenSEALUnavailable(
            "TenSEAL is not installed. Install requirements-he.txt in Python 3.10/3.11."
        ) from exc

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
    encrypted_rows = [
        ts.ckks_vector(context, [float(value) for value in row]) for row in genotypes
    ]
    encryption_seconds = perf_counter() - start

    start = perf_counter()
    encrypted_scores = [encrypted_row.dot(weights) for encrypted_row in encrypted_rows]
    evaluation_seconds = perf_counter() - start

    start = perf_counter()
    scores = [float(encrypted_score.decrypt()[0]) for encrypted_score in encrypted_scores]
    decryption_seconds = perf_counter() - start

    start = perf_counter()
    public_context_bytes = len(context.serialize(save_secret_key=False))
    input_ciphertext_bytes = sum(len(encrypted_row.serialize()) for encrypted_row in encrypted_rows)
    result_ciphertext_bytes = sum(
        len(encrypted_score.serialize()) for encrypted_score in encrypted_scores
    )
    serialization_seconds = perf_counter() - start

    return EncryptedPRSResult(
        scores=scores,
        context_seconds=context_seconds,
        encryption_seconds=encryption_seconds,
        evaluation_seconds=evaluation_seconds,
        decryption_seconds=decryption_seconds,
        serialization_seconds=serialization_seconds,
        public_context_bytes=public_context_bytes,
        input_ciphertext_bytes=input_ciphertext_bytes,
        result_ciphertext_bytes=result_ciphertext_bytes,
        poly_modulus_degree=poly_modulus_degree,
        scale=scale,
        coeff_mod_bit_sizes=coeff_mod_bit_sizes,
    )


def encrypted_prs_scores(
    genotypes: list[list[int]],
    weights: list[float],
    poly_modulus_degree: int = DEFAULT_POLY_MODULUS_DEGREE,
    scale: float = DEFAULT_SCALE,
) -> list[float]:
    """Compute PRS scores with CKKS-encrypted genotype vectors.

    This simple backend encrypts one sample vector at a time. It is suitable for
    the KCI pilot validation, not yet the optimized SCI-scale implementation.
    """

    return encrypted_prs_benchmark(
        genotypes=genotypes,
        weights=weights,
        poly_modulus_degree=poly_modulus_degree,
        scale=scale,
    ).scores
