"""OpenFHE same-backend benchmarks for CKKS, BFV, and BGV PRS dot products."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from tempfile import TemporaryDirectory
from time import perf_counter

from .bfv_backend import integer_dot_bound, quantize_weights
from .chunked_ckks import chunk_ranges
from .openfhe_bgv_backend import _center_mod


DEFAULT_OPENFHE_RING_DIM = 16_384
DEFAULT_OPENFHE_WINDOW_SNPS = 4_096
DEFAULT_OPENFHE_PLAIN_MODULUS = 2_147_352_577
DEFAULT_OPENFHE_FIXED_POINT_SCALE = 10_000_000
DEFAULT_OPENFHE_CKKS_SCALING_MOD_SIZE = 50


class OpenFHEUnavailable(RuntimeError):
    """Raised when OpenFHE-Python is unavailable."""


@dataclass(frozen=True)
class OpenFHEControlledPRSResult:
    scheme: str
    scores: list[float]
    context_seconds: float
    encryption_seconds: float
    evaluation_seconds: float
    decryption_seconds: float
    serialization_seconds: float
    crypto_context_bytes: int
    public_key_bytes: int
    secret_key_bytes: int
    multiplication_key_bytes: int
    rotation_key_bytes: int
    total_key_material_bytes: int
    public_context_bytes: int
    input_ciphertext_bytes: int
    result_ciphertext_bytes: int
    window_snps: int
    windows_per_sample: int
    input_ciphertexts: int
    result_ciphertexts: int
    ring_dim: int
    batch_size: int
    plain_modulus: int | None
    fixed_point_scale: int | None
    ckks_scaling_mod_size: int | None
    ckks_first_modulus_size: int | None
    integer_weight_abs_sum: int | None
    integer_dot_abs_bound: int | None
    safe_modulus_margin: int | None
    aggregation_strategy: str
    reduction_strategy: str
    rotation_indices: tuple[int, ...]
    security_profile: str
    key_switching_technique: str
    scaling_technique: str | None
    key_generation_seconds: float
    multiplication_key_generation_seconds: float
    rotation_key_generation_seconds: float
    encoding_seconds: float
    decoding_seconds: float


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


def _load_openfhe():
    try:
        from openfhe import (  # type: ignore[import-not-found]
            BINARY,
            CCParamsBFVRNS,
            CCParamsBGVRNS,
            CCParamsCKKSRNS,
            GenCryptoContext,
            KeySwitchTechnique,
            PKESchemeFeature,
            ScalingTechnique,
            SecurityLevel,
            Serialize,
        )
    except ImportError as exc:
        raise OpenFHEUnavailable(
            "OpenFHE-Python is not importable. Run this benchmark in the "
            "OpenFHE Docker image."
        ) from exc
    return (
        BINARY,
        CCParamsBFVRNS,
        CCParamsBGVRNS,
        CCParamsCKKSRNS,
        GenCryptoContext,
        KeySwitchTechnique,
        PKESchemeFeature,
        ScalingTechnique,
        SecurityLevel,
        Serialize,
    )


def _enable_prs_features(context, feature_enum) -> None:
    for name in ("PKE", "KEYSWITCH", "LEVELEDSHE", "ADVANCEDSHE"):
        context.Enable(getattr(feature_enum, name))


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


def _serialize_key_material(context, key_pair, serialize, binary) -> tuple[int, int, int, int, int, int]:
    """Return independently measured context and key-material sizes.

    OpenFHE exposes public/secret keys as ordinary serializable objects, while
    multiplication and automorphism keys live in context-managed caches and
    must be serialized through the context API.
    """
    context_bytes = len(serialize(context, binary))
    public_key_bytes = len(serialize(key_pair.publicKey, binary))
    secret_key_bytes = len(serialize(key_pair.secretKey, binary))
    with TemporaryDirectory(prefix="autohe-prs-keys-") as raw:
        directory = Path(raw)
        multiplication_path = directory / "multiplication.keys"
        rotation_path = directory / "rotation.keys"
        if not context.SerializeEvalMultKey(str(multiplication_path), binary):
            raise RuntimeError("OpenFHE failed to serialize multiplication keys")
        if not context.SerializeEvalAutomorphismKey(str(rotation_path), binary):
            raise RuntimeError("OpenFHE failed to serialize rotation keys")
        multiplication_key_bytes = multiplication_path.stat().st_size
        rotation_key_bytes = rotation_path.stat().st_size
    total_key_material_bytes = (
        public_key_bytes
        + secret_key_bytes
        + multiplication_key_bytes
        + rotation_key_bytes
    )
    return (
        context_bytes,
        public_key_bytes,
        secret_key_bytes,
        multiplication_key_bytes,
        rotation_key_bytes,
        total_key_material_bytes,
    )


def _reduce_slots(context, ciphertext, active_slots: int, strategy: str):
    if strategy == "eval_sum":
        return context.EvalSum(ciphertext, active_slots)
    if strategy != "binary_tree":
        raise ValueError(f"unsupported reduction strategy: {strategy}")
    reduced = ciphertext
    rotation = 1
    while rotation < active_slots:
        reduced = context.EvalAdd(reduced, context.EvalRotate(reduced, rotation))
        rotation <<= 1
    return reduced


def _required_rotations(active_slots: int) -> tuple[int, ...]:
    rotations: list[int] = []
    value = 1
    while value < active_slots:
        rotations.append(value)
        value <<= 1
    return tuple(rotations)


def openfhe_ckks_prs_benchmark(
    genotypes: list[list[int]],
    weights: list[float],
    *,
    window_snps: int = DEFAULT_OPENFHE_WINDOW_SNPS,
    ring_dim: int = DEFAULT_OPENFHE_RING_DIM,
    scaling_mod_size: int = DEFAULT_OPENFHE_CKKS_SCALING_MOD_SIZE,
    first_modulus_size: int = 60,
    aggregation_strategy: str = "reduce_then_add",
    reduction_strategy: str = "binary_tree",
    rotation_indices: tuple[int, ...] | None = None,
    key_switching_technique: str = "HYBRID",
    scaling_technique: str = "FLEXIBLEAUTO",
) -> OpenFHEControlledPRSResult:
    n_snps = _validate_rectangular(genotypes, weights)
    if window_snps <= 0:
        raise ValueError("window_snps must be positive")
    if window_snps > ring_dim // 2:
        raise ValueError("window_snps exceeds CKKS slot capacity")
    if aggregation_strategy not in {"reduce_then_add", "add_then_reduce"}:
        raise ValueError("unsupported chunk aggregation strategy")
    if reduction_strategy not in {"binary_tree", "eval_sum"}:
        raise ValueError("unsupported reduction strategy")
    required_rotations = _required_rotations(min(n_snps, window_snps))
    rotation_indices = tuple(rotation_indices or required_rotations)
    missing_rotations = set(required_rotations) - set(rotation_indices)
    if missing_rotations:
        raise ValueError(f"missing rotation indices: {sorted(missing_rotations)}")
    windows = chunk_ranges(n_snps, window_snps)

    (
        BINARY,
        _CCParamsBFVRNS,
        _CCParamsBGVRNS,
        CCParamsCKKSRNS,
        GenCryptoContext,
        KeySwitchTechnique,
        PKESchemeFeature,
        ScalingTechnique,
        SecurityLevel,
        Serialize,
    ) = _load_openfhe()

    start = perf_counter()
    params = CCParamsCKKSRNS()
    params.SetMultiplicativeDepth(1)
    params.SetScalingModSize(scaling_mod_size)
    params.SetFirstModSize(first_modulus_size)
    params.SetRingDim(ring_dim)
    params.SetBatchSize(window_snps)
    params.SetSecurityLevel(SecurityLevel.HEStd_128_classic)
    try:
        params.SetKeySwitchTechnique(
            getattr(KeySwitchTechnique, key_switching_technique.upper())
        )
    except AttributeError as exc:
        raise ValueError(
            f"unsupported key-switching technique: {key_switching_technique}"
        ) from exc
    try:
        params.SetScalingTechnique(
            getattr(ScalingTechnique, scaling_technique.upper())
        )
    except AttributeError as exc:
        raise ValueError(f"unsupported CKKS scaling technique: {scaling_technique}") from exc
    context = GenCryptoContext(params)
    _enable_prs_features(context, PKESchemeFeature)
    context_seconds = perf_counter() - start

    start = perf_counter()
    key_pair = context.KeyGen()
    key_generation_seconds = perf_counter() - start

    start = perf_counter()
    context.EvalMultKeyGen(key_pair.secretKey)
    multiplication_key_generation_seconds = perf_counter() - start

    start = perf_counter()
    if reduction_strategy == "eval_sum":
        context.EvalSumKeyGen(key_pair.secretKey)
    else:
        context.EvalRotateKeyGen(key_pair.secretKey, list(rotation_indices))
    rotation_key_generation_seconds = perf_counter() - start

    start = perf_counter()
    weight_plaintexts = [
        context.MakeCKKSPackedPlaintext([float(value) for value in weights[start_index:end_index]])
        for start_index, end_index in windows
    ]
    encoding_seconds = perf_counter() - start

    start = perf_counter()
    encrypted_sample_chunks = [
        [
            context.Encrypt(
                key_pair.publicKey,
                context.MakeCKKSPackedPlaintext([float(value) for value in row[start_index:end_index]]),
            )
            for start_index, end_index in windows
        ]
        for row in genotypes
    ]
    encryption_seconds = perf_counter() - start

    start = perf_counter()
    encrypted_scores = []
    for sample_chunks in encrypted_sample_chunks:
        encrypted_total = None
        for encrypted_chunk, weight_plaintext, (start_index, end_index) in zip(
            sample_chunks, weight_plaintexts, windows
        ):
            encrypted_partial = context.EvalMult(encrypted_chunk, weight_plaintext)
            if aggregation_strategy == "reduce_then_add":
                encrypted_partial = _reduce_slots(
                    context,
                    encrypted_partial,
                    end_index - start_index,
                    reduction_strategy,
                )
            encrypted_total = (
                encrypted_partial
                if encrypted_total is None
                else context.EvalAdd(encrypted_total, encrypted_partial)
            )
        if encrypted_total is None:
            raise ValueError("no encrypted chunks were produced")
        if aggregation_strategy == "add_then_reduce":
            encrypted_total = _reduce_slots(
                context,
                encrypted_total,
                min(n_snps, window_snps),
                reduction_strategy,
            )
        encrypted_scores.append(encrypted_total)
    evaluation_seconds = perf_counter() - start

    start = perf_counter()
    decrypted_scores = [
        context.Decrypt(encrypted_score, key_pair.secretKey)
        for encrypted_score in encrypted_scores
    ]
    decryption_seconds = perf_counter() - start

    start = perf_counter()
    scores: list[float] = []
    for decrypted in decrypted_scores:
        decrypted.SetLength(1)
        scores.append(float(decrypted.GetCKKSPackedValue()[0].real))
    decoding_seconds = perf_counter() - start

    start = perf_counter()
    (
        crypto_context_bytes,
        public_key_bytes,
        secret_key_bytes,
        multiplication_key_bytes,
        rotation_key_bytes,
        total_key_material_bytes,
    ) = _serialize_key_material(context, key_pair, Serialize, BINARY)
    public_context_bytes = crypto_context_bytes + public_key_bytes
    input_ciphertext_bytes = sum(
        len(Serialize(encrypted_chunk, BINARY))
        for sample_chunks in encrypted_sample_chunks
        for encrypted_chunk in sample_chunks
    )
    result_ciphertext_bytes = sum(len(Serialize(encrypted_score, BINARY)) for encrypted_score in encrypted_scores)
    serialization_seconds = perf_counter() - start

    return OpenFHEControlledPRSResult(
        scheme="CKKS",
        scores=scores,
        context_seconds=context_seconds,
        encryption_seconds=encryption_seconds,
        evaluation_seconds=evaluation_seconds,
        decryption_seconds=decryption_seconds,
        serialization_seconds=serialization_seconds,
        crypto_context_bytes=crypto_context_bytes,
        public_key_bytes=public_key_bytes,
        secret_key_bytes=secret_key_bytes,
        multiplication_key_bytes=multiplication_key_bytes,
        rotation_key_bytes=rotation_key_bytes,
        total_key_material_bytes=total_key_material_bytes,
        public_context_bytes=public_context_bytes,
        input_ciphertext_bytes=input_ciphertext_bytes,
        result_ciphertext_bytes=result_ciphertext_bytes,
        window_snps=window_snps,
        windows_per_sample=len(windows),
        input_ciphertexts=len(genotypes) * len(windows),
        result_ciphertexts=len(encrypted_scores),
        ring_dim=ring_dim,
        batch_size=window_snps,
        plain_modulus=None,
        fixed_point_scale=None,
        ckks_scaling_mod_size=scaling_mod_size,
        ckks_first_modulus_size=first_modulus_size,
        integer_weight_abs_sum=None,
        integer_dot_abs_bound=None,
        safe_modulus_margin=None,
        aggregation_strategy=aggregation_strategy,
        reduction_strategy=reduction_strategy,
        rotation_indices=rotation_indices,
        security_profile="HEStd_128_classic",
        key_switching_technique=key_switching_technique.upper(),
        scaling_technique=scaling_technique.upper(),
        key_generation_seconds=key_generation_seconds,
        multiplication_key_generation_seconds=multiplication_key_generation_seconds,
        rotation_key_generation_seconds=rotation_key_generation_seconds,
        encoding_seconds=encoding_seconds,
        decoding_seconds=decoding_seconds,
    )


def openfhe_fixed_point_prs_benchmark(
    genotypes: list[list[int]],
    weights: list[float],
    *,
    scheme: str,
    window_snps: int = DEFAULT_OPENFHE_WINDOW_SNPS,
    ring_dim: int = DEFAULT_OPENFHE_RING_DIM,
    plain_modulus: int = DEFAULT_OPENFHE_PLAIN_MODULUS,
    fixed_point_scale: int = DEFAULT_OPENFHE_FIXED_POINT_SCALE,
    aggregation_strategy: str = "reduce_then_add",
    reduction_strategy: str = "binary_tree",
    rotation_indices: tuple[int, ...] | None = None,
    key_switching_technique: str = "HYBRID",
) -> OpenFHEControlledPRSResult:
    scheme = scheme.upper()
    if scheme not in {"BFV", "BGV"}:
        raise ValueError("scheme must be BFV or BGV")
    n_snps = _validate_rectangular(genotypes, weights)
    if window_snps <= 0:
        raise ValueError("window_snps must be positive")
    if window_snps > ring_dim:
        raise ValueError("window_snps exceeds OpenFHE packed plaintext capacity")
    if aggregation_strategy not in {"reduce_then_add", "add_then_reduce"}:
        raise ValueError("unsupported chunk aggregation strategy")
    if reduction_strategy not in {"binary_tree", "eval_sum"}:
        raise ValueError("unsupported reduction strategy")
    required_rotations = _required_rotations(min(n_snps, window_snps))
    rotation_indices = tuple(rotation_indices or required_rotations)
    missing_rotations = set(required_rotations) - set(rotation_indices)
    if missing_rotations:
        raise ValueError(f"missing rotation indices: {sorted(missing_rotations)}")
    windows = chunk_ranges(n_snps, window_snps)
    integer_weights, integer_abs_sum, dot_bound, safe_margin = _fixed_point_metadata(
        genotypes,
        weights,
        fixed_point_scale,
        plain_modulus,
    )

    (
        BINARY,
        CCParamsBFVRNS,
        CCParamsBGVRNS,
        _CCParamsCKKSRNS,
        GenCryptoContext,
        KeySwitchTechnique,
        PKESchemeFeature,
        _ScalingTechnique,
        SecurityLevel,
        Serialize,
    ) = _load_openfhe()

    start = perf_counter()
    params = CCParamsBFVRNS() if scheme == "BFV" else CCParamsBGVRNS()
    params.SetPlaintextModulus(plain_modulus)
    params.SetMultiplicativeDepth(1)
    params.SetRingDim(ring_dim)
    params.SetBatchSize(window_snps)
    params.SetSecurityLevel(SecurityLevel.HEStd_128_classic)
    try:
        params.SetKeySwitchTechnique(
            getattr(KeySwitchTechnique, key_switching_technique.upper())
        )
    except AttributeError as exc:
        raise ValueError(
            f"unsupported key-switching technique: {key_switching_technique}"
        ) from exc
    context = GenCryptoContext(params)
    _enable_prs_features(context, PKESchemeFeature)
    context_seconds = perf_counter() - start

    start = perf_counter()
    key_pair = context.KeyGen()
    key_generation_seconds = perf_counter() - start

    start = perf_counter()
    context.EvalMultKeyGen(key_pair.secretKey)
    multiplication_key_generation_seconds = perf_counter() - start

    start = perf_counter()
    if reduction_strategy == "eval_sum":
        context.EvalSumKeyGen(key_pair.secretKey)
    else:
        context.EvalRotateKeyGen(key_pair.secretKey, list(rotation_indices))
    rotation_key_generation_seconds = perf_counter() - start

    start = perf_counter()
    weight_plaintexts = [
        context.MakePackedPlaintext(integer_weights[start_index:end_index])
        for start_index, end_index in windows
    ]
    encoding_seconds = perf_counter() - start

    start = perf_counter()
    encrypted_sample_chunks = [
        [
            context.Encrypt(key_pair.publicKey, context.MakePackedPlaintext(row[start_index:end_index]))
            for start_index, end_index in windows
        ]
        for row in genotypes
    ]
    encryption_seconds = perf_counter() - start

    start = perf_counter()
    encrypted_scores = []
    for sample_chunks in encrypted_sample_chunks:
        encrypted_total = None
        for encrypted_chunk, weight_plaintext, (start_index, end_index) in zip(
            sample_chunks, weight_plaintexts, windows
        ):
            encrypted_partial = context.EvalMult(encrypted_chunk, weight_plaintext)
            if aggregation_strategy == "reduce_then_add":
                encrypted_partial = _reduce_slots(
                    context,
                    encrypted_partial,
                    end_index - start_index,
                    reduction_strategy,
                )
            encrypted_total = (
                encrypted_partial
                if encrypted_total is None
                else context.EvalAdd(encrypted_total, encrypted_partial)
            )
        if encrypted_total is None:
            raise ValueError("no encrypted chunks were produced")
        if aggregation_strategy == "add_then_reduce":
            encrypted_total = _reduce_slots(
                context,
                encrypted_total,
                min(n_snps, window_snps),
                reduction_strategy,
            )
        encrypted_scores.append(encrypted_total)
    evaluation_seconds = perf_counter() - start

    start = perf_counter()
    decrypted_scores = [
        context.Decrypt(encrypted_score, key_pair.secretKey)
        for encrypted_score in encrypted_scores
    ]
    decryption_seconds = perf_counter() - start

    start = perf_counter()
    integer_scores: list[int] = []
    for decrypted in decrypted_scores:
        decrypted.SetLength(1)
        raw_value = int(decrypted.GetPackedValue()[0])
        integer_scores.append(_center_mod(raw_value, plain_modulus))
    scores = [integer_score / fixed_point_scale for integer_score in integer_scores]
    decoding_seconds = perf_counter() - start

    start = perf_counter()
    (
        crypto_context_bytes,
        public_key_bytes,
        secret_key_bytes,
        multiplication_key_bytes,
        rotation_key_bytes,
        total_key_material_bytes,
    ) = _serialize_key_material(context, key_pair, Serialize, BINARY)
    public_context_bytes = crypto_context_bytes + public_key_bytes
    input_ciphertext_bytes = sum(
        len(Serialize(encrypted_chunk, BINARY))
        for sample_chunks in encrypted_sample_chunks
        for encrypted_chunk in sample_chunks
    )
    result_ciphertext_bytes = sum(len(Serialize(encrypted_score, BINARY)) for encrypted_score in encrypted_scores)
    serialization_seconds = perf_counter() - start

    return OpenFHEControlledPRSResult(
        scheme=scheme,
        scores=scores,
        context_seconds=context_seconds,
        encryption_seconds=encryption_seconds,
        evaluation_seconds=evaluation_seconds,
        decryption_seconds=decryption_seconds,
        serialization_seconds=serialization_seconds,
        crypto_context_bytes=crypto_context_bytes,
        public_key_bytes=public_key_bytes,
        secret_key_bytes=secret_key_bytes,
        multiplication_key_bytes=multiplication_key_bytes,
        rotation_key_bytes=rotation_key_bytes,
        total_key_material_bytes=total_key_material_bytes,
        public_context_bytes=public_context_bytes,
        input_ciphertext_bytes=input_ciphertext_bytes,
        result_ciphertext_bytes=result_ciphertext_bytes,
        window_snps=window_snps,
        windows_per_sample=len(windows),
        input_ciphertexts=len(genotypes) * len(windows),
        result_ciphertexts=len(encrypted_scores),
        ring_dim=ring_dim,
        batch_size=window_snps,
        plain_modulus=plain_modulus,
        fixed_point_scale=fixed_point_scale,
        ckks_scaling_mod_size=None,
        ckks_first_modulus_size=None,
        integer_weight_abs_sum=integer_abs_sum,
        integer_dot_abs_bound=dot_bound,
        safe_modulus_margin=safe_margin,
        aggregation_strategy=aggregation_strategy,
        reduction_strategy=reduction_strategy,
        rotation_indices=rotation_indices,
        security_profile="HEStd_128_classic",
        key_switching_technique=key_switching_technique.upper(),
        scaling_technique=None,
        key_generation_seconds=key_generation_seconds,
        multiplication_key_generation_seconds=multiplication_key_generation_seconds,
        rotation_key_generation_seconds=rotation_key_generation_seconds,
        encoding_seconds=encoding_seconds,
        decoding_seconds=decoding_seconds,
    )
