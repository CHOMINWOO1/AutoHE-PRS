"""Packed CKKS prototypes for SCI extension experiments."""

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
class PackedPRSResult:
    scores: list[float]
    context_seconds: float
    encryption_seconds: float
    evaluation_seconds: float
    decryption_seconds: float
    serialization_seconds: float
    public_context_bytes: int
    input_ciphertext_bytes: int
    result_ciphertext_bytes: int
    block_size: int
    blocks: int
    slots_per_ciphertext: int
    poly_modulus_degree: int
    scale: float
    coeff_mod_bit_sizes: tuple[int, ...]


@dataclass(frozen=True)
class PackedAggregateResult:
    block_sums: list[float]
    context_seconds: float
    encryption_seconds: float
    evaluation_seconds: float
    decryption_seconds: float
    serialization_seconds: float
    public_context_bytes: int
    input_ciphertext_bytes: int
    result_ciphertext_bytes: int
    block_size: int
    blocks: int
    slots_per_ciphertext: int
    poly_modulus_degree: int
    scale: float
    coeff_mod_bit_sizes: tuple[int, ...]


@dataclass(frozen=True)
class WindowedPackedAggregateResult:
    site_sums: list[float]
    context_seconds: float
    encryption_seconds: float
    evaluation_seconds: float
    decryption_seconds: float
    serialization_seconds: float
    public_context_bytes: int
    input_ciphertext_bytes: int
    result_ciphertext_bytes: int
    block_size: int
    blocks: int
    windows_per_site: int
    ciphertexts: int
    calls: int
    slots_per_ciphertext: int
    window_snps: int
    poly_modulus_degree: int
    scale: float
    coeff_mod_bit_sizes: tuple[int, ...]


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


def _block_weight_matrix(weights: list[float], block_size: int) -> list[list[float]]:
    matrix: list[list[float]] = []
    for sample_index in range(block_size):
        for weight in weights:
            row = [0.0] * block_size
            row[sample_index] = float(weight)
            matrix.append(row)
    return matrix


def _block_weight_vector(weights: list[float], block_size: int) -> list[float]:
    values: list[float] = []
    for _sample_index in range(block_size):
        values.extend(float(weight) for weight in weights)
    return values


def _flatten_block(block: list[list[int]], block_size: int, n_snps: int) -> list[float]:
    values: list[float] = []
    for row in block:
        values.extend(float(value) for value in row)
    missing_rows = block_size - len(block)
    if missing_rows > 0:
        values.extend([0.0] * missing_rows * n_snps)
    return values


def aggregate_window_ranges(n_snps: int, window_snps: int) -> list[tuple[int, int]]:
    if n_snps <= 0:
        raise ValueError("n_snps must be positive")
    if window_snps <= 0:
        raise ValueError("window_snps must be positive")
    return [
        (start, min(start + window_snps, n_snps))
        for start in range(0, n_snps, window_snps)
    ]


def _slice_block_window(
    block: list[list[int]],
    start: int,
    end: int,
) -> list[list[int]]:
    return [row[start:end] for row in block]


def packed_prs_benchmark(
    genotypes: list[list[int]],
    weights: list[float],
    block_size: int = 4,
    poly_modulus_degree: int = DEFAULT_POLY_MODULUS_DEGREE,
    scale: float = DEFAULT_SCALE,
    coeff_mod_bit_sizes: tuple[int, ...] = DEFAULT_COEFF_MOD_BIT_SIZES,
) -> PackedPRSResult:
    """Run a packed multi-sample CKKS PRS benchmark.

    This SCI prototype concatenates multiple samples into one CKKS vector and
    uses a plaintext block-diagonal matrix to recover one PRS score per packed
    sample. It is a baseline packed design, not yet the final optimized method.
    """

    if block_size <= 0:
        raise ValueError("block_size must be positive")
    n_snps = _validate_rectangular(genotypes)
    if len(weights) != n_snps:
        raise ValueError("weights length must match genotype SNP count")

    slots_per_ciphertext = poly_modulus_degree // 2
    packed_slots = block_size * n_snps
    if packed_slots > slots_per_ciphertext:
        raise ValueError(
            "block_size * SNP count exceeds available CKKS slots for one ciphertext"
        )

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

    matrix = _block_weight_matrix(weights, block_size)
    blocks = [
        genotypes[index:index + block_size]
        for index in range(0, len(genotypes), block_size)
    ]

    start = perf_counter()
    encrypted_blocks = [
        ts.ckks_vector(context, _flatten_block(block, block_size, n_snps))
        for block in blocks
    ]
    encryption_seconds = perf_counter() - start

    start = perf_counter()
    encrypted_scores = [encrypted_block.mm(matrix) for encrypted_block in encrypted_blocks]
    evaluation_seconds = perf_counter() - start

    start = perf_counter()
    scores: list[float] = []
    for block, encrypted_score in zip(blocks, encrypted_scores):
        scores.extend(float(value) for value in encrypted_score.decrypt()[:len(block)])
    decryption_seconds = perf_counter() - start

    start = perf_counter()
    public_context_bytes = len(context.serialize(save_secret_key=False))
    input_ciphertext_bytes = sum(len(encrypted_block.serialize()) for encrypted_block in encrypted_blocks)
    result_ciphertext_bytes = sum(len(encrypted_score.serialize()) for encrypted_score in encrypted_scores)
    serialization_seconds = perf_counter() - start

    return PackedPRSResult(
        scores=scores,
        context_seconds=context_seconds,
        encryption_seconds=encryption_seconds,
        evaluation_seconds=evaluation_seconds,
        decryption_seconds=decryption_seconds,
        serialization_seconds=serialization_seconds,
        public_context_bytes=public_context_bytes,
        input_ciphertext_bytes=input_ciphertext_bytes,
        result_ciphertext_bytes=result_ciphertext_bytes,
        block_size=block_size,
        blocks=len(blocks),
        slots_per_ciphertext=slots_per_ciphertext,
        poly_modulus_degree=poly_modulus_degree,
        scale=scale,
        coeff_mod_bit_sizes=coeff_mod_bit_sizes,
    )


def packed_prs_aggregate_benchmark(
    genotypes: list[list[int]],
    weights: list[float],
    block_size: int = 4,
    poly_modulus_degree: int = DEFAULT_POLY_MODULUS_DEGREE,
    scale: float = DEFAULT_SCALE,
    coeff_mod_bit_sizes: tuple[int, ...] = DEFAULT_COEFF_MOD_BIT_SIZES,
) -> PackedAggregateResult:
    """Run a packed CKKS benchmark for block-level PRS sums.

    This computes the sum of PRS scores within each packed sample block with one
    encrypted dot product per block. It is useful for cohort/site-level
    aggregate statistics and as an SCI follow-up direction distinct from
    per-sample PRS recovery.
    """

    if block_size <= 0:
        raise ValueError("block_size must be positive")
    n_snps = _validate_rectangular(genotypes)
    if len(weights) != n_snps:
        raise ValueError("weights length must match genotype SNP count")

    slots_per_ciphertext = poly_modulus_degree // 2
    packed_slots = block_size * n_snps
    if packed_slots > slots_per_ciphertext:
        raise ValueError(
            "block_size * SNP count exceeds available CKKS slots for one ciphertext"
        )

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

    weight_vector = _block_weight_vector(weights, block_size)
    blocks = [
        genotypes[index:index + block_size]
        for index in range(0, len(genotypes), block_size)
    ]

    start = perf_counter()
    encrypted_blocks = [
        ts.ckks_vector(context, _flatten_block(block, block_size, n_snps))
        for block in blocks
    ]
    encryption_seconds = perf_counter() - start

    start = perf_counter()
    encrypted_sums = [encrypted_block.dot(weight_vector) for encrypted_block in encrypted_blocks]
    evaluation_seconds = perf_counter() - start

    start = perf_counter()
    block_sums = [float(encrypted_sum.decrypt()[0]) for encrypted_sum in encrypted_sums]
    decryption_seconds = perf_counter() - start

    start = perf_counter()
    public_context_bytes = len(context.serialize(save_secret_key=False))
    input_ciphertext_bytes = sum(len(encrypted_block.serialize()) for encrypted_block in encrypted_blocks)
    result_ciphertext_bytes = sum(len(encrypted_sum.serialize()) for encrypted_sum in encrypted_sums)
    serialization_seconds = perf_counter() - start

    return PackedAggregateResult(
        block_sums=block_sums,
        context_seconds=context_seconds,
        encryption_seconds=encryption_seconds,
        evaluation_seconds=evaluation_seconds,
        decryption_seconds=decryption_seconds,
        serialization_seconds=serialization_seconds,
        public_context_bytes=public_context_bytes,
        input_ciphertext_bytes=input_ciphertext_bytes,
        result_ciphertext_bytes=result_ciphertext_bytes,
        block_size=block_size,
        blocks=len(blocks),
        slots_per_ciphertext=slots_per_ciphertext,
        poly_modulus_degree=poly_modulus_degree,
        scale=scale,
        coeff_mod_bit_sizes=coeff_mod_bit_sizes,
    )


def packed_windowed_prs_aggregate_benchmark(
    site_genotypes: list[list[list[int]]],
    weights: list[float],
    block_size: int = 4,
    window_snps: int | None = None,
    poly_modulus_degree: int = DEFAULT_POLY_MODULUS_DEGREE,
    scale: float = DEFAULT_SCALE,
    coeff_mod_bit_sizes: tuple[int, ...] = DEFAULT_COEFF_MOD_BIT_SIZES,
) -> WindowedPackedAggregateResult:
    """Run a batched, context-reusing windowed packed aggregate benchmark.

    All site/window ciphertexts are encrypted under one CKKS context and then
    evaluated as one benchmark batch. The output is a site-level PRS sum, so this
    is for federated aggregate statistics rather than individual PRS recovery.
    """

    if block_size <= 0:
        raise ValueError("block_size must be positive")
    if not site_genotypes:
        raise ValueError("site_genotypes must not be empty")
    if not weights:
        raise ValueError("weights must not be empty")

    n_snps = len(weights)
    for site_rows in site_genotypes:
        if not site_rows:
            raise ValueError("site genotype partitions must not be empty")
        site_snps = _validate_rectangular(site_rows)
        if site_snps != n_snps:
            raise ValueError("weights length must match genotype SNP count")

    slots_per_ciphertext = poly_modulus_degree // 2
    if window_snps is None:
        window_snps = max(1, slots_per_ciphertext // block_size)
    if window_snps <= 0:
        raise ValueError("window_snps must be positive")
    if block_size * window_snps > slots_per_ciphertext:
        raise ValueError(
            "block_size * window SNP count exceeds available CKKS slots for one ciphertext"
        )

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

    windows = aggregate_window_ranges(n_snps, window_snps)
    weight_vectors = {
        (start, end): _block_weight_vector(weights[start:end], block_size)
        for start, end in windows
    }

    encrypted_tasks = []
    blocks_total = 0
    start_time = perf_counter()
    for site_index, site_rows in enumerate(site_genotypes):
        blocks = [
            site_rows[index:index + block_size]
            for index in range(0, len(site_rows), block_size)
        ]
        for window_start, window_end in windows:
            n_window_snps = window_end - window_start
            for block in blocks:
                encrypted = ts.ckks_vector(
                    context,
                    _flatten_block(
                        _slice_block_window(block, window_start, window_end),
                        block_size,
                        n_window_snps,
                    ),
                )
                encrypted_tasks.append(
                    (site_index, encrypted, weight_vectors[(window_start, window_end)])
                )
                blocks_total += 1
    encryption_seconds = perf_counter() - start_time

    start_time = perf_counter()
    encrypted_sums = [
        (site_index, encrypted.dot(weight_vector))
        for site_index, encrypted, weight_vector in encrypted_tasks
    ]
    evaluation_seconds = perf_counter() - start_time

    site_sums = [0.0 for _ in site_genotypes]
    start_time = perf_counter()
    for site_index, encrypted_sum in encrypted_sums:
        site_sums[site_index] += float(encrypted_sum.decrypt()[0])
    decryption_seconds = perf_counter() - start_time

    start_time = perf_counter()
    public_context_bytes = len(context.serialize(save_secret_key=False))
    input_ciphertext_bytes = sum(
        len(encrypted.serialize())
        for _site_index, encrypted, _weight_vector in encrypted_tasks
    )
    result_ciphertext_bytes = sum(
        len(encrypted_sum.serialize())
        for _site_index, encrypted_sum in encrypted_sums
    )
    serialization_seconds = perf_counter() - start_time

    return WindowedPackedAggregateResult(
        site_sums=site_sums,
        context_seconds=context_seconds,
        encryption_seconds=encryption_seconds,
        evaluation_seconds=evaluation_seconds,
        decryption_seconds=decryption_seconds,
        serialization_seconds=serialization_seconds,
        public_context_bytes=public_context_bytes,
        input_ciphertext_bytes=input_ciphertext_bytes,
        result_ciphertext_bytes=result_ciphertext_bytes,
        block_size=block_size,
        blocks=blocks_total,
        windows_per_site=len(windows),
        ciphertexts=len(encrypted_tasks),
        calls=1,
        slots_per_ciphertext=slots_per_ciphertext,
        window_snps=window_snps,
        poly_modulus_degree=poly_modulus_degree,
        scale=scale,
        coeff_mod_bit_sizes=coeff_mod_bit_sizes,
    )
