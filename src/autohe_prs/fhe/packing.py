"""Slot, chunk, and packing calculations with no backend side effects."""

from __future__ import annotations

import math


def scheme_slot_capacity(scheme: str, ring_dimension: int) -> int:
    scheme = scheme.upper()
    if scheme == "CKKS":
        return ring_dimension // 2
    if scheme in {"BFV", "BGV"}:
        return ring_dimension
    raise ValueError(f"unsupported scheme: {scheme}")


def chunk_ranges(variants: int, variants_per_ciphertext: int) -> tuple[tuple[int, int], ...]:
    if variants < 0 or variants_per_ciphertext <= 0:
        raise ValueError("variants must be non-negative and chunk size positive")
    return tuple(
        (start, min(start + variants_per_ciphertext, variants))
        for start in range(0, variants, variants_per_ciphertext)
    )


def chunk_count(variants: int, variants_per_ciphertext: int) -> int:
    return math.ceil(variants / variants_per_ciphertext) if variants else 0


def padding_mask(active_slots: int, slot_count: int) -> tuple[int, ...]:
    if active_slots < 0 or slot_count < active_slots:
        raise ValueError("active slots must be within slot capacity")
    return tuple([1] * active_slots + [0] * (slot_count - active_slots))


def reduction_rotation_indices(active_slots: int) -> tuple[int, ...]:
    if active_slots <= 0:
        return ()
    rotations: list[int] = []
    offset = 1
    while offset < active_slots:
        rotations.append(offset)
        offset <<= 1
    return tuple(rotations)


def operation_counts(
    *,
    samples: int,
    variants: int,
    variants_per_ciphertext: int,
    aggregate_strategy: str,
) -> dict[str, int]:
    chunks = chunk_count(variants, variants_per_ciphertext)
    ranges = chunk_ranges(variants, variants_per_ciphertext)
    rotations_per_sample = sum(
        len(reduction_rotation_indices(end - start)) for start, end in ranges
    )
    if aggregate_strategy == "add_then_reduce" and chunks:
        rotations_per_sample = len(reduction_rotation_indices(variants_per_ciphertext))
    return {
        "ciphertext_count": samples * chunks,
        "multiplication_count": samples * chunks,
        "rotation_count": samples * rotations_per_sample,
        "addition_count": samples * (rotations_per_sample + max(0, chunks - 1)),
    }
