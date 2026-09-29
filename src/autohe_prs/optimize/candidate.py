"""Reproducible Cartesian candidate generation."""

from __future__ import annotations

from dataclasses import dataclass
from itertools import product
import math

from autohe_prs.config import FHEPlan
from autohe_prs.fhe.packing import reduction_rotation_indices, scheme_slot_capacity


@dataclass(frozen=True)
class CandidateSpace:
    schemes: tuple[str, ...] = ("CKKS", "BFV", "BGV")
    ring_dimensions: tuple[int, ...] = (8192, 16384, 32768)
    window_sizes: tuple[int, ...] = (512, 1024, 2048, 4096)
    scaling_modulus_sizes: tuple[int, ...] = (40, 50)
    first_modulus_sizes: tuple[int, ...] = (60,)
    scaling_techniques: tuple[str, ...] = ("FLEXIBLEAUTO",)
    fixed_point_scales: tuple[int, ...] = (100_000, 1_000_000, 10_000_000)
    plaintext_moduli: tuple[int, ...] = (2_147_352_577,)
    aggregation_strategies: tuple[str, ...] = ("reduce_then_add", "add_then_reduce")
    reduction_strategies: tuple[str, ...] = ("binary_tree",)
    key_switching_techniques: tuple[str, ...] = ("HYBRID",)


def generate_candidates(variants: int, space: CandidateSpace | None = None) -> list[FHEPlan]:
    space = space or CandidateSpace()
    candidates: list[FHEPlan] = []
    for scheme, ring, window, strategy, reduction, key_switching in product(
        space.schemes,
        space.ring_dimensions,
        space.window_sizes,
        space.aggregation_strategies,
        space.reduction_strategies,
        space.key_switching_techniques,
    ):
        scheme = scheme.upper()
        capacity = scheme_slot_capacity(scheme, ring)
        if window > capacity:
            continue
        rotations = reduction_rotation_indices(min(variants, window))
        common = dict(
            scheme=scheme,
            ring_dimension=ring,
            slot_count=capacity,
            batch_size=window,
            window_size=window,
            variants_per_ciphertext=window,
            chunk_count=max(1, math.ceil(variants / window)),
            rotation_indices=rotations,
            chunk_aggregation_strategy=strategy,
            reduction_tree_strategy=reduction,
            key_switching_technique=key_switching,
        )
        if scheme == "CKKS":
            for scale_bits, first_modulus, scaling_technique in product(
                space.scaling_modulus_sizes,
                space.first_modulus_sizes,
                space.scaling_techniques,
            ):
                candidates.append(FHEPlan(
                    **common, scaling_modulus_size=scale_bits,
                    first_modulus_size=first_modulus, scale_bits=scale_bits,
                    scaling_technique=scaling_technique,
                    rescale_strategy="automatic",
                ))
        else:
            for fixed_scale, modulus in product(space.fixed_point_scales, space.plaintext_moduli):
                candidates.append(FHEPlan(
                    **common, plaintext_modulus=modulus,
                    fixed_point_scaling_factor=fixed_scale,
                    encoding_precision=max(0, round(math.log10(fixed_scale))),
                ))
    return candidates
