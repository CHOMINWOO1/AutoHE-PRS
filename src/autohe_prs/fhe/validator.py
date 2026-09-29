"""Deterministic rejection of unsafe or structurally invalid FHE plans."""

from __future__ import annotations

from dataclasses import dataclass
import math

from autohe_prs.config import Constraints, FHEPlan
from .packing import reduction_rotation_indices, scheme_slot_capacity


# Conservative HomomorphicEncryption.org-style upper bounds for ~128-bit
# classical security. OpenFHE remains authoritative at context generation.
MAX_TOTAL_MODULUS_BITS = {
    1024: 27,
    2048: 54,
    4096: 109,
    8192: 218,
    16384: 438,
    32768: 881,
    65536: 1761,
}


@dataclass(frozen=True)
class ValidationIssue:
    stage: str
    failure_type: str
    reason: str
    severity: str = "error"


@dataclass(frozen=True)
class ValidationResult:
    valid: bool
    issues: tuple[ValidationIssue, ...]
    expected_score_bound: float
    integer_dot_bound: int | None
    estimated_memory_mb: float
    authoritative_backend_check_required: bool = True


def _power_of_two(value: int) -> bool:
    return value > 0 and value & (value - 1) == 0


def _is_prime(value: int) -> bool:
    """Deterministic Miller-Rabin for unsigned 64-bit integers."""

    if value < 2:
        return False
    small = (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37)
    if value in small:
        return True
    if any(value % prime == 0 for prime in small):
        return False
    divisor, shifts = value - 1, 0
    while divisor % 2 == 0:
        divisor //= 2
        shifts += 1
    for base in (2, 325, 9375, 28178, 450775, 9780504, 1795265022):
        if base % value == 0:
            continue
        current = pow(base, divisor, value)
        if current in {1, value - 1}:
            continue
        for _ in range(shifts - 1):
            current = pow(current, 2, value)
            if current == value - 1:
                break
        else:
            return False
    return True


def estimate_memory_mb(plan: FHEPlan, samples: int, variants: int) -> float:
    chunks = max(1, math.ceil(variants / max(1, plan.variants_per_ciphertext)))
    total_modulus_bits = _total_modulus_bits(plan)
    # Two/three polynomial ciphertext components plus evaluation/key material.
    component_bytes = plan.ring_dimension * max(1, total_modulus_bits) / 8
    ciphertext_bytes = 3 * component_bytes
    key_bytes = component_bytes * (8 + 2 * len(plan.rotation_indices))
    return (samples * chunks * ciphertext_bytes + key_bytes) / (1024 * 1024)


def _total_modulus_bits(plan: FHEPlan) -> int:
    if plan.scheme.upper() == "CKKS":
        first = plan.first_modulus_size or 60
        scaling = plan.scaling_modulus_size or plan.scale_bits or 50
        return first + scaling * max(1, plan.multiplicative_depth + 1)
    # OpenFHE selects RNS ciphertext moduli from depth; this is a conservative
    # estimate used only for early rejection.
    return 60 * max(2, plan.multiplicative_depth + 1)


def validate_plan(
    plan: FHEPlan,
    *,
    samples: int,
    variants: int,
    genotype_min: float,
    genotype_max: float,
    sum_abs_weights: float,
    constraints: Constraints | None = None,
    imputed_dosage: bool = False,
) -> ValidationResult:
    constraints = constraints or Constraints()
    plan = plan.normalized()
    issues: list[ValidationIssue] = []
    if plan.security_bits != constraints.security_bits or constraints.security_bits != 128:
        issues.append(ValidationIssue("security", "unsupported_security_profile",
                                      "MVP requires the OpenFHE 128-bit classic profile"))
    if not _power_of_two(plan.ring_dimension) or plan.ring_dimension < 1024:
        issues.append(ValidationIssue("context", "invalid_ring_dimension",
                                      "ring dimension must be a power of two >= 1024"))
    slot_capacity = scheme_slot_capacity(plan.scheme, plan.ring_dimension)
    if plan.slot_count not in {0, slot_capacity}:
        issues.append(ValidationIssue("packing", "slot_count_mismatch",
                                      f"{plan.scheme} capacity is {slot_capacity}"))
    if plan.window_size <= 0 or plan.window_size > slot_capacity:
        issues.append(ValidationIssue("packing", "window_exceeds_slots",
                                      f"window {plan.window_size} exceeds {slot_capacity} slots"))
    if plan.variants_per_ciphertext <= 0 or plan.variants_per_ciphertext > slot_capacity:
        issues.append(ValidationIssue("packing", "variants_per_ciphertext_invalid",
                                      "variants per ciphertext must fit slot capacity"))
    if plan.batch_size not in {0, plan.window_size}:
        issues.append(ValidationIssue(
            "packing", "batch_window_mismatch",
            "MVP batch size must equal the SNP window size",
        ))
    expected_chunks = max(1, math.ceil(variants / max(1, plan.variants_per_ciphertext)))
    if plan.chunk_count != expected_chunks:
        issues.append(ValidationIssue(
            "packing", "chunk_count_mismatch",
            f"plan has {plan.chunk_count}, expected {expected_chunks}",
        ))
    if not plan.feature_packing:
        issues.append(ValidationIssue(
            "packing", "unsupported_non_feature_packing",
            "deterministic MVP runner requires feature packing",
        ))
    if plan.sample_packing or plan.samples_per_ciphertext != 1:
        issues.append(ValidationIssue(
            "packing", "unsupported_sample_packing",
            "sample packing is represented in the IR but not implemented by the MVP runner",
        ))
    if plan.score_packing or plan.scores_per_ciphertext != 1:
        issues.append(ValidationIssue(
            "packing", "unsupported_score_packing",
            "multi-score packing is not implemented by the MVP runner",
        ))
    required_rotations = reduction_rotation_indices(min(variants, plan.variants_per_ciphertext))
    absent = set(required_rotations) - set(plan.rotation_indices)
    if absent:
        issues.append(ValidationIssue("keys", "missing_rotation_keys",
                                      f"missing rotation indices: {sorted(absent)}"))
    invalid_rotations = [
        index for index in plan.rotation_indices
        if index <= 0 or index >= slot_capacity or not _power_of_two(index)
    ]
    if invalid_rotations:
        issues.append(ValidationIssue(
            "keys", "invalid_rotation_index",
            f"rotation indices must be positive powers of two below capacity: {invalid_rotations}",
        ))
    if plan.multiplicative_depth < 1:
        issues.append(ValidationIssue("levels", "insufficient_multiplicative_depth",
                                      "MulPlain PRS requires depth >= 1"))
    if plan.chunk_aggregation_strategy not in {"reduce_then_add", "add_then_reduce"}:
        issues.append(ValidationIssue("graph", "unsupported_chunk_aggregation",
                                      "aggregation must be reduce_then_add or add_then_reduce"))
    if plan.reduction_tree_strategy not in {"binary_tree", "eval_sum"}:
        issues.append(ValidationIssue(
            "graph",
            "unsupported_reduction_strategy",
            "reduction must be binary_tree or eval_sum",
        ))
    if plan.key_switching_technique.upper() not in {"HYBRID", "BV"}:
        issues.append(ValidationIssue(
            "keys",
            "unsupported_key_switching_technique",
            "OpenFHE runner supports HYBRID or BV key switching",
        ))
    total_bits = _total_modulus_bits(plan)
    bit_limit = MAX_TOTAL_MODULUS_BITS.get(plan.ring_dimension)
    if bit_limit is None or total_bits > bit_limit:
        issues.append(ValidationIssue("security", "modulus_budget_exceeded",
                                      f"estimated {total_bits} bits exceeds limit {bit_limit}"))

    expected_bound = max(abs(genotype_min), abs(genotype_max)) * sum_abs_weights
    integer_bound: int | None = None
    if plan.scheme == "CKKS":
        if (plan.scaling_technique or "FLEXIBLEAUTO").upper() not in {
            "FIXEDMANUAL",
            "FIXEDAUTO",
            "FLEXIBLEAUTO",
            "FLEXIBLEAUTOEXT",
        }:
            issues.append(ValidationIssue(
                "encoding",
                "unsupported_scaling_technique",
                "unsupported CKKS scaling technique",
            ))
        scale_bits = plan.scale_bits or plan.scaling_modulus_size
        if scale_bits is None or scale_bits < 20:
            issues.append(ValidationIssue("encoding", "ckks_scale_too_small",
                                          "CKKS requires at least 20 scale bits"))
        if scale_bits and plan.scaling_modulus_size and scale_bits > plan.scaling_modulus_size:
            issues.append(ValidationIssue("encoding", "ckks_scale_mismatch",
                                          "scale bits exceed scaling modulus size"))
        if (
            plan.first_modulus_size is not None
            and scale_bits is not None
            and plan.first_modulus_size < scale_bits
        ):
            issues.append(ValidationIssue(
                "encoding", "ckks_first_modulus_too_small",
                "first modulus must be at least as large as the initial scale",
            ))
    else:
        if imputed_dosage:
            issues.append(ValidationIssue(
                "encoding", "unsupported_imputed_fixed_point",
                "MVP BFV/BGV accepts hard-call integer dosages only",
            ))
        scale = plan.fixed_point_scaling_factor
        modulus = plan.plaintext_modulus
        if not scale or scale <= 0:
            issues.append(ValidationIssue("encoding", "invalid_fixed_point_scale",
                                          "positive fixed-point scale is required"))
        if not modulus or modulus <= 2:
            issues.append(ValidationIssue("encoding", "invalid_plaintext_modulus",
                                          "plaintext modulus must be greater than 2"))
        elif not _is_prime(modulus):
            issues.append(ValidationIssue("encoding", "plaintext_modulus_not_prime",
                                          "packed BFV/BGV plaintext modulus must be prime"))
        elif modulus % (2 * plan.ring_dimension) != 1:
            issues.append(ValidationIssue(
                "encoding", "plaintext_modulus_not_batching_compatible",
                f"p mod 2N must equal 1, got {modulus % (2 * plan.ring_dimension)}",
            ))
        if (
            plan.scheme == "BGV"
            and modulus
            and modulus >= 2**30
            and plan.ring_dimension < 16_384
        ):
            issues.append(ValidationIssue(
                "security",
                "bgv_ring_dimension_below_openfhe_minimum",
                "OpenFHE HEStd_128_classic requires ring dimension >= 16384 "
                "for the configured large BGV plaintext modulus",
            ))
        if scale and modulus:
            integer_bound = math.ceil(expected_bound * scale + variants * 0.5)
            if 2 * integer_bound >= modulus:
                issues.append(ValidationIssue("encoding", "fixed_point_wraparound",
                                              f"2*bound={2 * integer_bound} >= p={modulus}"))

    memory = estimate_memory_mb(plan, samples, variants)
    if constraints.max_memory_mb is not None and memory > constraints.max_memory_mb:
        issues.append(ValidationIssue("resources", "estimated_memory_exceeded",
                                      f"{memory:.1f} MB > {constraints.max_memory_mb:.1f} MB"))
    return ValidationResult(
        valid=not any(issue.severity == "error" for issue in issues),
        issues=tuple(issues),
        expected_score_bound=expected_bound,
        integer_dot_bound=integer_bound,
        estimated_memory_mb=memory,
    )
