"""OpenFHE execution adapter with explicit provenance and failure records."""

from __future__ import annotations

from dataclasses import asdict, dataclass
from datetime import datetime, timezone
import math
import platform
import sys
from time import perf_counter
import traceback
from typing import Any

from autohe_prs.config import Constraints, FHEPlan
from autohe_prs.benchmark.hardware import hardware_metadata
from autohe_prs.genomics.plaintext import compare, prs
from .validator import ValidationResult, validate_plan


@dataclass(frozen=True)
class ExecutionResult:
    timestamp_utc: str
    scheme: str
    backend: str
    measurement_kind: str
    status: str
    failure_type: str | None
    failure_reason: str | None
    samples: int
    variants: int
    scores: tuple[float, ...]
    reference_scores: tuple[float, ...]
    context_seconds: float | None
    key_generation_seconds: float | None
    multiplication_key_generation_seconds: float | None
    rotation_key_generation_seconds: float | None
    encoding_seconds: float | None
    encryption_seconds: float | None
    evaluation_seconds: float | None
    decryption_seconds: float | None
    decoding_seconds: float | None
    serialization_seconds: float | None
    total_seconds: float
    peak_memory_mb: float | None
    ciphertext_bytes: int | None
    communication_bytes: int | None
    key_bytes: int | None
    key_bytes_scope: str | None
    crypto_context_bytes: int | None
    public_key_bytes: int | None
    secret_key_bytes: int | None
    multiplication_key_bytes: int | None
    rotation_key_bytes: int | None
    max_absolute_error: float | None
    max_relative_error: float | None
    error_constraint_feasible: bool | None
    decryption_success: bool
    plan: dict[str, Any]
    validation: dict[str, Any]
    python_version: str
    platform: str
    cpu_model: str
    physical_core_count: int | None
    logical_core_count: int | None
    ram_bytes: int | None
    operating_system: str
    docker_or_wsl: str
    openfhe_version: str | None


def _peak_memory_mb() -> float | None:
    try:
        if sys.platform == "win32":
            from genome_he_pilot.memory import process_memory
            value = process_memory().peak_rss_bytes
            return value / (1024 * 1024) if value else None
        import resource
        value = float(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
        # Linux reports KiB, macOS bytes.
        return value / (1024 if sys.platform != "darwin" else 1024 * 1024)
    except Exception:
        return None


def _hard_calls(dosages: list[list[float]]) -> list[list[int]]:
    values: list[list[int]] = []
    for row in dosages:
        integer_row: list[int] = []
        for value in row:
            rounded = int(round(value))
            if not math.isclose(value, rounded) or rounded not in {0, 1, 2}:
                raise ValueError("BFV/BGV MVP requires integer hard-call dosage in {0,1,2}")
            integer_row.append(rounded)
        values.append(integer_row)
    return values


def run_openfhe_plan(
    plan: FHEPlan,
    dosages: list[list[float]],
    weights: list[float],
    *,
    constraints: Constraints | None = None,
) -> ExecutionResult:
    started = perf_counter()
    timestamp = datetime.now(timezone.utc).isoformat()
    reference = prs(dosages, weights)
    flat = [value for row in dosages for value in row]
    imputed = any(not math.isclose(value, round(value)) for value in flat)
    validation = validate_plan(
        plan,
        samples=len(dosages),
        variants=len(weights),
        genotype_min=min(flat, default=0.0),
        genotype_max=max(flat, default=0.0),
        sum_abs_weights=math.fsum(abs(weight) for weight in weights),
        constraints=constraints,
        imputed_dosage=imputed,
    )
    hardware = hardware_metadata()
    base = {
        "timestamp_utc": timestamp,
        "scheme": plan.scheme.upper(),
        "backend": "openfhe-python",
        "measurement_kind": "measured",
        "samples": len(dosages),
        "variants": len(weights),
        "reference_scores": tuple(reference),
        "plan": asdict(plan),
        "validation": asdict(validation),
        "python_version": str(hardware["python_version"]),
        "platform": platform.platform(),
        "cpu_model": str(hardware["cpu_model"]),
        "physical_core_count": hardware["physical_core_count"],
        "logical_core_count": hardware["logical_core_count"],
        "ram_bytes": hardware["ram_bytes"],
        "operating_system": str(hardware["operating_system"]),
        "docker_or_wsl": str(hardware["docker_or_wsl"]),
        "openfhe_version": hardware["openfhe_version"],
    }
    if not validation.valid:
        return ExecutionResult(
            **{**base, "measurement_kind": "validation_only"},
            status="rejected", failure_type="constraint_validation",
            failure_reason="; ".join(issue.reason for issue in validation.issues),
            scores=(), context_seconds=None, encryption_seconds=None,
            key_generation_seconds=None, multiplication_key_generation_seconds=None,
            rotation_key_generation_seconds=None,
            encoding_seconds=None,
            evaluation_seconds=None, decryption_seconds=None, decoding_seconds=None,
            serialization_seconds=None,
            total_seconds=perf_counter() - started, peak_memory_mb=None,
            ciphertext_bytes=None, communication_bytes=None, key_bytes=None,
            key_bytes_scope=None, crypto_context_bytes=None,
            public_key_bytes=None, secret_key_bytes=None,
            multiplication_key_bytes=None, rotation_key_bytes=None,
            max_absolute_error=None,
            max_relative_error=None, error_constraint_feasible=None,
            decryption_success=False,
        )
    try:
        from genome_he_pilot.openfhe_controlled_backend import (
            openfhe_ckks_prs_benchmark,
            openfhe_fixed_point_prs_benchmark,
        )
        if plan.scheme.upper() == "CKKS":
            result = openfhe_ckks_prs_benchmark(
                _hard_calls(dosages) if not imputed else dosages,  # type: ignore[arg-type]
                weights,
                window_snps=plan.window_size,
                ring_dim=plan.ring_dimension,
                scaling_mod_size=plan.scaling_modulus_size or 50,
                first_modulus_size=plan.first_modulus_size or 60,
                aggregation_strategy=plan.chunk_aggregation_strategy,
                reduction_strategy=plan.reduction_tree_strategy,
                rotation_indices=plan.rotation_indices,
                key_switching_technique=plan.key_switching_technique,
                scaling_technique=plan.scaling_technique or "FLEXIBLEAUTO",
            )
        else:
            result = openfhe_fixed_point_prs_benchmark(
                _hard_calls(dosages),
                weights,
                scheme=plan.scheme,
                window_snps=plan.window_size,
                ring_dim=plan.ring_dimension,
                plain_modulus=plan.plaintext_modulus or 2_147_352_577,
                fixed_point_scale=plan.fixed_point_scaling_factor or 10_000_000,
                aggregation_strategy=plan.chunk_aggregation_strategy,
                reduction_strategy=plan.reduction_tree_strategy,
                rotation_indices=plan.rotation_indices,
                key_switching_technique=plan.key_switching_technique,
            )
        errors = compare(reference, result.scores)
        constraints = constraints or Constraints()
        total_seconds = perf_counter() - started
        peak_memory_mb = _peak_memory_mb()
        constraint_failures: list[str] = []
        error_constraint_feasible = True
        if (
            constraints.max_absolute_error is not None
            and errors.max_absolute_error > constraints.max_absolute_error
        ):
            constraint_failures.append("maximum absolute error exceeded")
            error_constraint_feasible = False
        if (
            constraints.max_relative_error is not None
            and errors.max_relative_error > constraints.max_relative_error
        ):
            constraint_failures.append("maximum relative error exceeded")
            error_constraint_feasible = False
        if (
            constraints.max_latency_seconds is not None
            and total_seconds > constraints.max_latency_seconds
        ):
            constraint_failures.append("maximum end-to-end latency exceeded")
        if (
            constraints.max_memory_mb is not None
            and peak_memory_mb is not None
            and peak_memory_mb > constraints.max_memory_mb
        ):
            constraint_failures.append("maximum peak memory exceeded")
        constraints_ok = not constraint_failures
        return ExecutionResult(
            **base,
            status="ok" if constraints_ok else "failed",
            failure_type=None if constraints_ok else "measured_hard_constraint",
            failure_reason=None if constraints_ok else "; ".join(constraint_failures),
            scores=tuple(result.scores),
            context_seconds=result.context_seconds,
            key_generation_seconds=result.key_generation_seconds,
            multiplication_key_generation_seconds=result.multiplication_key_generation_seconds,
            rotation_key_generation_seconds=result.rotation_key_generation_seconds,
            encoding_seconds=result.encoding_seconds,
            encryption_seconds=result.encryption_seconds,
            evaluation_seconds=result.evaluation_seconds,
            decryption_seconds=result.decryption_seconds,
            decoding_seconds=result.decoding_seconds,
            serialization_seconds=result.serialization_seconds,
            total_seconds=total_seconds,
            peak_memory_mb=peak_memory_mb,
            ciphertext_bytes=result.input_ciphertext_bytes + result.result_ciphertext_bytes,
            communication_bytes=result.input_ciphertext_bytes + result.result_ciphertext_bytes,
            key_bytes=result.total_key_material_bytes,
            key_bytes_scope="public+secret+multiplication+rotation",
            crypto_context_bytes=result.crypto_context_bytes,
            public_key_bytes=result.public_key_bytes,
            secret_key_bytes=result.secret_key_bytes,
            multiplication_key_bytes=result.multiplication_key_bytes,
            rotation_key_bytes=result.rotation_key_bytes,
            max_absolute_error=errors.max_absolute_error,
            max_relative_error=errors.max_relative_error,
            error_constraint_feasible=error_constraint_feasible,
            decryption_success=True,
        )
    except Exception as exc:  # Every failed configuration is a data record.
        return ExecutionResult(
            **base, status="unavailable" if exc.__class__.__name__.endswith("Unavailable") else "failed",
            failure_type=exc.__class__.__name__,
            failure_reason=f"{exc}\n{traceback.format_exc(limit=3)}",
            scores=(), context_seconds=None, encryption_seconds=None,
            key_generation_seconds=None, multiplication_key_generation_seconds=None,
            rotation_key_generation_seconds=None,
            encoding_seconds=None,
            evaluation_seconds=None, decryption_seconds=None, decoding_seconds=None,
            serialization_seconds=None,
            total_seconds=perf_counter() - started, peak_memory_mb=_peak_memory_mb(),
            # Failure records retain the process peak observed before exit.
            ciphertext_bytes=None, communication_bytes=None, key_bytes=None,
            key_bytes_scope=None, crypto_context_bytes=None,
            public_key_bytes=None, secret_key_bytes=None,
            multiplication_key_bytes=None, rotation_key_bytes=None,
            max_absolute_error=None,
            max_relative_error=None, error_constraint_feasible=None,
            decryption_success=False,
        )
