"""Typed configuration shared by validation, benchmarking, and emission."""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
import json
from pathlib import Path
from typing import Any


@dataclass(frozen=True)
class WorkloadSpec:
    operation: str = "prs_weighted_sum"
    genotype_path: str = ""
    weight_path: str = ""
    genotype_encoding: str = "hard_call"
    variants: int | None = None
    samples: int | None = None
    scores: int = 1
    missing_policy: str = "error"
    genome_build: str | None = None


@dataclass(frozen=True)
class Constraints:
    security_bits: int = 128
    max_absolute_error: float | None = 1.0e-6
    max_relative_error: float | None = 1.0e-5
    max_memory_mb: float | None = 32768
    max_latency_seconds: float | None = None
    require_decryption_success: bool = True


@dataclass(frozen=True)
class FHEPlan:
    scheme: str
    ring_dimension: int
    security_bits: int = 128
    multiplicative_depth: int = 1
    batch_size: int = 0
    slot_count: int = 0
    window_size: int = 4096
    chunk_count: int = 1
    variants_per_ciphertext: int = 4096
    samples_per_ciphertext: int = 1
    scores_per_ciphertext: int = 1
    feature_packing: bool = True
    sample_packing: bool = False
    score_packing: bool = False
    chunk_aggregation_strategy: str = "reduce_then_add"
    reduction_tree_strategy: str = "binary_tree"
    rotation_indices: tuple[int, ...] = ()
    key_switching_technique: str = "HYBRID"
    scaling_modulus_size: int | None = None
    first_modulus_size: int | None = None
    scale_bits: int | None = None
    scaling_technique: str | None = None
    rescale_strategy: str | None = None
    plaintext_modulus: int | None = None
    fixed_point_scaling_factor: int | None = None
    encoding_precision: int | None = None
    value_bound_policy: str = "conservative_l1"
    metadata: dict[str, Any] = field(default_factory=dict)

    def normalized(self) -> "FHEPlan":
        scheme = self.scheme.upper()
        if scheme not in {"CKKS", "BFV", "BGV"}:
            raise ValueError(f"unsupported FHE scheme: {self.scheme}")
        return FHEPlan(**{**asdict(self), "scheme": scheme})


@dataclass(frozen=True)
class ExperimentSpec:
    workload: WorkloadSpec
    constraints: Constraints = field(default_factory=Constraints)
    objective: str = "evaluation_latency"
    allowed_schemes: tuple[str, ...] = ("CKKS", "BFV", "BGV")
    top_k_validation: int = 10
    use_bayesian_refinement: bool = False
    hardware_profile: str = "local_cpu"


def _scalar(value: str) -> Any:
    value = value.strip()
    if value in {"null", "Null", "NULL", "~"}:
        return None
    if value.lower() in {"true", "false"}:
        return value.lower() == "true"
    if value.startswith(("'", '"')) and value.endswith(value[0]):
        return value[1:-1]
    try:
        return int(value)
    except ValueError:
        try:
            return float(value)
        except ValueError:
            return value


def _minimal_yaml(text: str) -> dict[str, Any]:
    """Parse the small mapping/list subset used by AutoHE-PRS examples."""

    tokens: list[tuple[int, str]] = []
    for raw in text.splitlines():
        content = raw.split("#", 1)[0].rstrip()
        if not content.strip():
            continue
        indent = len(content) - len(content.lstrip())
        tokens.append((indent, content.strip()))

    def parse_block(index: int, indent: int) -> tuple[Any, int]:
        if index >= len(tokens):
            return {}, index
        is_list = tokens[index][1].startswith("- ")
        container: Any = [] if is_list else {}
        while index < len(tokens):
            current_indent, content = tokens[index]
            if current_indent < indent:
                break
            if current_indent > indent:
                raise ValueError(f"unexpected YAML indentation near: {content}")
            if is_list:
                if not content.startswith("- "):
                    break
                container.append(_scalar(content[2:]))
                index += 1
                continue
            if content.startswith("- ") or ":" not in content:
                break
            key, raw_value = content.split(":", 1)
            key, raw_value = key.strip(), raw_value.strip()
            index += 1
            if raw_value:
                container[key] = _scalar(raw_value)
            elif index < len(tokens) and tokens[index][0] > current_indent:
                container[key], index = parse_block(index, tokens[index][0])
            else:
                container[key] = {}
        return container, index

    parsed, consumed = parse_block(0, tokens[0][0] if tokens else 0)
    if consumed != len(tokens) or not isinstance(parsed, dict):
        raise ValueError("configuration root must be a YAML mapping")
    return parsed


def load_mapping(path: str | Path) -> dict[str, Any]:
    path = Path(path)
    text = path.read_text(encoding="utf-8-sig")
    if path.suffix.lower() == ".json":
        value = json.loads(text)
    else:
        try:
            import yaml  # type: ignore[import-not-found]
        except ImportError:
            value = _minimal_yaml(text)
        else:
            value = yaml.safe_load(text)
    if not isinstance(value, dict):
        raise ValueError("configuration root must be a mapping")
    return value


def load_experiment_spec(path: str | Path) -> ExperimentSpec:
    raw = load_mapping(path)
    workload = WorkloadSpec(**raw.get("workload", {}))
    constraints = Constraints(**raw.get("constraints", {}))
    optimization = raw.get("optimization", {})
    hardware = raw.get("hardware", {})
    spec = ExperimentSpec(
        workload=workload,
        constraints=constraints,
        objective=optimization.get("objective", "evaluation_latency"),
        allowed_schemes=tuple(optimization.get("allowed_schemes", ("CKKS", "BFV", "BGV"))),
        top_k_validation=int(optimization.get("top_k_validation", 10)),
        use_bayesian_refinement=bool(optimization.get("use_bayesian_refinement", False)),
        hardware_profile=hardware.get("profile", "local_cpu"),
    )
    if spec.workload.operation != "prs_weighted_sum":
        raise ValueError("MVP supports only workload.operation=prs_weighted_sum")
    if spec.workload.genotype_encoding not in {"hard_call", "imputed"}:
        raise ValueError("genotype_encoding must be hard_call or imputed")
    if spec.workload.missing_policy not in {
        "error", "exclude_variant", "zero", "mean_dosage"
    }:
        raise ValueError("unsupported workload missing_policy")
    if spec.workload.scores <= 0:
        raise ValueError("workload scores must be positive")
    if spec.constraints.security_bits != 128:
        raise ValueError("MVP requires constraints.security_bits=128")
    if spec.top_k_validation <= 0:
        raise ValueError("top_k_validation must be positive")
    valid_objectives = {
        "evaluation_latency",
        "total_latency",
        "memory",
        "ciphertext_size",
        "key_size",
        "balanced_multi_objective",
    }
    if spec.objective not in valid_objectives:
        raise ValueError(f"unsupported optimization objective: {spec.objective}")
    if not spec.allowed_schemes or any(
        scheme.upper() not in {"CKKS", "BFV", "BGV"}
        for scheme in spec.allowed_schemes
    ):
        raise ValueError("allowed_schemes must be a non-empty CKKS/BFV/BGV list")
    return spec


def dump_json(value: object, path: str | Path) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    if hasattr(value, "__dataclass_fields__"):
        value = asdict(value)  # type: ignore[arg-type]
    path.write_text(json.dumps(value, indent=2, sort_keys=True), encoding="utf-8")
