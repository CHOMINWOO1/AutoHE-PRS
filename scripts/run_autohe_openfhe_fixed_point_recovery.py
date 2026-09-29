"""Measure conservative BFV/BGV recovery plans after failed small-ring plans."""

from __future__ import annotations

import argparse
from dataclasses import asdict
import json
import math
from typing import Any

from autohe_prs.benchmark.runner import run_isolated_repeated
from autohe_prs.config import FHEPlan, load_experiment_spec
from autohe_prs.fhe.packing import reduction_rotation_indices, scheme_slot_capacity
from autohe_prs.fhe.validator import validate_plan
from autohe_prs.genomics.io import read_harmonized


def _plan(
    scheme: str,
    variants: int,
    *,
    ring_dimension: int,
    window_size: int,
    fixed_point_scale: int,
    plaintext_modulus: int,
) -> FHEPlan:
    return FHEPlan(
        scheme=scheme,
        ring_dimension=ring_dimension,
        slot_count=scheme_slot_capacity(scheme, ring_dimension),
        batch_size=window_size,
        window_size=window_size,
        variants_per_ciphertext=window_size,
        chunk_count=math.ceil(variants / window_size),
        rotation_indices=reduction_rotation_indices(window_size),
        chunk_aggregation_strategy="reduce_then_add",
        reduction_tree_strategy="binary_tree",
        plaintext_modulus=plaintext_modulus,
        fixed_point_scaling_factor=fixed_point_scale,
        encoding_precision=max(0, round(math.log10(fixed_point_scale))),
        key_switching_technique="HYBRID",
        multiplicative_depth=1,
    )


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--spec", required=True)
    parser.add_argument("--schemes", nargs="+", default=("BFV", "BGV"))
    parser.add_argument("--ring-dimension", type=int, default=16_384)
    parser.add_argument("--window-size", type=int, default=4_096)
    parser.add_argument("--fixed-point-scale", type=int, default=100_000_000)
    parser.add_argument("--plaintext-modulus", type=int, default=2_147_352_577)
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--timeout", type=float, default=300.0)
    parser.add_argument("--pgs-id", default="")
    parser.add_argument("--source-domain", default="")
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args(argv)

    spec = load_experiment_spec(args.spec)
    _, dosages, weights = read_harmonized(
        spec.workload.genotype_path, spec.workload.weight_path
    )
    flat = [value for row in dosages for value in row]
    records: list[dict[str, Any]] = []
    for scheme in args.schemes:
        plan = _plan(
            scheme,
            len(weights),
            ring_dimension=args.ring_dimension,
            window_size=args.window_size,
            fixed_point_scale=args.fixed_point_scale,
            plaintext_modulus=args.plaintext_modulus,
        )
        validation = validate_plan(
            plan,
            samples=len(dosages),
            variants=len(weights),
            genotype_min=min(flat),
            genotype_max=max(flat),
            sum_abs_weights=math.fsum(abs(value) for value in weights),
            constraints=spec.constraints,
            imputed_dosage=False,
        )
        if not validation.valid:
            records.append(
                {
                    "scheme": scheme,
                    "plan": asdict(plan),
                    "measurement_kind": "validation_only",
                    "status": "rejected",
                    "failure_type": "constraint_validation",
                    "failure_reason": "; ".join(
                        issue.reason for issue in validation.issues
                    ),
                    "validation": asdict(validation),
                }
            )
            continue
        if args.dry_run:
            records.append(
                {
                    "scheme": scheme,
                    "plan": asdict(plan),
                    "measurement_kind": "validation_only",
                    "status": "ok",
                    "validation": asdict(validation),
                }
            )
            continue
        summary = run_isolated_repeated(
            plan,
            dosages,
            weights,
            constraints=spec.constraints,
            repeats=args.repeats,
            warmup=False,
            timeout_seconds=args.timeout,
        )
        for record in summary.records:
            value = dict(record)
            value["search_stage"] = "fixed_point_recovery"
            if args.pgs_id:
                value["pgs_id"] = args.pgs_id
            if args.source_domain:
                value["source_domain"] = args.source_domain
            if args.pgs_id or args.source_domain:
                value["workload_group"] = ":".join(
                    item
                    for item in (args.source_domain, args.pgs_id)
                    if item
                )
            records.append(value)

    print(json.dumps(records, indent=2, sort_keys=True))
    return 0 if records and all(row.get("status") == "ok" for row in records) else 2


if __name__ == "__main__":
    raise SystemExit(main())
