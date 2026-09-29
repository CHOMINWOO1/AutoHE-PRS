"""Measure manual, random, and uninformed-Bayesian baselines in OpenFHE."""

from __future__ import annotations

import argparse
from dataclasses import asdict
import json
import statistics
from typing import Any

from autohe_prs.benchmark.runner import run_isolated_repeated
from autohe_prs.cli import _space, _workload_features
from autohe_prs.config import FHEPlan, load_experiment_spec
from autohe_prs.fhe.validator import validate_plan
from autohe_prs.genomics.io import read_harmonized
from autohe_prs.optimize.baselines import run_baseline_search
from autohe_prs.optimize.candidate import generate_candidates


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--spec", required=True)
    parser.add_argument("--search-space", required=True)
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--timeout", type=float, default=300.0)
    parser.add_argument("--random-budget", type=int, default=1)
    parser.add_argument("--bayesian-budget", type=int, default=2)
    parser.add_argument("--seed", type=int, default=20260728)
    parser.add_argument("--pgs-id", default="")
    parser.add_argument("--source-domain", default="")
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args(argv)

    spec = load_experiment_spec(args.spec)
    _, dosages, weights = read_harmonized(
        spec.workload.genotype_path, spec.workload.weight_path
    )
    features = _workload_features(dosages, weights)
    candidates = generate_candidates(
        len(weights),
        _space(args.search_space, spec.allowed_schemes),
    )
    feasible: list[FHEPlan] = []
    for plan in candidates:
        result = validate_plan(
            plan,
            samples=len(dosages),
            variants=len(weights),
            genotype_min=features["genotype_min"],
            genotype_max=features["genotype_max"],
            sum_abs_weights=features["sum_absolute_weights"],
            constraints=spec.constraints,
            imputed_dosage=features["imputed_dosage"],
        )
        if result.valid:
            feasible.append(plan)
    if not feasible:
        raise ValueError("the baseline search space has no valid plan")

    manual = next(
        (
            plan
            for plan in feasible
            if plan.scheme.upper() == "CKKS"
            and plan.ring_dimension == 16_384
            and plan.window_size == 4_096
            and plan.scaling_modulus_size == 50
            and plan.first_modulus_size == 60
            and plan.chunk_aggregation_strategy == "reduce_then_add"
            and plan.reduction_tree_strategy == "binary_tree"
            and plan.scaling_technique == "FLEXIBLEAUTO"
            and plan.key_switching_technique == "HYBRID"
        ),
        None,
    )
    if manual is None:
        raise ValueError("the configured space does not contain the manual CKKS baseline")
    if args.dry_run:
        print(
            json.dumps(
                {
                    "candidate_count": len(candidates),
                    "feasible_candidate_count": len(feasible),
                    "manual_plan": asdict(manual),
                },
                indent=2,
                sort_keys=True,
            )
        )
        return 0

    records: list[dict[str, Any]] = []
    search_results: list[dict[str, Any]] = []
    measurement_cache: dict[str, float | None] = {}

    def measure(method: str, plan: FHEPlan) -> float | None:
        cache_key = repr(plan)
        if cache_key in measurement_cache:
            return measurement_cache[cache_key]
        summary = run_isolated_repeated(
            plan,
            dosages,
            weights,
            constraints=spec.constraints,
            repeats=args.repeats,
            warmup=False,
            timeout_seconds=args.timeout,
        )
        successful: list[float] = []
        for record in summary.records:
            row = dict(record)
            row["search_stage"] = f"baseline:{method}"
            row["baseline_method"] = method
            if args.pgs_id:
                row["pgs_id"] = args.pgs_id
            if args.source_domain:
                row["source_domain"] = args.source_domain
            if args.pgs_id or args.source_domain:
                row["workload_group"] = ":".join(
                    value
                    for value in (args.source_domain, args.pgs_id)
                    if value
                )
            records.append(row)
            if row.get("status") == "ok" and row.get("evaluation_seconds") is not None:
                successful.append(float(row["evaluation_seconds"]))
        if len(successful) != args.repeats:
            measurement_cache[cache_key] = None
            return None
        measured = statistics.median(successful)
        measurement_cache[cache_key] = measured
        return measured

    for method, budget, fixed in (
        ("manual_ckks", 1, manual),
        ("random", args.random_budget, None),
        ("bayesian", args.bayesian_budget, None),
    ):
        search_method = "fixed" if fixed is not None else method
        result = run_baseline_search(
            feasible,
            lambda plan, label=method: measure(label, plan),
            method=search_method,
            budget=budget,
            seed=args.seed,
            fixed_plan=fixed,
        )
        payload = asdict(result)
        payload["method"] = method
        search_results.append(payload)

    print(
        json.dumps(
            {
                "candidate_count": len(candidates),
                "feasible_candidate_count": len(feasible),
                "objective": "evaluation_seconds",
                "unique_measured_configurations": len(measurement_cache),
                "records": records,
                "search_results": search_results,
            },
            indent=2,
            sort_keys=True,
        )
    )
    return 0 if all(item.get("best_value") is not None for item in search_results) else 2


if __name__ == "__main__":
    raise SystemExit(main())
