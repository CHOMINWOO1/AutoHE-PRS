"""Run the top-K valid ranked plans repeatedly in an OpenFHE environment."""

from __future__ import annotations

import argparse
from typing import Any
import json

from autohe_prs.benchmark.runner import run_isolated_repeated
from autohe_prs.benchmark.storage import read_records
from autohe_prs.config import FHEPlan, load_experiment_spec
from autohe_prs.genomics.io import read_harmonized


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--ranking", required=True)
    parser.add_argument("--spec", required=True)
    parser.add_argument("--top-k", type=int, default=3)
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--timeout", type=float, default=300.0)
    parser.add_argument("--pgs-id", default="")
    parser.add_argument("--source-domain", default="")
    args = parser.parse_args(argv)

    spec = load_experiment_spec(args.spec)
    _, dosages, weights = read_harmonized(
        spec.workload.genotype_path, spec.workload.weight_path
    )
    ranked = [row for row in read_records(args.ranking) if row.get("feasible")]
    selected = ranked[: args.top_k]
    if len(selected) != args.top_k:
        raise ValueError(f"ranking has only {len(selected)} feasible plans")

    records: list[dict[str, Any]] = []
    for candidate in selected:
        plan = FHEPlan(**candidate["plan"])
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
            value["predicted_rank"] = candidate["rank"]
            value["predicted_targets"] = candidate.get("predictions", {})
            value["prediction_target_scopes"] = candidate.get(
                "prediction_target_scopes", {}
            )
            value["search_stage"] = "top_k_actual_validation"
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
    return 0 if all(row.get("status") == "ok" for row in records) else 2


if __name__ == "__main__":
    raise SystemExit(main())
