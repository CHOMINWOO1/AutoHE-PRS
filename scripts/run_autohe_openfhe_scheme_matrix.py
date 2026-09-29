"""Run repeated BFV/BGV/CKKS plans from a ranking inside an OpenFHE runtime."""

from __future__ import annotations

import argparse
import json
from typing import Any

from autohe_prs.benchmark.runner import run_isolated_repeated
from autohe_prs.benchmark.storage import read_records
from autohe_prs.config import FHEPlan, load_experiment_spec
from autohe_prs.genomics.io import read_harmonized


def _first_feasible(
    ranking: list[dict[str, Any]], scheme: str
) -> tuple[int, FHEPlan]:
    row = next(
        (
            value
            for value in ranking
            if value.get("feasible")
            and isinstance(value.get("plan"), dict)
            and str(value["plan"].get("scheme", "")).upper() == scheme
        ),
        None,
    )
    if row is None:
        raise ValueError(f"ranking contains no feasible {scheme} plan")
    return int(row["rank"]), FHEPlan(**row["plan"])


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--ranking", required=True)
    parser.add_argument("--spec", required=True)
    parser.add_argument(
        "--schemes",
        nargs="+",
        default=("CKKS", "BFV", "BGV"),
        choices=("CKKS", "BFV", "BGV"),
    )
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--timeout", type=float, default=300.0)
    parser.add_argument("--warmup", action="store_true")
    parser.add_argument("--pgs-id", default="")
    parser.add_argument("--source-domain", default="")
    args = parser.parse_args(argv)

    spec = load_experiment_spec(args.spec)
    _, dosages, weights = read_harmonized(
        spec.workload.genotype_path, spec.workload.weight_path
    )
    ranking = read_records(args.ranking)
    records: list[dict[str, Any]] = []
    for scheme in args.schemes:
        predicted_rank, plan = _first_feasible(ranking, scheme)
        summary = run_isolated_repeated(
            plan,
            dosages,
            weights,
            constraints=spec.constraints,
            repeats=args.repeats,
            warmup=args.warmup,
            timeout_seconds=args.timeout,
        )
        if summary.warmup_record is not None:
            warmup = dict(summary.warmup_record)
            warmup["predicted_rank"] = predicted_rank
            warmup["search_stage"] = "scheme_matrix_warmup"
            records.append(warmup)
        for record in summary.records:
            value = dict(record)
            value["predicted_rank"] = predicted_rank
            value["search_stage"] = "same_input_scheme_matrix"
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
