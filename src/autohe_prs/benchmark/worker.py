"""One-plan worker process. Input/output are paths to avoid pipe truncation."""

from __future__ import annotations

from dataclasses import asdict
import json
from pathlib import Path
import sys

from autohe_prs.config import Constraints, FHEPlan
from autohe_prs.fhe.runner import run_openfhe_plan


def main(argv: list[str] | None = None) -> int:
    argv = argv or sys.argv[1:]
    if len(argv) != 2:
        raise SystemExit("usage: python -m autohe_prs.benchmark.worker REQUEST RESULT")
    request_path, result_path = map(Path, argv)
    request = json.loads(request_path.read_text(encoding="utf-8"))
    plan = FHEPlan(**request["plan"])
    constraints = Constraints(**request.get("constraints", {}))
    result = run_openfhe_plan(plan, request["dosages"], request["weights"], constraints=constraints)
    result_path.write_text(json.dumps(asdict(result), sort_keys=True), encoding="utf-8")
    return 0 if result.status == "ok" else 2


if __name__ == "__main__":
    raise SystemExit(main())
