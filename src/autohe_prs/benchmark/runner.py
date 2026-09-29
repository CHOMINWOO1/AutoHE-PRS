"""Timeout-safe repeated OpenFHE execution using one process per run."""

from __future__ import annotations

from dataclasses import asdict, dataclass
import json
from pathlib import Path
import subprocess
import sys
import tempfile
from typing import Any

from autohe_prs.config import Constraints, FHEPlan
from .hardware import hardware_metadata
from .metrics import SummaryStats, summarize


@dataclass(frozen=True)
class BenchmarkSummary:
    records: tuple[dict[str, Any], ...]
    statistics: dict[str, SummaryStats]
    warmup_status: str | None
    warmup_record: dict[str, Any] | None


def _failure_record(plan: FHEPlan, repeat: int, failure_type: str, reason: str) -> dict[str, Any]:
    return {
        "scheme": plan.scheme.upper(),
        "backend": "openfhe-python",
        "measurement_kind": "measured",
        "status": "failed",
        "failure_type": failure_type,
        "failure_reason": reason,
        "repeat": repeat,
        "plan": asdict(plan),
        **hardware_metadata(),
    }


def run_isolated_repeated(
    plan: FHEPlan,
    dosages: list[list[float]],
    weights: list[float],
    *,
    constraints: Constraints | None = None,
    repeats: int = 3,
    warmup: bool = True,
    timeout_seconds: float = 300.0,
) -> BenchmarkSummary:
    if repeats <= 0:
        raise ValueError("repeats must be positive")
    request = {
        "plan": asdict(plan),
        "constraints": asdict(constraints or Constraints()),
        "dosages": dosages,
        "weights": weights,
    }
    records: list[dict[str, Any]] = []
    warmup_status: str | None = None
    warmup_record: dict[str, Any] | None = None
    total_runs = repeats + int(warmup)
    with tempfile.TemporaryDirectory(prefix="autohe-prs-") as raw:
        directory = Path(raw)
        request_path = directory / "request.json"
        request_path.write_text(json.dumps(request), encoding="utf-8")
        for index in range(total_runs):
            is_warmup = warmup and index == 0
            repeat = index if warmup else index + 1
            result_path = directory / f"result-{index}.json"
            try:
                process = subprocess.run(
                    [sys.executable, "-m", "autohe_prs.benchmark.worker",
                     str(request_path), str(result_path)],
                    capture_output=True,
                    text=True,
                    timeout=timeout_seconds,
                    check=False,
                )
                if result_path.exists():
                    row = json.loads(result_path.read_text(encoding="utf-8"))
                    row["repeat"] = 0 if is_warmup else repeat
                    row["warmup"] = is_warmup
                    row.update(hardware_metadata())
                else:
                    row = _failure_record(
                        plan, repeat, "process_crash",
                        f"exit={process.returncode}; stderr={process.stderr[-2000:]}",
                    )
            except subprocess.TimeoutExpired as exc:
                row = _failure_record(plan, repeat, "timeout", str(exc))
            row["warmup"] = is_warmup
            if is_warmup:
                row["repeat"] = 0
            if is_warmup:
                warmup_status = str(row.get("status"))
                warmup_record = row
            else:
                records.append(row)
    metric_names = (
        "context_seconds", "encryption_seconds", "evaluation_seconds",
        "key_generation_seconds", "multiplication_key_generation_seconds",
        "rotation_key_generation_seconds",
        "encoding_seconds", "decryption_seconds", "decoding_seconds",
        "serialization_seconds", "total_seconds", "peak_memory_mb",
        "ciphertext_bytes", "communication_bytes", "key_bytes", "crypto_context_bytes",
        "public_key_bytes", "secret_key_bytes", "multiplication_key_bytes",
        "rotation_key_bytes", "max_absolute_error",
        "max_relative_error",
    )
    statistics: dict[str, SummaryStats] = {}
    for metric in metric_names:
        values = [
            float(row[metric]) for row in records
            if row.get("status") == "ok" and row.get(metric) is not None
        ]
        if values:
            statistics[metric] = summarize(values)
    return BenchmarkSummary(tuple(records), statistics, warmup_status, warmup_record)
