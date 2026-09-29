"""Auditable Markdown/JSON experiment reporting."""

from __future__ import annotations

from collections import Counter
from datetime import datetime, timezone
import json
from pathlib import Path
from typing import Any


def write_experiment_summary(
    records: list[dict[str, Any]], output_dir: str | Path
) -> tuple[Path, Path]:
    output = Path(output_dir)
    output.mkdir(parents=True, exist_ok=True)
    status = Counter(str(row.get("status", "unknown")) for row in records)
    provenance = Counter(str(row.get("measurement_kind", "unknown")) for row in records)
    failures = Counter(
        str(row.get("failure_type"))
        for row in records
        if row.get("failure_type")
    )
    payload = {
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "record_count": len(records),
        "status_counts": dict(status),
        "measurement_kind_counts": dict(provenance),
        "failure_type_counts": dict(failures),
        "successful_records": [row for row in records if row.get("status") == "ok"],
    }
    json_path = output / "experiment_summary.json"
    json_path.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    lines = [
        "# AutoHE-PRS experiment summary",
        "",
        f"- Records: {len(records)}",
        f"- Status: {dict(status)}",
        f"- Provenance: {dict(provenance)}",
        f"- Failures: {dict(failures)}",
        "",
        "Predicted, simulated, and measured values are never merged.",
    ]
    markdown_path = output / "experiment_summary.md"
    markdown_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return json_path, markdown_path
