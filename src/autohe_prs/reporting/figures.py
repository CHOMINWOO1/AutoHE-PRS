"""Optional paper figure generation without fabricating missing metrics."""

from __future__ import annotations

from pathlib import Path
from typing import Any


def plot_latency_memory(records: list[dict[str, Any]], path: str | Path) -> Path:
    rows = [
        row for row in records
        if row.get("evaluation_seconds") is not None and row.get("peak_memory_mb") is not None
        and str(row.get("measurement_kind", "")).startswith("measured")
    ]
    if not rows:
        raise ValueError("no rows contain both measured latency and memory")
    try:
        import matplotlib.pyplot as plt  # type: ignore[import-not-found]
    except ImportError as exc:
        raise RuntimeError("figure generation requires matplotlib") from exc
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    figure, axis = plt.subplots(figsize=(6.4, 4.2))
    for scheme in sorted({str(row.get("scheme")) for row in rows}):
        selected = [row for row in rows if str(row.get("scheme")) == scheme]
        axis.scatter(
            [float(row["evaluation_seconds"]) for row in selected],
            [float(row["peak_memory_mb"]) for row in selected],
            label=scheme,
        )
    latencies = [float(row["evaluation_seconds"]) for row in rows]
    positive_latencies = [value for value in latencies if value > 0]
    logarithmic = (
        positive_latencies
        and max(positive_latencies) / min(positive_latencies) >= 100
    )
    if logarithmic:
        axis.set_xscale("log")
    axis.set_xlabel(
        "Measured evaluation latency (s, log scale)"
        if logarithmic
        else "Measured evaluation latency (s)"
    )
    axis.set_ylabel("Measured peak memory (MB)")
    axis.grid(True, which="both", linewidth=0.4, alpha=0.35)
    axis.legend()
    figure.tight_layout()
    figure.savefig(path, dpi=300)
    plt.close(figure)
    return path
