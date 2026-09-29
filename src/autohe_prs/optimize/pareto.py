"""Pareto filtering for measured or predicted candidate records."""

from __future__ import annotations

from typing import Any, Iterable


def pareto_frontier(
    records: Iterable[dict[str, Any]],
    objectives: tuple[str, ...] = (
        "evaluation_seconds",
        "peak_memory_mb",
        "ciphertext_bytes",
        "key_bytes",
        "max_absolute_error",
    ),
) -> list[dict[str, Any]]:
    rows = [row for row in records if all(row.get(key) is not None for key in objectives)]
    frontier: list[dict[str, Any]] = []
    for candidate in rows:
        dominated = False
        for other in rows:
            if other is candidate:
                continue
            no_worse = all(float(other[key]) <= float(candidate[key]) for key in objectives)
            strictly_better = any(float(other[key]) < float(candidate[key]) for key in objectives)
            if no_worse and strictly_better:
                dominated = True
                break
        if not dominated:
            frontier.append(candidate)
    return sorted(frontier, key=lambda row: tuple(float(row[key]) for key in objectives))
