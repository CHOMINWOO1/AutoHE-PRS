"""Repeated-run descriptive statistics."""

from __future__ import annotations

from dataclasses import dataclass
import statistics


@dataclass(frozen=True)
class SummaryStats:
    median: float
    mean: float
    standard_deviation: float
    minimum: float
    maximum: float
    count: int


def summarize(values: list[float]) -> SummaryStats:
    if not values:
        raise ValueError("cannot summarize an empty metric")
    return SummaryStats(
        median=statistics.median(values),
        mean=statistics.fmean(values),
        standard_deviation=statistics.stdev(values) if len(values) > 1 else 0.0,
        minimum=min(values),
        maximum=max(values),
        count=len(values),
    )
