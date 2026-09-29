"""Reference slot reduction used to test generated reduction plans."""

from __future__ import annotations

from .packing import reduction_rotation_indices


def rotate_left(values: list[float], amount: int) -> list[float]:
    if not values:
        return []
    amount %= len(values)
    return values[amount:] + values[:amount]


def binary_tree_reduce(values: list[float], active_slots: int | None = None) -> list[float]:
    """Simulate rotate-and-add; slot zero equals the active-slot sum."""

    if active_slots is None:
        active_slots = len(values)
    if active_slots < 0 or active_slots > len(values):
        raise ValueError("active slot count is outside the vector")
    width = 1
    while width < active_slots:
        width <<= 1
    padded = list(values[:active_slots]) + [0.0] * (width - active_slots)
    for offset in reduction_rotation_indices(active_slots):
        shifted = rotate_left(padded, offset)
        padded = [left + right for left, right in zip(padded, shifted)]
    return padded
