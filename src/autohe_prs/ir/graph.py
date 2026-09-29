"""Construction and validation of the canonical encrypted PRS graph."""

from __future__ import annotations

from dataclasses import dataclass

from .nodes import IRNode, NodeKind


@dataclass(frozen=True)
class PRSGraph:
    nodes: tuple[IRNode, ...]

    def validate(self) -> None:
        seen: set[str] = set()
        for node in self.nodes:
            if node.node_id in seen:
                raise ValueError(f"duplicate IR node: {node.node_id}")
            missing = set(node.inputs) - seen
            if missing:
                raise ValueError(f"IR inputs must precede {node.node_id}: {sorted(missing)}")
            seen.add(node.node_id)


def reduction_rotations(slots_used: int) -> tuple[int, ...]:
    indices: list[int] = []
    value = 1
    while value < slots_used:
        indices.append(value)
        value *= 2
    return tuple(indices)


def build_prs_graph(
    samples: int,
    variants: int,
    window_size: int,
    *,
    scores: int = 1,
    scheme: str = "CKKS",
) -> PRSGraph:
    if samples < 0 or variants < 0 or scores <= 0 or window_size <= 0:
        raise ValueError("samples/variants must be non-negative and scores/window positive")
    chunks = max(1, (variants + window_size - 1) // window_size)
    rotations = reduction_rotations(min(variants, window_size))
    nodes: tuple[IRNode, ...] = (
        IRNode("genotypes", NodeKind.GENOTYPE_INPUT, shape=(samples, variants),
               numeric_range=(0.0, 2.0)),
        IRNode("weights", NodeKind.WEIGHT_INPUT, shape=(variants, scores),
               value_type="real"),
        IRNode("oriented", NodeKind.EFFECT_ALLELE_ORIENT, ("genotypes",),
               shape=(samples, variants)),
        IRNode("filtered", NodeKind.VARIANT_FILTER, ("oriented",),
               shape=(samples, variants)),
        IRNode("missing_mask", NodeKind.MISSING_MASK, ("filtered",),
               shape=(samples, variants), attributes={"policy": "configured"}),
        IRNode("chunks", NodeKind.CHUNK, ("missing_mask",),
               shape=(samples, chunks, window_size),
               attributes={"chunk_count": chunks, "window_size": window_size}),
        IRNode("encoded_weights", NodeKind.ENCODE_PLAIN, ("weights",),
               shape=(chunks, window_size, scores), slot_count=window_size,
               ckks_scale_bits=50 if scheme.upper() == "CKKS" else None),
        IRNode("encrypted", NodeKind.ENCRYPT, ("chunks",), encrypted=True,
               shape=(samples, chunks, window_size), slot_count=window_size),
        IRNode("weighted", NodeKind.MUL_PLAIN, ("encrypted", "encoded_weights"),
               encrypted=True, shape=(samples, chunks, window_size, scores),
               slot_count=window_size,
               ckks_level=1 if scheme.upper() == "CKKS" else None,
               multiplicative_depth=1, estimated_operation_count=samples * chunks),
        IRNode("reduced", NodeKind.REDUCE_SUM, ("weighted",), encrypted=True,
               shape=(samples, chunks, scores),
               estimated_operation_count=samples * chunks * scores * len(rotations),
               required_rotation_indices=rotations),
        IRNode("aggregated", NodeKind.CHUNK_AGGREGATE, ("reduced",), encrypted=True,
               shape=(samples, scores),
               estimated_operation_count=samples * scores * max(0, chunks - 1)),
        IRNode("output", NodeKind.OUTPUT, ("aggregated",), encrypted=True,
               shape=(samples, scores)),
    )
    graph = PRSGraph(nodes)
    graph.validate()
    return graph
