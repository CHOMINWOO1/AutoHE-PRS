"""Small, serializable PRS intermediate representation."""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from typing import Any


class NodeKind(str, Enum):
    GENOTYPE_INPUT = "GenotypeInput"
    WEIGHT_INPUT = "WeightInput"
    EFFECT_ALLELE_ORIENT = "EffectAlleleOrient"
    VARIANT_FILTER = "VariantFilter"
    MISSING_MASK = "MissingMask"
    ENCODE_PLAIN = "EncodePlain"
    ENCRYPT = "Encrypt"
    MUL_PLAIN = "MulPlain"
    MUL_CIPHER = "MulCipher"
    ADD = "Add"
    ROTATE = "Rotate"
    REDUCE_SUM = "ReduceSum"
    CHUNK = "Chunk"
    CHUNK_AGGREGATE = "ChunkAggregate"
    MULTI_SCORE = "MultiScore"
    DECRYPT = "Decrypt"
    OUTPUT = "Output"


@dataclass(frozen=True)
class IRNode:
    node_id: str
    kind: NodeKind
    inputs: tuple[str, ...] = ()
    value_type: str = "real"
    encrypted: bool = False
    shape: tuple[int, ...] = ()
    slot_count: int = 0
    numeric_range: tuple[float, float] | None = None
    ckks_level: int | None = None
    ckks_scale_bits: int | None = None
    multiplicative_depth: int = 0
    estimated_operation_count: int = 0
    required_rotation_indices: tuple[int, ...] = ()
    attributes: dict[str, Any] = field(default_factory=dict)
