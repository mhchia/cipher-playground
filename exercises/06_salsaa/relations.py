"""
SALSAA relation dataclasses (scaffolding).

Each sub-protocol is a reduction Σ_input → Σ_output. We separate the public
statement (Instance) from the witness (Witness) so verifier code can be
type-checked to never touch the witness.

NOTE: Not yet wired into rok_norm / rok_bar_sum — schema only.
Reference: SALSAA paper, Section 3-4.
"""
from dataclasses import dataclass
from typing import Any


# ============================================================
# Σ^lin: linear relation with infinity-norm bound
# ((H, F, Y, β), W) ∈ Σ^lin  iff  H · F · W = Y mod q  and  ‖W‖_∞ ≤ β
# ============================================================
@dataclass(frozen=True)
class LinInstance:
    H: Any   # n̂ × n in R_q
    F: Any   # n × m in R_q
    Y: Any   # n̂ × r in R_q
    beta: int


@dataclass(frozen=True)
class LinWitness:
    W: Any   # m × r in R_q


@dataclass(frozen=True)
class LinRelation:
    instance: LinInstance
    witness: LinWitness


# ============================================================
# Σ^norm: a-2 norm relation
# ((H, F, Y, ν), W) ∈ Σ^norm  iff  ((H,F,Y),W) ∈ Σ^lin  and
#   ‖w_i‖_{a,2} ≤ ν  ∀ i ∈ [r]
# Equivalently: Trace(⟨w_i, w̄_i⟩) ≤ ν²
# ============================================================
@dataclass(frozen=True)
class NormInstance:
    H: Any
    F: Any
    Y: Any
    nu: int   # a-2 norm bound


@dataclass(frozen=True)
class NormWitness:
    W: Any


@dataclass(frozen=True)
class NormRelation:
    instance: NormInstance
    witness: NormWitness


# ============================================================
# Σ^bar_sum (input to rok_bar_sum):
# Norm relation augmented by t = (⟨w_i, w̄_i⟩)_{i ∈ [r]} which has been
# revealed by Prover. Verifier checks Trace(t_i) ≤ ν² before sumcheck.
# ============================================================
@dataclass(frozen=True)
class BarSumInstance:
    norm: NormInstance
    t: tuple   # (t_0, ..., t_{r-1}) ∈ R_q^r


@dataclass(frozen=True)
class BarSumRelation:
    instance: BarSumInstance
    witness: NormWitness


# ============================================================
# Σ^lde⊗ (output of rok_bar_sum):
# After sumcheck, the claim is reduced to evaluations at challenge r̃.
# Verifier holds (s_0, s_1) ∈ R_q^r × R_q^r where:
#   s_0 = LDE[W](r̃)       and        s_1 = LDE[W](r̃̄)
# ============================================================
@dataclass(frozen=True)
class LDETensorInstance:
    norm: NormInstance
    s_0: tuple   # LDE[W](r̃)
    s_1: tuple   # LDE[W](r̃̄)
    r_T: tuple   # sumcheck challenges lifted into R_q


@dataclass(frozen=True)
class LDETensorRelation:
    instance: LDETensorInstance
    witness: NormWitness
