"""
SALSAA relation dataclasses (scaffolding).

Each sub-protocol is a reduction Σ_input → Σ_output. We separate the public
statement (Instance) from the witness (Witness) so verifier code can be
type-checked to never touch the witness.

NOTE: Not yet wired into rok_norm / rok_bar_sum — schema only.
Reference: SALSAA paper, Section 3-4.
"""

from sage.all import *

from dataclasses import dataclass, replace
from typing import Any

from ring import to_centered, Rq


# ============================================================
# Σ^lin: linear relation with column-wise L2 norm bound
# ((H, F, Y, ν), W) ∈ Σ^lin  iff  H · F · W = Y mod q  and
#   max_i ‖w_i‖₂² ≤ ν²    (paper notation: ‖W‖_{σ,2} ≤ ν)
# where w_i is the i-th column of W viewed as concatenated coefficient vector.
# ============================================================
@dataclass(frozen=True)
class LinInstance:
    H: Any   # n̂ × n in R_q
    # Upper F, commitment
    F_com: Any   # \bar n × m in R_q
    # Lower F, evaluation
    F_eval: Any  # \underbar n x m in R_q
    Y: Any   # n̂ × r in R_q
    v_square: int

    @property
    def F(self):
        if self.F_eval is None:
            return self.F_com
        return self.F_com.stack(self.F_eval)

    @property
    def hat_n(self):
        return self.H.nrows()

    @property
    def n(self):
        return self.H.ncols()

    @property
    def m(self):
        return self.F_com.ncols()

    @property
    def r(self):
        return self.Y.nrows()

    def with_extra_eval(self, new_F_rows, new_Y_rows):
        """
        Append eval rows
        """
        if self.F_eval is None:
            new_F_eval = new_F_rows
        else:
            new_F_eval = self.F_eval.stack(new_F_rows)
        # TODO: fix H. Now we just extend H with identity matrix.
        n_hat = self.H.nrows()
        n = self.F.nrows()
        new_H = block_matrix(Rq, [
            [self.H,                          zero_matrix(Rq, n_hat, 2)],
            [zero_matrix(Rq, 2, n),      identity_matrix(Rq, 2)   ],
        ])
        new_Y = self.Y.stack(new_Y_rows)
        return replace(
            self,
            F_eval=new_F_eval,
            H=new_H,
            Y=new_Y,
        )

    def __post_init__(self):
        if self.F_eval is not None:
            print(f"{self.F_eval=}")
            assert self.F_com.ncols() == self.F_eval.ncols()



@dataclass(frozen=True)
class LinWitness:
    W: Any   # m × r in R_q

    def with_extra_cols(self, new_w_cols):
        """
        Append eval rows
        """
        new_W = self.W.augment(new_w_cols)
        return replace(self, W=new_W)

    @property
    def m(self):
        return self.W.nrows()

    @property
    def r(self):
        return self.W.ncols()


@dataclass(frozen=True)
class LinRelation:
    instance: LinInstance
    witness: LinWitness


    def __post_init__(self):
        H, F, W, Y = self.instance.H, self.instance.F, self.witness.W, self.instance.Y
        assert H * (F * W) == Y

        max_col_norm_square = 0
        for i in range(W.ncols()):
            col_norm_square = 0
            for j in range(W.nrows()):
                for c in W[j][i].list():
                    cv = to_centered(c)
                    col_norm_square += cv * cv
            max_col_norm_square = max(col_norm_square, max_col_norm_square)
        assert max_col_norm_square <= self.instance.v_square

    @property
    def hat_n(self):
        return self.instance.hat_n

    @property
    def n(self):
        return self.instance.n

    @property
    def m(self):
        return self.instance.m

    @property
    def r(self):
        return self.instance.r

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
    v_square: int

    def __post_init__(self):
        H, F, W, Y = self.instance.H, self.instance.F, self.witness.W, self.instance.Y
        assert H * (F * W) == Y

        max_col_norm_square = 0
        for i in range(W.ncols()):
            col_norm_square = 0
            for j in range(W.nrows()):
                for c in W[j][i].list():
                    cv = to_centered(c)
                    col_norm_square += cv * cv
            max_col_norm_square = max(col_norm_square, max_col_norm_square)
        assert max_col_norm_square <= self.instance.v_square



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
    s_0: tuple    # LDE[W](r̃)        ∈ R_q^r  (r = W.ncols())
    s_1: tuple    # LDE[W](r̃̄)        ∈ R_q^r
    r_tilde: tuple   # sumcheck challenges (length μ), lifted into R_q


@dataclass(frozen=True)
class LDETensorRelation:
    instance: LDETensorInstance
    witness: NormWitness
