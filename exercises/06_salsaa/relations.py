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
        return self.Y.ncols()

    def with_extra_eval(self, new_F_rows, new_Y_rows):
        """
        Append eval rows
        """
        if self.F_eval is None:
            new_F_eval = new_F_rows
        else:
            new_F_eval = self.F_eval.stack(new_F_rows)
        # TODO: fix H. Now we just extend H with identity matrix.
        n_new_F_rows = new_F_rows.nrows()
        n_hat = self.H.nrows()
        n = self.F.nrows()
        new_H = block_matrix(Rq, [
            [self.H,                          zero_matrix(Rq, n_hat, n_new_F_rows)],
            [zero_matrix(Rq, n_new_F_rows, n),      identity_matrix(Rq, n_new_F_rows)   ],
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
