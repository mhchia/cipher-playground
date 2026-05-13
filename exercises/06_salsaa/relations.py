"""
SALSAA relation dataclasses (scaffolding).

Each sub-protocol is a reduction Σ_input → Σ_output. We separate the public
statement (Instance) from the witness (Witness) so verifier code can be
type-checked to never touch the witness.

NOTE: Not yet wired into rok_norm / rok_bar_sum — schema only.
Reference: SALSAA paper, Section 3-4.
"""

from sage.all import *
from sage.structure.element import Matrix as SageMatrix

from dataclasses import dataclass, replace
from typing import Any

from ring import to_centered, Rq


def _assert_matrix(name, obj):
    """Catch the common bug where a vector or list slips in where a matrix is expected.
    Sage `==` returns False for type mismatches silently, so an upfront check
    makes failures explicit instead of misleading 'relation doesn't hold'."""
    assert isinstance(obj, SageMatrix), \
        f"{name} must be a Sage Matrix, got {type(obj).__name__}: {obj!r}"


# ============================================================
# Σ^lin: linear relation with column-wise L2 norm bound
# ((H, F, Y, β), W) ∈ Σ^lin  iff  H · F · W = Y mod q  and
#   max_i ‖w_i‖₂ ≤ β    (paper notation; β is the column L2 norm bound)
# where w_i is the i-th column of W viewed as concatenated coefficient vector.
# β is stored unsquared; the dataclass invariant compares `max_col_norm² ≤ β²`.
# ============================================================
@dataclass(frozen=True)
class LinInstance:
    H: Any   # n̂ × n in R_q
    # Upper F, commitment
    F_com: Any   # \bar n × m in R_q
    # Lower F, evaluation
    F_eval: Any  # \underbar n x m in R_q
    Y: Any   # n̂ × r in R_q
    beta: int    # column L2 norm bound: max_i ‖w_i‖₂ ≤ beta (unsquared)

    @property
    def F(self):
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
        new_F_eval = self.F_eval.stack(new_F_rows)
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
        _assert_matrix("H", self.H)
        _assert_matrix("F_com", self.F_com)
        _assert_matrix("Y", self.Y)
        if self.F_eval is None:
            # Overwrite F_eval with a zero matrix
            object.__setattr__(
                self, 'F_eval',
                zero_matrix(Rq, 0, self.F_com.ncols()),
            )
        _assert_matrix("F_eval", self.F_eval)
        assert self.F_com.ncols() == self.F_eval.ncols(), \
            f"F_com/F_eval width mismatch: F_com cols={self.F_com.ncols()}, F_eval cols={self.F_eval.ncols()}"



@dataclass(frozen=True)
class LinWitness:
    W: Any   # m × r in R_q

    def __post_init__(self):
        _assert_matrix("W", self.W)

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

        # Dimension consistency BEFORE the math check — catches shape bugs with
        # a precise message rather than a silent False from `matrix == vector`.
        assert H.ncols() == F.nrows(), \
            f"H.ncols ({H.ncols()}) must equal F.nrows ({F.nrows()})"
        assert F.ncols() == W.nrows(), \
            f"F.ncols ({F.ncols()}) must equal W.nrows ({W.nrows()})"
        assert (H.nrows(), W.ncols()) == Y.dimensions(), \
            f"Y dimensions {Y.dimensions()} must equal (H.nrows, W.ncols) = ({H.nrows()}, {W.ncols()})"

        assert H * (F * W) == Y, f"relation doesn't hold: {H=}, {F=}, {W=}, {Y=}"

        max_col_norm_square = 0
        for i in range(W.ncols()):
            col_norm_square = 0
            for j in range(W.nrows()):
                for c in W[j][i].list():
                    cv = to_centered(c)
                    col_norm_square += cv * cv
            max_col_norm_square = max(col_norm_square, max_col_norm_square)
        assert max_col_norm_square <= self.instance.beta ** 2, \
            f"column l_2 norm^2 = {max_col_norm_square} > \\beta² = {self.instance.beta ** 2}"

    @property
    def hat_n(self):
        return self.instance.hat_n

    @property
    def n(self):
        return self.instance.n

    @property
    def n_top(self):
        return self.instance.F_com.nrows()

    @property
    def m(self):
        return self.instance.m

    @property
    def r(self):
        return self.instance.r

    @property
    def beta(self):
        return self.instance.beta
