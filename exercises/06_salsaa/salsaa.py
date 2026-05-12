"""
SALSAA — Lattice-based Folding Scheme. Top-level driver.

Module layout:
  ring.py        — R_q setup, parameters, conjugate
  relations.py   — Σ^lin, Σ^norm, Σ^bar_sum, Σ^lde⊗ dataclasses (scaffolding)
  lde.py         — LDE polynomial over [D]^l hypercube
  sumcheck.py    — multivariate sumcheck protocol
  norm_check.py  — rok_norm + rok_bar_sum sub-protocols
  salsaa.py      — main driver (this file)
"""
from sage.all import *
from ring import n_hat, n, d, m, r, _gen_random_low_norm_poly, Rq, beta
from norm_check import rok_norm
from join import rok_join
from rp import rok_rp
from fold import rok_fold
from relations import (
    LinInstance, LinWitness, LinRelation,
)


def gen_random_W(_m: int, _r: int):
    # # Hardcoded for debugging
    # W = matrix(Rq, [
    #     [16*x**2 + 16*x + 1,       x**2 + 16*x + 1],
    #     [16*x**3 + 16*x + 1,           16*x**2 + 1],
    # ])
    # assert W.nrows() == m
    # assert W.ncols() == r
    # return W
    MS = MatrixSpace(Rq, _m, _r)
    W = MS([ _gen_random_low_norm_poly(Rq, beta) for _ in range(_m * _r) ])
    return W


def gen_random_F(_n: int, _m: int):
    MS = MatrixSpace(Rq, _n, _m)
    F = MS([Rq.random_element() for _ in range(_n * _m) ])
    return F


def gen_H(_n_hat: int, _n: int):
    MS = MatrixSpace(Rq, _n_hat, _n)
    H = MS.identity_matrix()
    return H


def fold(lins: list[LinRelation]) -> LinRelation:
    L = len(lins)

    #
    # Join L instances to 1
    #
    lin_joined = lins[0]
    for i in range(1, L):
        lin_joined = rok_join(lin_joined, lins[i])


    #
    # Norm check
    #
    lin_normed = rok_norm(lin_joined)

    assert lin_normed.hat_n == lin_joined.hat_n + 2
    assert lin_normed.n == lin_joined.n + 2
    assert lin_normed.m == lin_joined.m
    assert lin_normed.r == lin_joined.r
    assert lin_normed.v_square == lin_joined.v_square

    #
    # ⊗RP: Perform Johnson-Lindenstrauss to improve soundness without using
    #   a subtractive set.
    #
    # FIXME: Hardcode n_rp and m_rp for now
    # n_rp should be a security parameter!
    n_rp = 1
    m_rp = lin_normed.r * n_rp
    assert m_rp == n_rp * lin_normed.r
    lin_orig, lin_w_hat = rok_rp(lin_normed, n_rp, m_rp)
    assert lin_orig.hat_n == lin_joined.hat_n + 3
    assert lin_orig.n == lin_joined.n + 3
    assert lin_orig.m == lin_joined.m
    assert lin_orig.r == lin_joined.r
    assert lin_orig.v_square == lin_joined.v_square

    # Check lin_w_hat
    assert lin_w_hat.hat_n == lin_joined.n_top + 1
    assert lin_w_hat.n == lin_joined.n_top + 1
    assert lin_w_hat.m == lin_joined.m
    assert lin_w_hat.r == 1
    assert lin_w_hat.v_square > lin_joined.v_square

    #
    # Fold
    #
    # Fold the witnesses of the main statements and output 1 relation
    r_out = 1
    lin_folded = rok_fold(lin_orig, r_out=r_out)

    assert lin_folded.hat_n == lin_joined.hat_n + 3
    assert lin_folded.n == lin_joined.n + 3
    assert lin_folded.m == lin_joined.m
    assert lin_folded.r == 1
    assert lin_folded.v_square > lin_orig.v_square



    return lin_normed


def main():
    # R_q should be small, like {-1,0,1}^d.
    # assume beta^2 = d
    beta_square = d
    # for each w_i \in R^m, |w_i| = |w_{i,1}|^2 + ... + |w_{i,m}|^2 <= beta
    v_square = m * beta_square

    # TODO: now we only have commitment but evaluation. should add evaluation
    # and thus H, F need to be changed.
    # W \in R_q^{m*r}
    H = gen_H(n_hat, n)
    F = gen_random_F(n, m)
    W = gen_random_W(m, r)
    Y = H * F * W
    print(f"{H=}")
    print(f"{F=}")
    print(f"{W=}")
    print(f"{Y=}")

    lin_r = LinRelation(
        instance=LinInstance(H=H, F_com=F, F_eval=None, Y=Y, v_square=v_square),
        witness=LinWitness(W),
    )
    print(f"{lin_r=}")

    #
    # Norm check
    #



    # ============================================================
    # 1. Challenge spaces
    # ============================================================

    # TODO: define C_small, sample challenges




    # ============================================================
    # 5. Decomposition
    # ============================================================

    # TODO: b-ary decomposition


    # ============================================================
    # 6. Norm check
    # ============================================================

    # TODO: post-folding norm verification


    # ============================================================
    # 7. Folding
    # ============================================================


    # ============================================================
    # 4. Support R1CS and CCS (Customizable Constraint System)
    # ============================================================
    # TODO: support R1CS
    # TODO: convert R1CS to CCS


    # TODO: fold two CCS instances

    # # Verify imports work
    # A = ajtai_setup(kappa, m, Rq)
    # z = _gen_random_low_norm_witness(Rq, m, beta)
    # c = ajtai_commit(A, z)
    # # print(f"Ajtai commit OK: {c=}")
