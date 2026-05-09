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
from relations import (
    LinInstance, LinWitness, LinRelation,
    NormInstance, NormWitness, NormRelation,
    BarSumInstance, BarSumRelation, LDETensorInstance, LDETensorRelation
)
from lde import test_lde_poly, tensor_product


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


def rok_join(lin_1: LinRelation, lin_2: LinRelation) -> LinRelation:
    pass


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

    lin_norm_r = rok_norm(lin_r, v_square)


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


def tests():
    test_lde_poly()


if __name__ == '__main__':
    tests()
    main()
