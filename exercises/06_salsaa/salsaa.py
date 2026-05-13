"""
SALSAA — Lattice-based Folding Scheme. Top-level driver.

Module layout:
  ring.py        — R_q setup, parameters, conjugate
  relations.py   — Σ^lin, Σ^norm, Σ^bar_sum, Σ^lde⊗ dataclasses (scaffolding)
  lde.py         — LDE polynomial over [D]^l hypercube
  sumcheck.py    — multivariate sumcheck protocol
  rok/           — RoK sub-protocols (join, norm, rp, fold, batch, decompose)
  salsaa.py      — main driver (this file)
"""
from sage.all import *
from ring import n_hat, n, d, m, r, _gen_random_low_norm_poly, Rq, beta
from rok import rok_norm, rok_join, rok_rp, rok_fold, rok_batch, rok_decompose, get_l
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


def fold(lins: list[LinRelation], b=2, r_out=1, n_rp=1) -> LinRelation:
    """
    The whole chain of folding schmeme in SALSAA paper.
    Π^join → Π^norm → Π^⊗RP → Π^fold → Π^join → Π^batch → Π^b-decomp
                             →
    Input:  L lin relations
    Output: 1 lin relations

    Params
    - lins: L instances to be folded
    - b (default=2): `b` in b-ary decomposition
    - r_out (default=1): `r_out` relation = 1
    - n_rp (default=1): the n_{rp} in rok rp. `n_rp` should be a security parameter
        - TODO: choose a proper number
    """
    #
    # Configs
    #

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
    # Fold the witnesses of the main statements and output `r_out` relation
    lin_folded = rok_fold(lin_orig, r_out=r_out)

    assert lin_folded.hat_n == lin_joined.hat_n + 3
    assert lin_folded.n == lin_joined.n + 3
    assert lin_folded.m == lin_joined.m
    assert lin_folded.r == 1
    assert lin_folded.v_square > lin_orig.v_square

    #
    # Merge (Join) the relation from ⊗RP (w_hat) and the one from Fold
    #
    lin_merged = rok_join(lin_w_hat, lin_folded)

    assert lin_merged.hat_n == lin_joined.hat_n + 4
    assert lin_merged.n == lin_joined.n + 4
    assert lin_merged.m == lin_joined.m
    assert lin_merged.r == 2
    # Unchanged
    assert lin_merged.v_square == lin_folded.v_square

    #
    # Batch evaluation statements (rows in HF) into less statements
    #
    # Pick c^T to be \hat n - \bar n (of `lin_joined`) rows so we get
    # H \in R^{\hat n \times (n + 4)}
    n_target_eval_rows = lin_joined.hat_n - lin_joined.n_top
    lin_batched = rok_batch(lin_merged, n_target_eval_rows=n_target_eval_rows)

    assert lin_batched.hat_n == lin_joined.hat_n
    assert lin_batched.n == lin_joined.n + 4
    assert lin_batched.m == lin_joined.m
    assert lin_batched.r == 2
    # Unchanged
    assert lin_batched.v_square == lin_merged.v_square

    #
    # b-ary decomposition to lower the norm bound back
    #
    lin_decomposed = rok_decompose(lin_batched, b)

    assert lin_decomposed.hat_n == lin_joined.hat_n
    assert lin_decomposed.n == lin_joined.n + 4
    assert lin_decomposed.m == lin_joined.m

    l = get_l(isqrt(lin_batched.v_square), b)
    # New Y is \tilde Z = [Z_0 || ... || Z_{l-1}]
    # and Z_i \in R^{m x r_old}, so r_new is `r_old * l`
    assert lin_decomposed.r == lin_batched.r * l


    # each coeff \in [-b/2, b/2], one Rq has d coeffs, and one column has m Rqs
    assert lin_decomposed.v_square == lin_decomposed.m * d * ((b // 2) ** 2)

    assert lin_decomposed.v_square <= lin_joined.v_square, \
        "norm budget must not exceed round-start (else next round won't compose)"

    print(f"salsaa.fold finished! Norm^2 goes down: {lin_merged.v_square} → {lin_decomposed.v_square}. Original norm^2: {lin_joined.v_square}")

    return lin_decomposed
