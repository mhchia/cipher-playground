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
from ring import d, m, r
from norm_check import get_W, rok_norm
from lde import test_lde_poly


def main():
    # R_q should be small, like {-1,0,1}^d.
    # assume beta^2 = d
    beta_square = d
    # for each w_i \in R^m, |w_i| = |w_{i,1}|^2 + ... + |w_{i,m}|^2 <= beta
    v_square = m * beta_square

    #
    # Norm check
    #

    # FIXME: Mock H, F, Y for now.
    # Should be fixed later when rok_sum_bar is working
    H, F, Y = None, None, None
    # W \in R_q^{m*r}
    W = get_W(m, r)
    # print(f"{W=}")

    rok_norm(H, F, Y, v_square, W)


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
