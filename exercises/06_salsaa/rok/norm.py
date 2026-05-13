"""
Norm-check sub-protocol for SALSAA.

  rok_norm  : prove ‖w_i‖_{a,2} ≤ ν via t_i = ⟨w_i, w̄_i⟩
  rok_bar_sum : sumcheck on RLC of CRT(LDE[W] · LDE[W̄]),
                output verified via s_0 = LDE[W](r̃), s_1 = LDE[W](r̃̄)
"""
import sys
from pathlib import Path
sys.path.append(str(Path(__file__).resolve().parent.parent.parent / "02_ntt"))

import random as rand
from sage.all import *
from ntt import ntt, intt

from ring import (
    q, d, Fq, Rq, x, e, D,
    conjugate
)
from relations import LinRelation
from lde import lde_poly, tensor_product
from sumcheck import sumcheck


def get_u_vec(Fq, r: int, d: int, e: int):
    u = Fq(rand.randint(1, q-1))
    u_vec = vector(Fq, [u ** i for i in range(r*d//e)])
    return u_vec


def rok_bar_sum(r: int, t: list[Rq], W: matrix) -> tuple[tuple[list[Rq], list[Rq]], tuple[list[Rq], list[Rq]]]:
    #
    # Verifier
    #

    # Challenges in RLC for all NTT slots
    # We have a LDE for each column w_i \in W=[w_1, ..., w_r] and
    # we split that LDE into d/e NTT slots f_\text{slot_0}, ..., f_\text{slot_3}
    # So in total there are r*d/e slots (F_{q^e})
    # u^T: [u^0, u^1, ..., u^{rd/e}]
    u_T = get_u_vec(Fq, r, d, e)

    # [t_{0,0}, ..., t_{0,d/e}, ..., t_{r-1,0}, ..., t_{r-1,d/e}]
    t_ntt = [
        t_i_s
        for t_i in t  # t_i \in R_q
        for t_i_s in ntt(t_i.list(), Rq)
    ]
    # Sanity check
    assert len(u_T) == len(t_ntt)
    # [u^0 * t_{0,0} + ... + u^{rd/e} * t_{r-1,d/e}]
    a_0 = sum([_u * _t for _u, _t in zip(u_T, t_ntt)])

    #
    # Prover
    #

    # We handle one column a time, and concatenate all of them in the end.
    def calculate_CRT_LDE_w_LDE_bar_w(col_idx: int):
        w = W.column(col_idx)
        w_bar = vector(Rq, [conjugate(w_i) for w_i in w])
        # Calculate LDE[w_i]
        lde_w, xs = lde_poly(w, D)
        # Calculate LDE[\bar w_i]
        lde_w_bar, _ = lde_poly(w_bar, D)
        lde_w_times_w_bar = lde_w * lde_w_bar

        # \tilde f (x_0, x_1) = A + B x_0 + C x_1 + D x_0 x_1  , A, B, C, D \in R_q
        # NOTE: poly.coefficients() and monomials() order mismatch!
        # poly.coefficients(): returns in the ascending order: [1, x_0, x_1, x_0 x_1]
        # poly.monomials(): returns in descending order: [x_0 x_1, x_1, x_0, 1], opposite to `coefficients()`
        # Timing [a*b for a, b zip(coefficients(), monomials())] would be wrong!
        # Use poly.dict() instead when traversing multivariates!
        monomials = []
        coeff_rqs = []
        for key, rq_coeff in lde_w_times_w_bar.items():
            # turn 1 into (1,) for simplicity
            try:
                monomial_exp_tuple = (int(key),)
            # already tuple
            except TypeError:
                monomial_exp_tuple = key
            # (1, 2) -> x_0^1 * x_1^2
            monomials.append(
                prod(xs[i]**int(monomial_exp_tuple[i]) for i in range(len(xs)))
            )
            coeff_rqs.append(rq_coeff)

        #   [NTT(A), NTT(B), NTT(C), NTT(D)]
        # = [(a_0, a_1, a_2, a_3), (b_0, b_1, b_2, b_3), (c_0, c_1, c_2, c_3), (d_0, d_1, d_2, d_3)]
        coeff_ntts = [ntt(rq.list(), Rq) for rq in coeff_rqs]

        # \tilde f_\text{slot_0} (x_0, x_1) = a_0 + b_0 x_0 + c_0 x_1 + d_0 x_0 x_1
        # ...
        # \tilde f_\text{slot_3} (x_0, x_1) = a_3 + b_3 x_0 + c_3 x_1 + d_3 x_0 x_1
        # tilde_fs = [\tilde f_{slot_0}, ..., \tilde f_{slot_{d//e-1}}]
        tilde_fs = []
        for i in range(d // e):
            # \tilde_f_i = a_0 + b_0 x_0 + c_0 x_1 + d_0 x_0 x_1
            coeff_ntts_slot_i = [slot[i] for slot in coeff_ntts]
            tilde_f_slot_i = sum([
                c * m
                for c, m in zip(coeff_ntts_slot_i, monomials)
            ])
            tilde_fs.append(tilde_f_slot_i)
        return lde_w, tilde_fs, xs

    # Get CRT(LDE[w]*LDE[\bar w]) for each column, and concatenate all of them into one list
    crt_LDE_W_LDE_bar_W = []
    lde_W = []
    for i in range(r):
        lde_w_i, tilde_fs_i, xs = calculate_CRT_LDE_w_LDE_bar_w(i)
        crt_LDE_W_LDE_bar_W.extend(tilde_fs_i)
        lde_W.append(lde_w_i)
    # xs is reused from last iteration, they're all the same
    assert len(crt_LDE_W_LDE_bar_W) != 0, "empty W and thus xs is not existing!"

    assert len(u_T) == len(crt_LDE_W_LDE_bar_W), "u_T length mismatch crt_LDE_W_LDE_bar_W"

    # RLC (Randomly Linear Combination) all r*d/e LDE per NTT slot
    # \tilde f = u^0 \tilde f_0 + ... + u^{r*d/e-1} \tilde f_{r*d/e-1}
    tilde_f = sum([
        u_i * tilde_f_i
        for u_i, tilde_f_i in zip(u_T, crt_LDE_W_LDE_bar_W)
    ])

    #
    # Prover <> Verifier on sumcheck
    #
    a_l, rands = sumcheck(tilde_f, xs, a_0, D)


    #
    # Prover
    #

    # P needs to prove a_l = \tilde f(r_0, ..., r_{l-1})
    # Since \tilde f is the RLC of (LDE[W] * LDE[\bar W]),
    # instead, we check a_l ?= u^T* CRT(s_0 * \bar s_1)
    # since (LDE[W] * LDE[\bar W])(r) = LDE[W](r) * LDE[\bar W](r) = s_0 * \bar s_1

    # Since all r*d/e slots uses the same challenge [r_0, ..., r_{l-1}] in sumcheck,
    # To evaluate them in Rq, we just get r_j_rq = iNTT([r_j, ..., r_j)]),
    # the coefficient form of r_j corresponding in Rq
    # Then we can use r_0_rq, ..., r_{l-1}_rq to evaluate LDE[W] at them.
    l = len(xs)
    r_T = []
    for j in range(l):
        # the coefficient form of the r
        r_coeffs = intt([rands[j]] * (d//e), Rq)
        # Rq coefficient form
        r_T.append(
            sum([
                c * x ** i
                for i, c in enumerate(r_coeffs)
            ])
        )

    # evaluate s_0 = LDE[W](r) \in R_q^r
    s0_rqs = [
        lde_w_i.subs({xs[j]: r_T[j] for j in range(l)})
        for lde_w_i in lde_W
    ]
    # s_0_ntts = [ntt(rq.list(), Rq) for rq in s_0_rqs]

    # evaluate s_1 = LDE[W](\bar r) \in R_q^r
    r_T_bar = [conjugate(r_T[j]) for j in range(l)]
    s1_rqs = [
        lde_w_i.subs({xs[j]: r_T_bar[j] for j in range(l)})
        for lde_w_i in lde_W
    ]
    # s_1_bar_ntts = [ntt(rq.list(), Rq) for rq in s_1_bar_rqs]
    assert len(s0_rqs) == len(s1_rqs)

    # Prover sends s_0, s_1 to Verifier

    #
    # Verifier
    #

    # pair-wise product of s_0 and \bar s_1
    s0_s1_bar_rqs = [
        a * conjugate(b)
        for a, b in zip(s0_rqs, s1_rqs)
    ]

    # CRT(s_0 * \bar s_1)
    s0_s1_bar_ntts = []
    for rq in s0_s1_bar_rqs:
        s0_s1_bar_ntts.extend(ntt(rq.list(), Rq))

    # rhs = u^T * CRT(s_0 * \bar s_1)
    assert len(u_T) == len(s0_s1_bar_ntts), "u^T length mismatch s0_s1_bar_ntts"
    rhs = sum([
        u * s
        for u, s in zip(u_T, s0_s1_bar_ntts)
    ])

    # a_\mu ?= u^T * CRT(s_0 * \bar s_1)
    assert a_l == rhs, f"a_\\mu mismatch u^T * CRT(s_0 * \\bar s_1): {a_l=}, {rhs=}"

    return (
        (r_T, s0_rqs),
        (r_T_bar, s1_rqs),
    )


def rok_norm(lin_r: LinRelation) -> LinRelation:
    W = lin_r.witness.W
    m = lin_r.m
    r = lin_r.r

    #
    # Prover
    #

    # t = (<w_i, \bar w_i>)_{i \in [r]}
    t = []
    for i in range(r):
        w_i = W.column(i)
        t_i = sum([
            w_i[j] * conjugate(w_i[j])
            for j in range(m)
        ])
        t.append(t_i)

    # Sends t to V

    #
    # Verifier: check <w_i, \bar w_i> <= \mu^2
    #

    # 1. Prover sends t and claims t_i = <w_i, \bar_i>
    # Verifier checks these t_i are <= \mu^2
    mu_square = d * (lin_r.beta ** 2)
    for i in range(r):
        # Trace(t_i) = d * ct(t_i)
        trace_t_i = d * int(t[i].list()[0])
        assert trace_t_i <= mu_square, f"Trace(t_i) is not <= \mu^2: {trace_t_i=}, {mu_square=}"

    # 2. P <-> V: Prover reduces "t_i really equals <w_i, \bar w_i>" to
    #  "LDE[W](r_0) = s_0 and LDE[W](r_1) = s_1" with `rok_bar_sum`
    (r_0, s_0), (r_1, s_1) = rok_bar_sum(r, t, W)

    # 3. Embed s_1, s_2 into the existing relation HFW = Y
    new_F_rows = matrix(Rq, [
        tensor_product(r_0, D),
        tensor_product(r_1, D),
    ])
    new_Y_rows = matrix(Rq, [
        s_0,
        s_1,
    ])
    lin_r_new = LinRelation(
        instance=lin_r.instance.with_extra_eval(
            new_F_rows=new_F_rows,
            new_Y_rows=new_Y_rows,
        ),
        witness=lin_r.witness,
    )
    # relation r: \hat n +2, n+2, m, r_\text{acc}+L, \beta
    return lin_r_new
