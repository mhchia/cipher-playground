"""
Multivariate sumcheck protocol over hypercube [D]^l.
"""
import itertools
import random as rand
from sage.all import *

from ring import Fq, q


def sum_over_hypercube(poly, xs, fixed: dict, D: int, start: int, end: int):
    """Substitute fixed values, then sum over xs[start:end] \\in [d]^(end-start)."""
    result = 0
    for b in itertools.product(range(D), repeat=end-start):
        subs = {**fixed, **{xs[start+i]: b[i] for i in range(end-start)}}
        result += poly.subs(subs)
    return result


def sumcheck(lde_poly, xs, a_0: int, D: int) -> tuple[Fq, Fq, list[Fq]]:
    # Claim: \sum_{b_0}...\sum_{b_{l-1}} f(b_0, ..., b_{l-1}) = a_j
    # P calculate the first a_j (a_0) and send it to V
    l = len(xs)
    # a_j, a at loop j
    # sum over [d]^l
    a = a_0

    received_randoms = []
    for j in range(l):
        #
        # Prover
        #
        # Let the variate j be X and sum the rest of the variates over [d].
        # let h_j(x) = \sum_{b_{j+1}} ... \sum_{b_{l-1}} f(r_0, ..., r_{j-1}, x, b_{j+1}, ..., b_{l-1})
        # P calculate h_j(x) and send it to V as g_j(x)
        rs = {xs[i]: v for i, v in enumerate(received_randoms)}

        # TODO: should use bookkeeping to avoid re-evaluate polynomials. I.e. next round we should reuse the multivariate polynomial
        # from the last round.
        # round i: f_i(r_0,...,r_{i-1}, x_i, x_{i+1}, ...)
        # round i+1: f_i(x_i=r_i)

        # h_j(x) = f(r_0, ..., r_{j-1}, x, 0, ..., 0) + f(r_0, ..., r_{j-1}, x, 0, ..., 1) +...
        g = sum_over_hypercube(lde_poly, xs, rs, D, start=j+1, end=l)

        # Send `g_j` to V

        #
        # Verifier
        #

        # V is not sure if g_j(x) = h_j(x) as P claimed
        #   and needs to verify
        #   1. a_j = g_j(0) + ... + g_j(d-1)
        #   2. g_j(r) ?= \sum_{b_{j+1}} ... \sum_{b_{l-1}} f(r_0, ..., r_j, b_{j+1}, ..., b_{l-1}), by SZPL
        #       - this is done by P running sumcheck again

        # Step 1. a_j = g_j(0) + ... + g_j(d-1)
        assert a == sum([g.subs({xs[j]: i}) for i in range(D)])

        # Step 2. g_j(r) ?= \sum_{b_{j+1}} ... \sum_{b_{l-1}} f(r_0, ..., r, b_{j+1}, ..., b_{l-1})
        # Sample r_j and calculate a_j = g_j(r)
        r = rand.randint(0, q-1)
        a = g.subs({xs[j]: r})
        # Send r_j to P

        # P: store r_j for next round
        received_randoms.append(r)
    # NOTE: After the loop, V has `a` (a_l = g_{l-1}(r_{l-1}))
    # Since V is not sure if g_{l-1}(x) = f(r_0,...,r_{l-2}, x),
    # V still needs to check `a ?= f(r_0, ..., r_{l-1})`
    return a, received_randoms
