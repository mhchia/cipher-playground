
import random as rand
import itertools

from sage.all import *

# ============================================================
# MLE + Sumcheck over F_q
# ============================================================

q = 17
Fq = GF(q)


def sum_over_hypercube(poly, xs, fixed: dict, d: int, start: int, end: int):
    """Substitute fixed values, then sum over xs[start:end] ∈ [d]^(end-start)."""
    result = 0
    for b in itertools.product(range(d), repeat=end-start):
        subs = {**fixed, **{xs[start+i]: b[i] for i in range(end-start)}}
        result += poly.subs(subs)
    return result


def sumcheck(f: vector, d: int = 2):
    mle_t, xs = lde(f)
    # Claim: \sum_{b_0}...\sum_{b_{l-1}} f(b_0, ..., b_{l-1}) = a_j
    # P calculate the first a_j (a_0) and send it to V
    l = len(xs)
    # a_j, a at loop j
    # sum over [d]^l
    a = sum_over_hypercube(mle_t, xs, {}, d, start=0, end=l)
    res = a

    received_randoms = []
    for j in range(l):
        #
        # Prover
        #
        # Let the variate j be X and sum the rest of the variates over [d].
        # let h_j(x) = \sum_{b_{j+1}} ... \sum_{b_{l-1}} f(r_0, ..., r_{j-1}, x, b_{j+1}, ..., b_{l-1})
        # P calculate h_j(x) and send it to V as g_j(x)
        rs = {xs[i]: v for i, v in enumerate(received_randoms)}
        # h_j(x) = f(r_0, ..., r_{j-1}, x, 0, ..., 0) + f(r_0, ..., r_{j-1}, x, 0, ..., 1) +...
        g = sum_over_hypercube(mle_t, xs, rs, d, start=j+1, end=l)

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
        assert a == sum([g.subs({xs[j]: i}) for i in range(d)])

        # Step 2. g_j(r) ?= \sum_{b_{j+1}} ... \sum_{b_{l-1}} f(r_0, ..., r, b_{j+1}, ..., b_{l-1})
        # Sample r_j and calculate a_j = g_j(r)
        r = rand.randint(0, q-1)
        a = g.subs({xs[j]: r})
        # Send r_j to P

        # P: store r_j for next round
        received_randoms.append(r)

    # V: oracle query to f(r_0, ..., r_{l-1}) and check if it matches a_{j+1}
    f_r = mle_t.subs({
        xs[i]: received_randoms[i]
        for i in range(l)
    })
    assert f_r == a, f"result mismatches oracle query: {f_r=}, {a=}, {received_randoms=}"
    return res


# Pad w to size 2^l to make it on hypercube
def pad_vec_to_d_exp(w: vector, d: int = 2):
    # First padding w to 8 elements so we have [d]^3
    len_w = len(w)
    l = 0
    while d**l < len_w:
        l += 1
    w_padded = vector(Fq, list(w) + [Fq(0)]*(d**l - len_w))
    return w_padded, l


def lde(f: vector, d: int = 2):
    f_pad, l = pad_vec_to_d_exp(f)

    # Prepare for MLE/sumcheck for fields
    P = PolynomialRing(Fq, [f"x{i}" for i in range(l)])
    xs = P.gens()

    #   \sum_{w \in {0,1}^l} eq(x, w) * f(w)
    # = \sum_{w \in {0,1}^l} {\prod_{i \in l} ((1-x_i)(1-w_i) + x_i*w_i)} f(w)
    tilde_f = 0
    for w in range(0, d**l):
        # i \in {0,1}^l. e.g. i = (b_1, .., b_l)
        eq = 1
        # eq(x, w): are x and w the same for all bits?
        # starting from bit 0 to bit l-1
        # MSB: we map b_0*2^2 + b_1*1 to (b_0, b_1)
        for i in range(l):
            w_i = (w >> (l-i-1)) & 1
            eq *= (1 - xs[i]) * (1 - w_i) + xs[i] * w_i
        tilde_f += eq * f_pad[w]
    return tilde_f, xs


def test_mle():

    # Test: mle of [1,2,8,10]
    t = [1,2,8,10]

    mle_t, xs = lde(t)
    l = len(xs)

    # check with known result
    assert mle_t == 1 * (1-xs[0])*(1-xs[1]) + 2 * (1-xs[0])*xs[1] + 8 * xs[0] * (1-xs[1]) + 10 * xs[0] * xs[1]
    # \sum_z mle[t][z] = sum(t)
    sum_mle = sum([
        mle_t.subs({xs[i]: b[i] for i in range(l)})
        for b in itertools.product(range(2), repeat=l)
    ])
    assert sum_mle == sum(t)


def test_sumcheck():
    w = vector(Fq, [1, 35, 3, 9, 27])
    res = sumcheck(w)
    assert res == sum(w), f"result from sumcheck doesn't match sum directly: {w=}, {res=}"


def tests():
    test_mle()
    test_sumcheck()

if __name__ == '__main__':
    tests()
