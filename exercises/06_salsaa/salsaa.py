
import random as rand
import itertools
import sys
from pathlib import Path
sys.path.append(str(Path(__file__).resolve().parent.parent / "04_ajtai"))
sys.path.append(str(Path(__file__).resolve().parent.parent / "02_ntt"))

from sage.all import *
from ajtai import setup as ajtai_setup, commit as ajtai_commit
from ajtai import _gen_random_low_norm_poly, _gen_random_low_norm_witness
from ajtai import _l_inf_norm_Rq, _l_inf_norm_vec
from ntt import ntt, intt, find_root_of_unity, find_generator_from_mul_group

# ============================================================
# SALSAA — Lattice-based Folding Scheme
# ============================================================

# Ring setup
# R_q = Z_q[X]/(X^d+1)
q = 17
Fq = GF(q)
d = 4
R = PolynomialRing(Fq, 'X')
X = R.gen()
Rq = R.quotient(X**d + 1, 'x')
x = Rq.gen()

# Ajtai parameters
kappa = 2
m = 4
beta = 1

# CRT slots degree, i.e. we split R_q into d/e fields.
# In NTT we fully split so e = 1.
# Seems salsaa is using 2 in section E in p.36. But probably we can
# use 1 for simplicity.
e = 1

# conductor
f = 2 * d

# H: \hat n \times n
# F: n \times m
# W: m \times r
# Y: \hat n \times r
# HFW = Y
m = 2
r = 2
# The d in [d]^l for LDE, to avoid reusing `d` from X^d+1
D = 2

print(f"SALSAA setup: Rq = Z_{q}[X]/(X^{d}+1), kappa={kappa}, m={m}, beta={beta}")


def get_W(m: int, r: int):
    # Hardcoded for debugging
    W = matrix(Rq, [
        [16*x**2 + 16*x + 1,       x**2 + 16*x + 1],
        [16*x**3 + 16*x + 1,           16*x**2 + 1],
    ])
    assert W.nrows() == m
    assert W.ncols() == r
    return W
    # MS = MatrixSpace(Rq, m, r)
    # W = MS([ _gen_random_low_norm_poly(Rq, beta) for _ in range(r * m) ])
    # return W


def get_u_vec(Fq, r: int, d: int, e: int):
    u = Fq(rand.randint(1, q-1))
    u_vec = vector(Fq, [u ** i for i in range(r*d//e)])
    return u_vec

def conjugate(r):
    return Rq(r.lift()(-x**(d-1)))

def rok_bar_sum(H, F, Y, t, W):
    # W = [w_0, ..., w_{r-1}], w_i \in R^m

    # We handle one column a time, and concatenate all of them in the end.
    def calculate_CRT_LDE_w_LDE_bar_w(col_idx: int):
        w = W.column(col_idx)
        print(f"{w=}")
        w_bar = vector(Rq, [conjugate(w_i) for w_i in w])
        print(f"{w_bar=}")
        lde_w, xs = lde_poly(w, D=D)
        lde_w_bar, _ = lde_poly(w_bar, D=D)
        print(f"{lde_w=}")
        print(f"{lde_w_bar=}")
        lde_w_times_w_bar = lde_w * lde_w_bar
        print(f"{lde_w_times_w_bar=}")

        # \tilde f (x_0, x_1) = A + B x_0 + C x_1 + D x_0 x_1  , A, B, C, D \in R_q
        coeff_Rqs = lde_w_times_w_bar.coefficients()
        #   [NTT(A), NTT(B), NTT(C), NTT(D)]
        # = [(a_0, a_1, a_2, a_3), (b_0, b_1, b_2, b_3), (c_0, c_1, c_2, c_3), (d_0, d_1, d_2, d_3)]
        coeff_Ntts = [ntt(rq.list(), Rq) for rq in coeff_Rqs]

        # \tilde f_\text{slot_0} (x_0, x_1) = a_0 + b_0 x_0 + c_0 x_1 + d_0 x_0 x_1
        # ...
        # \tilde f_\text{slot_3} (x_0, x_1) = a_3 + b_3 x_0 + c_3 x_1 + d_3 x_0 x_1
        # tilde_fs = [\tilde f_\text{slot_0}, ..., \tilde f_\text{slot_{d//e-1}}]
        tilde_fs = []
        for i in range(d // e):
            tilde_f_slot_i = sum([
                coeff_Ntts[j][i] * m
                for j, m in enumerate(lde_w_times_w_bar.monomials())
            ])
            tilde_fs.append(tilde_f_slot_i)
        return tilde_fs, xs

    # Get CRT(LDE[w]*LDE[\bar w]) for each column, and concatenate all of them into one list
    crt_LDE_W_LDE_bar_W = []
    for i in range(W.ncols()):
        tilde_fs_i, xs = calculate_CRT_LDE_w_LDE_bar_w(i)
        crt_LDE_W_LDE_bar_W.extend(tilde_fs_i)
    # xs is reused from last iteration, they're all the same
    assert len(crt_LDE_W_LDE_bar_W) != 0, "empty W and thus xs is not existing!"

    # u^T: [u^0, u^1, ..., u^{rd/e}
    u = get_u_vec(Fq, r, d, e)

    tilde_f = sum([
        u_i * tilde_f_i
        for u_i, tilde_f_i in zip(u, crt_LDE_W_LDE_bar_W)
    ])
    print(f"{tilde_f=}")

    res, rands = sumcheck(tilde_f, xs, D)
    print(f"{res=}, {rands=}")

    # TODO: verify result from sumcheck
    #   1. get r^T from CRT
    #   2. check a_\mu ?= u^T * CRT(s_0 \dot \bar s_1)

    # TODO: check the lde-tensor relation.

    return

def rok_norm(H, F, Y, v, W):


    #
    # Prover
    #

    # Since P has witness, just sum the square of all coeffs in the poly.
    def get_norm_square_Rq(w: vector):
        """
        w: a column of W, \in R^m
        """
        res = 0
        for w_i in w:
            res += sum([c*c for c in w_i.list()])
        return res

    # t = (<w_i, \bar w_i>)_{i \in [r]}
    t = [get_norm_square_Rq(W.column(i)) for i in range(r)]
    print(f"{t=}")

    # Sends t to V

    #
    # Verifier
    #
    # TODO: check Trace(t_i) <= v^2 \forall i \in [r]

    # P and V go on to rok \bar sum
    rok_bar_sum(H, F, Y, t, W)


def main():
    beta = 1
    # it's l_\infty though
    v = None

    #
    # Norm check
    #

    # FIXME: Mock H, F, Y, v for now.
    # Should be fixed later when rok_sum_bar is working
    H, F, Y, v = (None, None, None, None)
    # W \in R_q^{m*r}
    W = get_W(m, r)
    # print(f"{W=}")

    # TODO: V still needs to check `a ?= f(r_0, ..., r_{l-1})`!
    # a, rs = sumcheck(w, d=2)
    # print(f"{a=}, {rs=}")
    rok_norm(H, F, Y, v, W)


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


def lde_poly(w: vector, D: int):
    """
    Build LDE[w] as a symbolic multivariate polynomial over [D]^l using Lagrange basis.
    Returns (polynomial, xs) where xs are the symbolic variables.
    """
    # Pad w to size=d^l
    w_pad, l = pad_vec_to_d_exp(w, D)

    P = PolynomialRing(w.base_ring(), [f"x{i}" for i in range(l)])
    xs = P.gens()

    #
    # Prepare for eqs
    # [
    #   [eq(x_0, 0), ..., eq(x_0, d-1)],
    #   ...,
    #   [eq(x_{l-1}, 0), ..., eq(x_{l-1}, d-1)],
    # ]
    eqs = []
    for j in range(l):
        # in the end of this loop
        # xi_eq_k = (eq(x_j, 0), ..., eq(x_j, d-1))
        xi_eq_k = []
        for k in range(D):
            lagrange = prod([
                (xs[j] - k_prime) / (Fq(k) - k_prime)
                for k_prime in range(D) if k_prime != k
            ])
            xi_eq_k.append(lagrange)
        eqs.append(xi_eq_k)

    # Do tensor product on our own

    # eqs_flat = [eq(x, (0,..., 0)), eq(x, (d-1, ..., d-1)]
    eqs_flat = []
    for idx in itertools.product(range(D), repeat=l):
        # e.g. idx = (0, ..., 0)
        # eq(x, (0,..., 0)) = eq(x_0, 0) * ... * eq(x_{l-1}, 0)
        eq = prod([eqs[j][idx[j]] for j in range(l)])
        eqs_flat.append(eq)

    # <w, tensor(r)>
    # res = <w, eqs_flat>
    #     = w_0*eq(x, (0,..., 0)) + ... + w_{d^l-1}*eq(x, (d-1, ..., d-1)
    tilde_f = sum(w_pad[i] * eqs_flat[i] for i in range(D**l))
    return tilde_f, xs


def sum_over_hypercube(poly, xs, fixed: dict, D: int, start: int, end: int):
    """Substitute fixed values, then sum over xs[start:end] \in [d]^(end-start)."""
    result = 0
    for b in itertools.product(range(D), repeat=end-start):
        subs = {**fixed, **{xs[start+i]: b[i] for i in range(end-start)}}
        result += poly.subs(subs)
    return result


def sumcheck(lde_poly, xs, D: int) -> tuple[Fq, list[Fq]]:
    # Claim: \sum_{b_0}...\sum_{b_{l-1}} f(b_0, ..., b_{l-1}) = a_j
    # P calculate the first a_j (a_0) and send it to V
    l = len(xs)
    # a_j, a at loop j
    # sum over [d]^l
    a = sum_over_hypercube(lde_poly, xs, {}, D, start=0, end=l)

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


# Pad w to size d^l to make it on hypercube
def pad_vec_to_d_exp(w: vector, D: int):
    # First padding w to 8 elements so we have [d]^3
    len_w = len(w)
    l = 0
    while D**l < len_w:
        l += 1
    w_padded = vector(w.base_ring(), list(w) + [Fq(0)]*(D**l - len_w))
    return w_padded, l


def test_lde_poly():
    # Test with multiple d values on the same w
    w = vector(Fq, [1, 2, 8, 10, 3, 5])

    for D in [2, 3, 4]:
        w_pad, l = pad_vec_to_d_exp(w, D)
        poly, xs = lde_poly(w, D)

        # Every point in [d]^l should recover the padded value
        for bit_repr in itertools.product(range(D), repeat=l):
            bit_repr_reversed = bit_repr[::-1]
            idx = sum([bit * D**i for i, bit in enumerate(bit_repr_reversed)])
            val = poly.subs({xs[i]: bit_repr[i] for i in range(l)})
            assert w_pad[idx] == val, \
                f"lde_poly mismatch at {bit_repr=}, {D=}: expected {w_pad[idx]}, got {val}"

    # Test: size already equals d^l (no padding needed)
    w_exact = vector(Fq, [1, 2, 3, 4])
    for D in [2, 4]:
        w_pad, l = pad_vec_to_d_exp(w_exact, D)
        poly, xs = lde_poly(w_exact, D)
        assert len(w_pad) == len(w_exact) or w_pad[len(w_exact):] == vector(Fq, [0]*(len(w_pad)-len(w_exact)))
        for bit_repr in itertools.product(range(D), repeat=l):
            bit_repr_reversed = bit_repr[::-1]
            idx = sum([bit * D**i for i, bit in enumerate(bit_repr_reversed)])
            assert w_pad[idx] == poly.subs({xs[i]: bit_repr[i] for i in range(l)})

    # Test: single element
    w_single = vector(Fq, [7])
    for D in [2, 3]:
        poly, xs = lde_poly(w_single, D)
        if len(xs) == 0:
            assert poly == Fq(7)
        else:
            assert poly.subs({xs[0]: 0}) == Fq(7)


def tests():
    test_lde_poly()

if __name__ == '__main__':
    tests()
    main()
