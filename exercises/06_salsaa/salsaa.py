
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
from ntt import ntt, intt

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

print(f"SALSAA setup: Rq = Z_{q}[X]/(X^{d}+1), kappa={kappa}, m={m}, beta={beta}")



def main():
    # 0. Use R1CS witness from ../13_groth16/groth16.py

    # --- R1CS example from Groth16 ---
    # Arithmetic circuit: x^3 + x + 5 = out
    # Quadratic constraints:
    #   1. ( x ) * ( x ) = v_1
    #   2. (v_1) * ( x ) = v_2
    #   3. (v_2 + x + 5) * (1) = out

    # w = (1, public I/O, private inputs, intermediates)
    #              (1, out, x, v1, v2)
    w = vector(Fq, [1, 35, 3, 9, 27])
    n = len(w)
    L = matrix(Fq, [
        [0, 0, 1, 0, 0],
        [0, 0, 0, 1, 0],
        [5, 0, 1, 0, 1]
    ])
    R_mat = matrix(Fq, [
        [0, 0, 1, 0, 0],
        [0, 0, 1, 0, 0],
        [1, 0, 0, 0, 0],
    ])
    O_mat = matrix(Fq, [
        [0, 0, 0, 1, 0],
        [0, 0, 0, 0, 1],
        [0, 1, 0, 0, 0],
    ])
    num_constraints = L.nrows()
    Lw = L * w
    Rw = R_mat * w
    Ow = O_mat * w
    for i in range(num_constraints):
        assert Lw[i] * Rw[i] == Ow[i], f"constraint {i} doesn't hold, {Lw[i]=}, {Rw[i]=}, {Ow[i]=}"
    print(f"R1CS satisfied: {num_constraints} constraints, {n} variables")


    # ============================================================
    # Sumcheck protocol
    # ============================================================

    # TODO: V still needs to check `a ?= f(r_0, ..., r_{l-1})`!
    a, rs = sumcheck(w, d=2)
    print(f"{a=}, {rs=}")


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
    # 4. CCS (Customizable Constraint System)
    # ============================================================

    # TODO: convert R1CS to CCS


    # TODO: fold two CCS instances

    # # Verify imports work
    # A = ajtai_setup(kappa, m, Rq)
    # z = _gen_random_low_norm_witness(Rq, m, beta)
    # c = ajtai_commit(A, z)
    # # print(f"Ajtai commit OK: {c=}")


def lde_poly(w: vector, d: int):
    """
    Build LDE[w] as a symbolic multivariate polynomial over [d]^l using Lagrange basis.
    Returns (polynomial, xs) where xs are the symbolic variables.
    """
    # Pad w to size=d^l
    w_pad, l = pad_vec_to_d_exp(w, d)

    P = PolynomialRing(Fq, [f"x{i}" for i in range(l)])
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
        for k in range(d):
            lagrange = prod([
                (xs[j] - k_prime) / (Fq(k) - k_prime)
                for k_prime in range(d) if k_prime != k
            ])
            xi_eq_k.append(lagrange)
        eqs.append(xi_eq_k)

    # Do tensor product on our own

    # eqs_flat = [eq(x, (0,..., 0)), eq(x, (d-1, ..., d-1)]
    eqs_flat = []
    for idx in itertools.product(range(d), repeat=l):
        # e.g. idx = (0, ..., 0)
        # eq(x, (0,..., 0)) = eq(x_0, 0) * ... * eq(x_{l-1}, 0)
        eq = prod([eqs[j][idx[j]] for j in range(l)])
        eqs_flat.append(eq)

    # <w, tensor(r)>
    # res = <w, eqs_flat>
    #     = w_0*eq(x, (0,..., 0)) + ... + w_{d^l-1}*eq(x, (d-1, ..., d-1)
    tilde_f = sum(w_pad[i] * eqs_flat[i] for i in range(d**l))
    return tilde_f, xs


def sum_over_hypercube(poly, xs, fixed: dict, d: int, start: int, end: int):
    """Substitute fixed values, then sum over xs[start:end] \in [d]^(end-start)."""
    result = 0
    for b in itertools.product(range(d), repeat=end-start):
        subs = {**fixed, **{xs[start+i]: b[i] for i in range(end-start)}}
        result += poly.subs(subs)
    return result


def sumcheck(f: vector, d: int):
    mle_t, xs = lde_poly(f, d)
    # Claim: \sum_{b_0}...\sum_{b_{l-1}} f(b_0, ..., b_{l-1}) = a_j
    # P calculate the first a_j (a_0) and send it to V
    l = len(xs)
    # a_j, a at loop j
    # sum over [d]^l
    a = sum_over_hypercube(mle_t, xs, {}, d, start=0, end=l)

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
    # NOTE: After the loop, V has `a` (a_l = g_{l-1}(r_{l-1}))
    # Since V is not sure if g_{l-1}(x) = f(r_0,...,r_{l-2}, x),
    # V still needs to check `a ?= f(r_0, ..., r_{l-1})`
    return a, received_randoms


# Pad w to size d^l to make it on hypercube
def pad_vec_to_d_exp(w: vector, d: int):
    # First padding w to 8 elements so we have [d]^3
    len_w = len(w)
    l = 0
    while d**l < len_w:
        l += 1
    w_padded = vector(Fq, list(w) + [Fq(0)]*(d**l - len_w))
    return w_padded, l


def test_lde_poly():
    # Test with multiple d values on the same w
    w = vector(Fq, [1, 2, 8, 10, 3, 5])

    for d in [2, 3, 4]:
        w_pad, l = pad_vec_to_d_exp(w, d)
        poly, xs = lde_poly(w, d)

        # Every point in [d]^l should recover the padded value
        for bit_repr in itertools.product(range(d), repeat=l):
            bit_repr_reversed = bit_repr[::-1]
            idx = sum([bit * d**i for i, bit in enumerate(bit_repr_reversed)])
            val = poly.subs({xs[i]: bit_repr[i] for i in range(l)})
            assert w_pad[idx] == val, \
                f"lde_poly mismatch at {bit_repr=}, {d=}: expected {w_pad[idx]}, got {val}"

    # Test: size already equals d^l (no padding needed)
    w_exact = vector(Fq, [1, 2, 3, 4])
    for d in [2, 4]:
        w_pad, l = pad_vec_to_d_exp(w_exact, d)
        poly, xs = lde_poly(w_exact, d)
        assert len(w_pad) == len(w_exact) or w_pad[len(w_exact):] == vector(Fq, [0]*(len(w_pad)-len(w_exact)))
        for bit_repr in itertools.product(range(d), repeat=l):
            bit_repr_reversed = bit_repr[::-1]
            idx = sum([bit * d**i for i, bit in enumerate(bit_repr_reversed)])
            assert w_pad[idx] == poly.subs({xs[i]: bit_repr[i] for i in range(l)})

    # Test: single element
    w_single = vector(Fq, [7])
    for d in [2, 3]:
        poly, xs = lde_poly(w_single, d)
        if len(xs) == 0:
            assert poly == Fq(7)
        else:
            assert poly.subs({xs[0]: 0}) == Fq(7)


def tests():
    test_lde_poly()

if __name__ == '__main__':
    tests()
    main()
