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
q = 17
d = 4
Fq = GF(q)
R = PolynomialRing(Fq, 'X')
X = R.gen()
Rq = R.quotient(X**d + 1, 'x')
x = Rq.gen()

# Ajtai parameters
kappa = 2
m = 4
beta = 1

print(f"SALSAA setup: Rq = Z_{q}[X]/(X^{d}+1), kappa={kappa}, m={m}, beta={beta}")

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


# Practice a bit of MLE
# First padding w to 8 elements so we have {0,1}^3
l = 0
while 2**l < len(w):
    l += 1
w_padded = vector(Fq, list(w) + [Fq(0)]*(2**l - len(w)))
print(f"{w_padded=}")

P = PolynomialRing(Fq, [f"x{i}" for i in range(l)])
xs = P.gens()
print(f"{xs=}")


def mle(f: vector, l: int):
    # TODO: check len(f) is 2**l
    #   \sum_{w \in {0,1}^l} eq(x, w) * f(w)
    # = \sum_{w \in {0,1}^l} {\prod_{i \in l} ((1-x_i)(1-w_i) + x_i*w_i)} f(w)
    tilde_f = 0
    for w in range(0, 2**l):
        # i \in {0,1}^l. e.g. i = (b_1, .., b_l)
        eq = 1
        # eq(x, w): are x and w the same for all bits?
        # starting from bit 0 to bit l-1
        # MSB: we map b_0*2^2 + b_1*1 to (b_0, b_1)
        for i in range(l):
            w_i = (w >> (l-i-1)) & 1
            eq *= (1 - xs[i]) * (1 - w_i) + xs[i] * w_i
        tilde_f += eq * f[w]
    return tilde_f


# Test: mle of [1,2,8,10]
t = [1,2,8,10]
l_t = 2
mle_t = mle(t, l_t)
assert mle_t == 1 * (1-xs[0])*(1-xs[1]) + 2 * (1-xs[0])*xs[1] + 8 * xs[0] * (1-xs[1]) + 10 * xs[0] * xs[1]
print(f"{mle_t=}")
c = sum(t)


def sumcheck(mle, xs, l: int, d: int = 2):
    import random
    from itertools import product
    # Claim: \sum_{b_0}...\sum_{b_{l-1}} f(b_0, ..., b_{l-1}) = a_j
    # P calculate the first a_j (a_0) and send it to V
    a = 0
    for b in product([0,1], repeat=l):
        a += mle.subs({xs[i]: b[i] for i in range(l)})

    received_randoms = []
    for j in range(l):
        #
        # Prover
        #
        # Let the variate i be X and sum the rest of the variates over hypercube.
        # let h_0(x) = \sum_{b_1} ... \sum_{b_{l-1}} f(x, b_1, ..., b_{l-1})
        # P calculate h_0(x) and send it to V as g_0(x)
        rs = {
            xs[i]: v for i, v in enumerate(received_randoms)
        }
        # Sum of all g evaluation from position j+1 to l-1
        # g_j(x) = f(r_0, ..., r_{j-1}, x, 0, ..., 0) + f(r_0, ..., r_{j-1}, x, 0, ..., 1) +...
        g = 0
        for b in product([0,1], repeat=l-1-j):
            bs = {
                xs[i+j+1]: b[i] for i in range(l-1-j)
            }
            g += mle.subs({**rs, **bs})

        # Send `g_j` to V

        #
        # Verifier
        #

        # V is not sure if g_0(x) = h_0(x) as P claimed
        #   and needs to verify by
        #   1. c = g_0(0) + g_0(1)
        #   2. g_j(r) ?= \sum_{b_2} ... \sum_{b_{l-1}} f(r, b_2, ..., b_{l-1}), by SZPL
        #       - this is done by P running sumcheck again

        # Step 1, c ?= g_0(0) + g_0(1)
        assert a == sum([g.subs({xs[j]: i}) for i in range(d)])

        # Step 2. g_j(r) ?= \sum_{b_{j+1}} ... \sum_{b_{l-1}} f(r_0,...,r_{j-1}, r, b_{j+1}, ..., b_{l-1})
        r = random.randint(0, q-1)
        a = g.subs({xs[j]: r})
        # Send r to P

        # P: store r for next round
        received_randoms.append(r)
    # NOTE: After the loop, V has `a` (a_l = g_{l-1}(r_{l-1}))
    # Since V is not sure if g_{l-1}(x) = f(r_0,...,r_{l-2}, x),
    # V still needs to check `a ?= f(r_0, ..., r_{l-1})`
    return a

# TODO: V still needs to check `a ?= f(r_0, ..., r_{l-1})`!
a = sumcheck(mle_t, xs, l_t)
print(f"{a=}")


# mle_w = mle(w_padded, l)
# print(f"{mle_w=}")




# ============================================================
# 1. Challenge spaces
# ============================================================

# TODO: define C_small, sample challenges


# ============================================================
# 2. MLE (Multilinear Extension)
# ============================================================

# TODO: construct MLE from truth table, evaluate


# ============================================================
# 3. Sumcheck protocol
# ============================================================

# TODO: implement sumcheck over F_q, then over R_q


# ============================================================
# 4. CCS (Customizable Constraint System)
# ============================================================

# TODO: convert R1CS to CCS


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

# TODO: fold two CCS instances

def main():

    # Verify imports work
    A = ajtai_setup(kappa, m, Rq)
    z = _gen_random_low_norm_witness(Rq, m, beta)
    c = ajtai_commit(A, z)
    # print(f"Ajtai commit OK: {c=}")


if __name__ == '__main__':
    main()
