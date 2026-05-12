"""
Ring setup and protocol parameters for SALSAA.

R_q = Z_q[X]/(X^d + 1)  with the Galois conjugate σ_{-1}: X ↦ X^{-1} = -X^{d-1}.
"""
import sys
from pathlib import Path
sys.path.append(str(Path(__file__).resolve().parent.parent / "04_ajtai"))

from sage.all import *
from ajtai import setup as ajtai_setup, commit as ajtai_commit
from ajtai import _gen_random_low_norm_poly, _gen_random_low_norm_witness
from ajtai import _l_inf_norm_Rq, _l_inf_norm_vec


# ============================================================
# Ring: R_q = Z_q[X]/(X^d+1)
# ============================================================
q = 17
Fq = GF(q)
d = 4
R = PolynomialRing(Fq, 'X')
X = R.gen()
Rq = R.quotient(X**d + 1, 'x')
x = Rq.gen()

# CRT slots degree, i.e. we split R_q into d/e fields.
# In NTT we fully split so e = 1.
# Seems salsaa is using 2 in section E in p.36. But probably we can
# use 1 for simplicity.
e = 1

# conductor
f = 2 * d


# ============================================================
# Protocol parameters
# ============================================================

# Ajtai parameters
kappa = 2
beta = 1

# H: \hat n \times n
# F: n \times m
# W: m \times r
# Y: \hat n \times r
# HFW = Y
n_hat = 2
n = 2
m = 4
r = 2

# The d in [d]^l for LDE, to avoid reusing `d` from X^d+1
D = 2


print(f"SALSAA setup: Rq = Z_{q}[X]/(X^{d}+1), kappa={kappa}, m={m}, beta={beta}")


# ============================================================
# Ring helpers
# ============================================================

def conjugate(r):
    return Rq(r.lift()(-x**(d-1)))



def to_centered(c):
    """
    Represent Z_q as [-\ceil{q/2}, \floor{q/2}]
    """
    v = int(c)
    return v - q if v > q // 2 else v


def get_l2_norm_square(w: vector):
    l2_norm_s = 0
    for w_i in w:
        coeffs = w_i.list()
        l2_norm_s += sum([c*c for c in coeffs])
    return l2_norm_s
