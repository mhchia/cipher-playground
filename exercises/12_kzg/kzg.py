from sage.all import *
import sys
from pathlib import Path
sys.path.append(str(Path(__file__).resolve().parent.parent / "11_pairings"))
from pairing import find_G1, find_G2

# ============================================================
# Pairing setup (same as Groth16)
# E: y^2 = x^3 + 1 over F_257, embedding degree k=2, r=43
# ============================================================

p = 257
F = GF(p)
E = EllipticCurve(F, [0, 1])
n = E.order()
r = factor(n)[-1][0]
assert r == 43

F2 = GF(p**2, 'u')
E2 = EllipticCurve(F2, [0, 1])
n2 = E2.order()
assert n2 % (r**2) == 0
assert (p**2 - 1) % r == 0

G1 = find_G1(E, r)
G2 = find_G2(E, E2, r)
G1_lifted = E2(G1)

# sanity check: non-trivial pairing
e_g = G1_lifted.weil_pairing(G2, r)
assert G1.order() == r and G2.order() == r
assert e_g**r == 1 and e_g != 1
print(f"Pairing setup OK: {p=}, {r=}, {G1=}, {G2=}")

# Fr: the scalar field for all circuit arithmetic
# The size of F_r is a tradeoff between **DLP hardness**
# and **execution efficiency**.
Fr = GF(r)
print(f"{Fr=}")

# ============================================================
# Trusted setup (SRS from random tau)
# ============================================================

_R = PolynomialRing(Fr, 'X')
X = _R.gen()

# Set d as the max degree of polynomials we can commit to
# This affects
d = 10
# tau must not be 0
tau = Fr.random_element()
while tau == 0:
    tau = Fr.random_element()

# Since we verify f(\tau)-v in G1, we need up to \tau^d G1
# For G2, we need at most \tau-z so it's up to \tau G2
# SRS:
#   [G1, \tau G1, ..., \tau^d G1]
#   [G2, \tau G2]

SRS_G1 = [(tau**i)* G1 for i in range(d+1)]
SRS_G2 = [(tau**i)* G2 for i in range(2)]

def eval_on_tau(poly: _R, srs) -> _R:
    """
    Evaluate polynomial at tau in the exponent: f(tau)*G = sum(f_i * [tau^i]) through throguh
    addition homomorphism
    """
    coeff = poly.list()
    assert len(coeff) <= len(srs), f"SRS too short: {len(coeff)=}, {len(srs)=}"
    return sum(c * s for c, s in zip(coeff, srs[:len(coeff)]))


print(f"SRS ready: {d=}, {tau=}")

# ============================================================
# Q1: KZG Commit + Open + Verify
# ============================================================

# Define f(x) = 2x^2 + 3x + 1 over Fr
f = 2*X**2 + 3*X + 1

# 1.1. P commits to f
# C = [f(\tau)]_1
C = eval_on_tau(f, SRS_G1)

# 1.2. V sends z to P
z = 5

# 1.3. P proves f(z) = v by sending (\pi, v)
def prove(_z: int) -> tuple[G1, int]:
    v = f(_z)
    # f(z) = v <-> f(x) - v = q(x)(x-z)
    _q, remainder = (f-v).quo_rem(X-_z)
    if remainder != 0:
        raise Exception("this should never happen")
    # \pi=[q(\tau)]_1
    _pi = eval_on_tau(_q, SRS_G1)
    return _pi, v

# P proves
pi, v = prove(z)

# 1.4. V verifies \pi is valid or not with pairing
def verify(_C: int, _pi: G1, z: int, v: int) -> bool:
    """
    we know f(z) = v <-> f(x) - v = q(x)(x-z)
    With commitment and pairing, we get
    -> e(com(f)-com(v), 1) ?= e(com(q), com(x-z))
    -> e(C-[v]_1, G2) ?= e(\pi, [tau-z]_2)
    """
    l = _C - eval_on_tau(v * X**0, SRS_G1)
    lhs = E2(l).weil_pairing(
        G2,
        r,
    )
    rhs = E2(_pi).weil_pairing(
        eval_on_tau(X - z, SRS_G2),
        r,
    )
    print(f"verify: {lhs=}, {rhs=}")
    return lhs == rhs


# V get \pi and v, verifies with _C from commitment and SRS
is_proof_valid = verify(C, pi, z, v)
assert is_proof_valid, f"proof is invalid, {C=}, {pi=}, {z=}, {v=}"


# ============================================================
# Verify failure case (prover lies about v)
# ============================================================

wrong_v = v + 1

assert not verify(C, pi, z, wrong_v)


# ============================================================
# Batching
# ============================================================

# Approach 1: opening 2 points need 4 pairing. we can use RLC on V side to decrease to 3 pairing
# (1+r) com(f) - com(v1+ r v2) = com(q1(x)*(x-z1) + r q2(x) * (x-z2))
# -> e((1+r)com(f)-com(v1+r v2), 1) = e(com(q1), com(x-z1)) * e(r*com(q2), com(x-z2))


# Approach 2: opening k points, use vanishing polynomial.
#
# P wants to open z_k, ..., z_k at once.
# So f(z_1) = v_1, ..., f(z_k) = v_k
# Let I be the lagrange interpolation polynomial with [(z_1,v_1), ..., (z_k,v_k)]
# Then f(x) - I(x) = 0 forall x \in [z_1, ..., z_k], i.e. they are roots of f(x) - I(x)
# let t(x) = (x-z_1)...(x-z_k) be the vaninishing polynomial
# we know f(x) - I(x) = q(x) t(x)
# -> com(f-I) = com(q*t)
# -> e(com(f)-com(I), 1) = e(com(q), com(t))

SRS_G2_batch = [(tau**i)* G2 for i in range(d+1)]

# V picks k points
k = 2
zs = list(range(1, k+1))

# P calculates I, t, and get \pi and v1...vk
def batch_prove(_zs: list[int]) -> tuple[G1, list[int]]:
    vs = [f(z) for z in _zs]
    # vanishing polynomial
    t = prod([X-z for z in _zs])
    I = _R.lagrange_polynomial([(z, f(z)) for z in _zs])

    # f(z) = v <-> f(x) - v = q(x)(x-z)
    _q, remainder = (f-I).quo_rem(t)
    if remainder != 0:
        raise Exception("this should never happen")
    # \pi=[q(\tau)]_1
    print(f"{_q=}")
    _pi = eval_on_tau(_q, SRS_G1)
    return _pi, vs

# P proves with k points and generate single proof \pi
pi_batch, vs = batch_prove(zs)
print(f"{pi=}, {vs=}")

def batch_verify(_C: G1, _pi: G1, _zs: list[int], _vs: list[int]) -> bool:
    # vanishing polynomial
    t = prod([X-z for z in _zs])
    I = _R.lagrange_polynomial(list(zip(_zs, _vs)))
    l = _C - eval_on_tau(I, SRS_G1)
    lhs = E2(l).weil_pairing(
        G2,
        r,
    )
    rhs = E2(_pi).weil_pairing(
        eval_on_tau(t, SRS_G2_batch),
        r,
    )
    print(f"batch_verify: {lhs=}, {rhs=}")
    return lhs == rhs


# V verifies with single \pi
# But, here V needs to interpolate I and calculate vanishing polynomial t
# time becomes O(d^2)?
is_all_points_valid = batch_verify(C, pi_batch, zs, vs)
assert is_all_points_valid, "batch opening failed"
