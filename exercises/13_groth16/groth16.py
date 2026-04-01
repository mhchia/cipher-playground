from sage.all import *
import sys
from pathlib import Path
sys.path.append(str(Path(__file__).resolve().parent.parent / "11_pairings"))
from pairing import find_G1, find_G2

# ============================================================
# Pairing setup
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
Fr = GF(r)

# =======
# Prover proves they know x s.t. x^3 + x + 5 = out, where out is
# a public output.
#

# ============================================================
# Step 1: Arithmetic circuit
# x^3 + x + 5 = out
# ============================================================

# A: turn into quadratic expression
#   1. ( x ) * ( x ) = v_1
#   2. (v_1) * ( x ) = v_2
#   3. (v_2 + x + 5) * (1) = out


# ============================================================
# Step 2: R1CS (L, R, O matrices + witness vector w)
# ============================================================

# w = (1, public I/O, private inputs, intermediates)
#              (1, out, x, v1, v2)
w = vector(Fr, [1, 35, 3, 9, 27])
n = len(w)
L = matrix(Fr, [
    [0, 0, 1, 0, 0],
    [0, 0, 0, 1, 0],
    [5, 0, 1, 0, 1]
])
m = L.nrows()
R = matrix(Fr, [
    [0, 0, 1, 0, 0],
    [0, 0, 1, 0, 0],
    [1, 0, 0, 0, 0],
])
O = matrix(Fr, [
    [0, 0, 0, 1, 0],
    [0, 0, 0, 0, 1],
    [0, 1, 0, 0, 0],
])
Lw = L * w
Rw = R * w
Ow = O * w
for i in range(m):
    assert Lw[i] * Rw[i] == Ow[i], f"constraint {i} doesn't hold, {Lw[i]=}, {Rw[i]=}, {Ow[i]=}"


# ============================================================
# Step 3: QAP (Lagrange interpolation -> L(x), R(x), O(x), h(x))
# ============================================================

_R = PolynomialRing(Fr, 'X')
X = _R.gen()

# 1. find roots of unity: {\omega^0, \omega^1, ...}
# 1.1. recall: F_r is mul group = {1, ..., r-1}, it's cyclic since r is prime.
# we need to find a generator g first, and we know g^{r-1} = 1
g = Fr.multiplicative_generator()
assert g.multiplicative_order() == r-1

# 1.2. (Skipped) find primitive d-th roots of unity.
# TODO: Now we use normal points and lagrange interpolation to simplify implementation.
# Can revisit later. Set one point for each row
roots_x = list(map(Fr, range(m)))

# 2. t(x) = \prod_i (x-\omega^i)
tx = prod([(X - omega_i) for i, omega_i in enumerate(roots_x)])
print(f"{tx=}")

# 3. L(x) = \sum_j w_j * l_j(x), R(x), O(x)

# L(w^i) is Lw[i]
def get_poly(mat: matrix, w: vector) -> _R:
    # l_j(x) for all columns
    ljs = []
    for j in range(n):
        roots_xy = [(x, mat[idx][j]) for idx, x in enumerate(roots_x)]
        lj = _R.lagrange_polynomial(roots_xy)
        ljs.append(lj)
    # \sum_j w_j l_j(x)
    return sum([_w * lj for _w, lj in zip(w, ljs)])

Lx, Rx, Ox = get_poly(L, w), get_poly(R, w), get_poly(O, w)
# Sanity check: Lx(\omega^i) = Lw[i]
for i in range(m):
    omega_i = roots_x[i]
    assert Lx(omega_i) == Lw[i], f"Lx(\omega^i) != Lw[i] at {i=}, {omega_i=}, {Lx(omega_i)=}, {Lw[i]=}"

print(f"{Lx=}")
print(f"{Rx=}")
print(f"{Ox=}")

# 4. h(x) = (L * R - O)/t(x)
hx, remainder = (Lx * Rx - Ox).quo_rem(tx)
assert remainder == 0
print(f"{hx=}")

# ============================================================
# Step 4: Trusted setup (SRS from random tau)
# ============================================================

# NOTE: Here we implement with the no-blinding version, so it's not depending on the circuit,
# but depending only on the number of constraints.
# For real groth16, SRS contains [l_j(\tau)]_1, [r_j(\tau)]_1, [o_j(\tau)]_1, etc
# so it depends on circuits.

# FIXME: tau should be random and must not be one of the roots {1, w^1, w^2, ...}
# otherwise L(\tau)*R(\tau) - O(\tau) is always be 0.
tau = r - 1
# since deg(L)=deg(R)=deg(O) = m-1, deg(t) = m, we need SRS `m` long
# [1]_1, [\tau]_1, ..., [\tau^m]_1
tau_G1 = [(tau**i)*G1 for i in range(m+1)]
# [1]_2, [\tau]_2, ..., [\tau^m]_2
tau_G2 = [(tau**i)*G2 for i in range(m+1)]

# p(\tau)G = (p0+p1 \tau^1 + ...)G = p0 tau_G1[0] + p1 tau_G1[1]+...
def eval_on_tau(poly: _R, srs: list):
    coeff = poly.list()
    assert len(coeff) <= len(srs), f"SRS is not long enough: {len(coeff)=}, {len(srs)=}, {coeff=}"
    return sum([
        c_i * srs_i
        for c_i, srs_i in zip(
            coeff,
            srs[:len(coeff)]
        )
    ])

# Commit t(x) (evaluate t(\tau) on G2)
t_tau = eval_on_tau(tx, tau_G2)

# ============================================================
# Step 5: Prove (A, B, C from SRS)
# ============================================================
# P proves L, R, O, h (since these depends on witness w)
# Commit L(x) (evaluate L(\tau) on G1)
A = L_tau = eval_on_tau(Lx, tau_G1)
# Commit R(x) (evaluate R(\tau) on G2)
B = R_tau = eval_on_tau(Rx, tau_G2)
C = h_tau = eval_on_tau(hx, tau_G1)
O_tau = eval_on_tau(Ox, tau_G1)

# ============================================================
# Step 6: Verify (pairing check)
# ============================================================

# FIXME: this version V doesn't check if public inputs match the constraint.
# Should be further addressed.

# Recall: L(x) * R(x) - O(x) = h(x) * t(x)
# -> with pairing we verify: e(com(L), com(R)) = e(com(O), 1) * e(com(h), com(t))

# GT = G1_lifted.weil_pairing(G2, r)
lhs = E2(A).weil_pairing(B, r)
rhs = E2(O_tau).weil_pairing(G2, r) * E2(C).weil_pairing(t_tau, r)
assert lhs == rhs, f"pairing check mismatch: {lhs=}, {rhs=}"
