from sage.all import *

# ============================================================
# Pairing Exercises
# Bilinear map: e(P, Q) where P ∈ G1, Q ∈ G2, result ∈ GT
# Key property: e(aP, bQ) = e(P, Q)^{ab}
# ============================================================

def find_G1(E: EllipticCurve, r: int):
    n = E.order()
    if n % r != 0:
        raise ValueError(f"n must be divided by r: {n=}, {r=}")
    while True:
        P = E.random_point()
        h = n / r
        R = h*P
        print(f"{R.order()}")
        if R.order() == r:
            return R


def find_G2(E: EllipticCurve, E2: EllipticCurve, r: int) -> EllipticCurve:
    def is_in_G1(P: EllipticCurve) -> bool:
        if P.is_zero():
            return True
        try:
            # If P's coordinate fits in E, P must be the r subgroup in G1.
            # What we're looking for is G2 where G1 != G2
            E(P[0], P[1])
            return True
        except (TypeError, KeyError, ValueError):
            return False

    n = E2.order()
    r_square = r**2
    if n % r_square != 0:
        raise ValueError(f"n must be divided by r: {n=}, {r_square=}")
    while True:
        P = E2.random_point()
        if P.is_zero():
            continue
        R = (n / r_square)*P
        print(f"{R.order()=}")
        if R.order() == r and not is_in_G1(R):
            return R


#
# Find G_1 \subset E(F_p), G_2 \subset E(F_{p^2}), and GT \subset F_{p^k}
# where o(G_1) = o(G_2) = o(GT), s.t. e: G_1 \times G_2 \rightarrow GT.
# Step 1. r | |E(F_p)|       : so G_1 = prime order r subgroup
# Step 2. r^2 | |E(F_{p^2})| : so E(F_{p^2}) has a prime order r subgroup **other than** the lifted G_1.
# Step 3. r | p^k - 1        : k is the embedding degree
#
p = 59
# embedding deg = 2
k = 2
F = GF(p)
E = EllipticCurve(F, [1,0])  # y^2 = x^3+1x
n = E.order()
factors = factor(n)
# Step 1. r | |E(F_p)|       : so G_1 = prime order r subgroup
# get the largest prime factor r
r = factors[-1][0]
G1 = find_G1(E, r)

F2 = GF(p**2, 'u')
E2 = EllipticCurve(F2, [1,0])  # same curve but different field
n2 = E2.order()
print(f"|E(F2)|={n2}")

# Step 2. r^2 | |E(F_{p^2})| : so E(F_{p^2}) has a prime order r subgroup **other than** the lifted G_1.
assert n2 % (r** 2) == 0, f"r^2 does not divide |E(Fp^k)|, {r**2=}, {n2=}"
r2 = r**2
h2 = n2 / r2

G2 = find_G2(E, E2, r)
print(f"{G1.order()=}, {G2.order()=}, {G2=}")

# Step 3. r | p^k - 1        : k is the embedding degree
assert (p**k - 1) % r == 0
# Need to lift G1 into E(F_{p^2}) first, since both points need
# to be on the same curve object
G1_lifted = E2(G1)
# do pairing: e(G1, G2) = g, where g is the generator of GT.
GT = G1_lifted.weil_pairing(G2, r)
assert GT**r == 1, f"g is not a generator of GT, {GT=}, {GT**r=}"



# ============================================================
# Q1: Bilinearity verification (pen & paper, write answers here)
#
# Given e(P, Q) = g (a generator in GT):
#
# (a) e(3P, 5Q) = ?

assert (3*G1_lifted).weil_pairing(5*G2, r) == (GT ** 15)

# A: g^{15}

#
# (b) e(aP, bQ) * e(cP, dQ) = ?  (simplify to one pairing)
#
# ============================================================

a, b, c, d = 2, 3, 5, 7
assert (a*G1_lifted).weil_pairing(b*G2, r) * (c*G1_lifted).weil_pairing(d*G2, r) == (GT ** (a*b+c*d))

# A: g^{ab + cd}

# ============================================================
# Q2: Groth16 intuition (pen & paper)
#
# Prover knows a, b such that a * b = c (c is public).
# Prover sends A = aG1, B = bG2.
# Verifier knows C = e(G1, G2)^c.
#
# How does the verifier check a * b = c using a pairing?
#
# ============================================================

a, b, c = 2, 3, 6
A = a*G1_lifted
B = b*G2
# V performs pairing on A,B and C
lhs = A.weil_pairing(B, r)
rhs = G1_lifted.weil_pairing(G2, r) ** c
assert lhs == rhs, f"lhs rhs mismatch, {lhs=}, {rhs=}"


# A:
#   V knows A = aG1, B = bG2, c, then C = e(G1, G2)^c
#   V put everything on pairing
#   lhs = e(A, B) = e(a G1, b G2) = e(G1, G2)^{ab}
#   rhs = C                       = e(G1, G2)^c
#   V checks lhs ?= rhs


# ============================================================
# Q3: Why two different groups? (pen & paper)
#
# If G1 = G2 (same group), what security problem arises?
# Hint: related to the DDH assumption.
#
# ============================================================


# A:
#   recall DDH: given G, xG, yG, Z, we can decide Z ?= xyG.
#   if G1 = G2, we can solve DDH on G1 if we're given X=xG, Y=yG, Z=xyG by
#   lhs = e(X, Y), rhs = e(G, Z), checking lhs ?= rhs
#
# Q (from Kevin): does this mean, when we construct a pairing, we must find a GT s.t.
# it's pairing doesn't work with G1 x G1 \rightarrow GT?
# A (from Claude): no, all we have to do it to pick G1 \subset E(F_p), G2 \subset (F_{p^k}^*) s.t. G1, G2 don't have known
# homomorphism. Also, G1 != G2 and o(G1) = g(G2) = r
