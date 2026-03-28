from sage.all import *


# ECC: y^2 = x^3 + ax + b
# Define E: y^2 = x^3 + 2 mod 37
p = 37
F = GF(p)
# E(F_{37})
E = EllipticCurve(F, [0, 2])  # a=0, b=2
print(f"{E=}")

# by Hasse's theorem, q=|E(F_{37}| is close to p

# 1. find all elements on the curve by enumerating x, y (since x, y \in F_{37})
elements: list[tuple[F, F]] = []
for x in map(F, range(p)):
    rhs = x**3 + 2
    if rhs.is_square():
        y_plus = rhs.sqrt()
        elements.append((x, y_plus))
        # Don't forget -y!
        elements.append((x, -y_plus))
# Don't forget infinite point O
order_E = len(elements) + 1
print(f"{order_E=}")
# get prime factors
factors = factor(order_E)
# get all possible orders from divisors (by Lagrange Theoerem)
# [1, 7, 49]
_divisors = divisors(order_E)
print(f"Possible orders of the subgroups: {_divisors}")
# {7}. Get the largest prime factor
r = factors[-1][0]
# n (order_E) = h * r  <- subgroup with prime order `r`
#               ^cofactor
h = order_E / r

# 2. Find a group element P in the subgroup with prime order r
# \forall g \in G, n*g = 0 -> (hr)g = 0 -> r(hg) = 0
#   -> o(hg) = r

def find_subgroup_gen() -> E:
    for ele in elements:
        G = E(ele[0], ele[1])
        P = h * G
        if P.order() == r:
            return P
    raise Exception(f"couldn't find a group element with order {r=}")


G = find_subgroup_gen()
print(f"Found point G {G} with order {r=}")


# 3. Given arbitrary point Q, check if it's in the prime order `r` subgroup <G>.
import random
G_generated = [i*G for i in range(r)]

def is_in_subgroup(Q: E):
    # since r is **prime** and <G> is cyclic, every element in G is a generator.
    # rQ = 0 <-> o(Q) = r <-> Q \in <G>
    return (r*Q).is_zero()

Q = E(*random.choice(elements))
assert (Q in G_generated) == is_in_subgroup(Q)

# if cofactor = 1, the order of ec group is prime and thus every group element
# G is a generator of the entire group. order(G) = n = r*1
