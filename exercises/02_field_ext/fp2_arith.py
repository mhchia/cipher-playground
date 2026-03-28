from sage.all import *

# ============================================================
# Fp² Arithmetic
# Extension field F_7^2 = F_7[u] / (u^2 + 1), so u^2 = -1
# ============================================================

p = 7
F = GF(p)
R = PolynomialRing(F, 'x')
x = R.gen()
irreducible_poly = x**2+1
# extended field F_{p^2} with k=2
F2 = F.extension(irreducible_poly, 'u')
print(f"{F2=}")
# generator of F_{p^2}
u = F2.gen()

print(f"F2 = {F2}")
print(f"u^2 = {u**2}")

# ============================================================
# Q1: Basic arithmetic in F_7^2
# Compute by hand first, then verify with Sage.
# ============================================================

a = 3 + 2*u
b = 5 + 4*u

# remember F_{p^2} = {a + bu | a, b \in F_p}
# (a) addition
assert a + b == 1 + 6*u
# (b) multiplication
assert a * b == u
# (c) inverse of a — verify a * a^(-1) == 1
assert a**(-1) == 4 + 2*u


# ============================================================
# Q2: Structure of F_7^2
# - How many elements in F_7^2?

# F_7 has 7 elements (prime, cyclic). F_{7^2} has 7^2 elements.
assert F2.order() == p**2


# - What is the order of (F_7^2)*?

# We find it by ourself. Should be the same as |F_{7^2}|-1 (0 excluded).
F2_eles = F2.list()

def find_mul_subgroup():
    max_subgroup = [F2.one()]
    max_subgroup_size = 1
    for ele in F2_eles:
        # skip 0
        if ele.is_zero():
            continue
        cur = ele
        # ele_gen_subgroup = <ele>
        ele_gen_subgroup = [cur]
        # stop when it reaches 1 (full cycle)
        while cur != F2.one():
            cur *= ele
            ele_gen_subgroup.append(cur)
        cur_subgroup_size = len(ele_gen_subgroup)
        if cur_subgroup_size > max_subgroup_size:
            max_subgroup_size = cur_subgroup_size
            max_subgroup = ele_gen_subgroup
    return max_subgroup

F2_mul_eles = find_mul_subgroup()
F2_mul_order = len(F2_mul_eles)
print(f"{F2_mul_order=}")
## Check with sage's F.multiplicative_generator() and multiplicative_order()
F2_mul_sage = F2.multiplicative_generator()
assert F2_mul_order == F2_mul_sage.multiplicative_order(), f"mismatch: {F2_mul_order=}, {F2_mul_sage.order()=}"


# - Factor that order.
factors = factor(F2_mul_order)
print(f"factored multiplicative order: {factors=}")
# - Find a generator of (F_7^2)*.
F2_mul_gen = F2_mul_eles[0]
print(f"F_7^2 multiplicative generator: {F2_mul_gen=}")

# ============================================================




# ============================================================
# Q3: ECC over F_7^2
# E: y^2 = x^3 + 1
# ============================================================

# Build E over F_7 and F_7^2
# E: y^2 = x^3 + ax + b, a=0, b=1
E = EllipticCurve(F, [0, 1])
E2 = EllipticCurve(F2, [0, 1])


# Q: Compare |E(F_7)| vs |E(F_7^2)|
print(f"|E(F_7)|={E.order()}, |E(F_7^2)|={E2.order()}")

# Q: Are all points of E(F_7) also in E(F_7^2)? Why?
# A: It's trivial, every solution in E(F_7) holds in E(F_7^2)


# ============================================================
# Q4: Tower extensions (no code needed, just think)
#
# BN254 builds Fp12 via tower: Fp -> Fp2 -> Fp6 -> Fp12
# Extension degrees: 2, 3, 2
#
# - Total degree = ?
# A: 12.

# - Naive Fp2 multiplication: how many Fp muls?
# A: 4.

# - With Karatsuba?
# A: 3.


# - Rough cost comparison: tower vs direct Fp12?
# A:
#   - direct Fp12: 12 * 12 = 144
#       - since F_{p^12} = {a_0 +...+ a_11 w^11 | a_i \in F_p}. a * b -> 12 * 12 F_p muls
#   - tower: 3 * 9 * 3 = 81
#       - F_{p^12} = {a_0 + a_1 u | a_0, a_1 \in F_{p^6}}. a * b -> 2*2=4 F_{p^6} muls --Karatsuba--> 3 F_{p^6} muls
#       - F_{p^6} = {a_0 + a_1 u + a_2 u^2 | a_i \in F_{p^2}}. a * b -> 3 * 3 = 9 F_{p^2} muls
#       - F_{p^2} = {a_0 + a_1 u | a_i \in F_p}. a * b -> 2*2=4 --Karatsuba--> 3 F_p muls
# ============================================================
