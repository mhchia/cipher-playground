"""
SageMath Cheatsheet — common operations for lattice crypto exercises.

Run with: sage -python sage_cheatsheet.py

IMPORTANT: In .py files, use ** for exponentiation, NOT ^.
    ^ is XOR in Python and will raise RuntimeError with sage objects.
    In .sage files, ^ works because sage preprocesses it to **.
"""
from sage.all import *

# ============================================================
# 1. Finite Fields
# ============================================================

p = 101
Fp = GF(p)                          # prime field F_101
print(f"Fp = {Fp}, order = {Fp.order()}")

a = Fp(42)
print(f"42 inverse mod {p} = {a**-1}")  # multiplicative inverse
print(f"42^10 mod {p} = {a**10}")

# Extension field F_{p^2} with irreducible poly over F_p
Fp2 = GF(p**2, 'a')
a_gen = Fp2.gen()                    # generator of the extension
print(f"Fp2 = {Fp2}, generator = {a_gen}, minimal poly = {a_gen.minpoly()}")

# ============================================================
# 2. Integers mod q  (Z/qZ)
# ============================================================

q = 12289                            # common lattice prime (q ≡ 1 mod 512)
Zq = Integers(q)                     # Z/12289Z (ring, not field if q composite)
# or equivalently: Zq = Zmod(q)
print(f"\nZq = {Zq}")
print(f"3001 * 4096 mod q = {Zq(3001) * Zq(4096)}")

# ============================================================
# 3. Polynomial Rings — Z_q[X]
# ============================================================

R = PolynomialRing(Zq, 'X')
X = R.gen()

f = 3*X**3 + 5*X + 7
g = X**2 + 2*X + 1
print(f"\nf = {f}")
print(f"g = {g}")
print(f"f * g = {f * g}")
print(f"f mod g = {f % g}")
print(f"f(42) = {f(Zq(42))}")        # evaluate at a point

# Coefficients access
print(f"f.list() = {f.list()}")       # [7, 5, 0, 3] (constant term first)
print(f"f.degree() = {f.degree()}")

# ============================================================
# 4. Quotient Rings — Z_q[X] / (X^d + 1)  (cyclotomic ring R_q)
# ============================================================

d = 8
Rq = R.quotient(X**d + 1, 'x')      # R_q = Z_q[X] / (X^8 + 1)
x = Rq.gen()
print(f"\nRq = {Rq}")

# In quotient ring, X^d ≡ -1
h = x**7 + 3*x**2 + 1
k = x**6 + x + 5
product = h * k                      # automatically reduced mod X^8 + 1
print(f"h = {h}")
print(f"k = {k}")
print(f"h * k (mod X^8+1) = {product}")

# Lift back to polynomial ring to access coefficients
product_poly = R(product.lift())
print(f"coefficients = {product_poly.list()}")

# ============================================================
# 5. Finding Roots of Unity (needed for NTT)
# ============================================================

# For NTT over Z_q[X]/(X^d+1), we need a primitive 2d-th root of unity
# Requirement: q ≡ 1 (mod 2d)
two_d = 2 * d
assert (q - 1) % two_d == 0, f"q-1 must be divisible by 2d={two_d}"

# Find primitive 2d-th root of unity
g_unit = Zq.unit_group_exponent()     # order of multiplicative group = q-1
# A generator of the full multiplicative group:
gen_mult = Zq(primitive_root(q))
# 2d-th root of unity:
omega_2d = gen_mult ** ((q - 1) // two_d)
print(f"\nprimitive 2d-th root of unity: omega = {omega_2d}")
print(f"omega^(2d) = {omega_2d**two_d}")   # should be 1
print(f"omega^d = {omega_2d**d}")           # should be -1 (i.e., q-1)

# ============================================================
# 6. Matrices and Vectors over Z_q
# ============================================================

n, m = 4, 8
MS = MatrixSpace(Zq, n, m)
A = MS.random_element()
print(f"\nRandom {n}x{m} matrix A over Z_{q}:")
print(A)

# Vector spaces
VS = MatrixSpace(Zq, 1, m)           # row vectors
v = VS.random_element()
print(f"\nRandom vector v: {v}")

# Matrix-vector product (v as column)
v_col = vector(Zq, [Zq.random_element() for _ in range(m)])
result = A * v_col
print(f"A * v = {result}")

# ============================================================
# 7. Norms (for lattice crypto)
# ============================================================

coeffs = [Zq(c) for c in [3, -5, 12, 0, -1, 7, 2, -8]]
# Centered representation: map to {-⌊q/2⌋, ..., ⌊q/2⌋}
def centered_mod(x, q):
    """Convert from {0,...,q-1} to {-⌊q/2⌋,...,⌊q/2⌋}."""
    r = int(x) % q
    if r > q // 2:
        r -= q
    return r

centered = [centered_mod(c, q) for c in coeffs]
print(f"\nOriginal coeffs (mod q): {[int(c) for c in coeffs]}")
print(f"Centered repr: {centered}")

# ℓ∞ norm (max absolute coefficient)
linf = max(abs(c) for c in centered)
print(f"ℓ∞ norm = {linf}")

# ℓ2 norm
l2 = sqrt(sum(c**2 for c in centered))
print(f"ℓ2 norm = {float(l2):.4f}")

# ============================================================
# 8. Useful Sage Functions
# ============================================================

print(f"\nis_prime(12289) = {is_prime(12289)}")
print(f"factor(12288) = {factor(12288)}")
print(f"euler_phi(12289) = {euler_phi(12289)}")
print(f"primitive_root(12289) = {primitive_root(12289)}")
print(f"next_prime(100) = {next_prime(100)}")
print(f"random_prime(2^16) = {random_prime(2**16)}")
print(f"CRT([2, 3], [5, 7]) = {CRT([2, 3], [5, 7])}")  # Chinese Remainder Theorem
