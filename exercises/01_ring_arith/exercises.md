# Exercise 1: Ring Arithmetic in Z_q[X]/(X^d + 1)

Parameters: q = 17, d = 4 (ring: Z_17[X]/(X^4 + 1))

## 1a. Polynomial multiplication in the quotient ring

Define a = 1 + 2X + 3X^2 + 4X^3 and b = 5 + 6X + 7X^2 + 8X^3.
Compute a * b in Z_17[X]/(X^4 + 1).
Verify that X^4 ≡ -1 takes effect (i.e., higher-degree terms wrap around with a sign flip).

## 1b. Norm computation with centered representation

Compute the ℓ∞ norm and ℓ2 norm of a.

Use centered representation: map coefficients from {0,...,16} to {-8,...,8}.
Mapping: if coeff > q//2, subtract q. So 16 → -1, 15 → -2, ..., 9 → -8.

- ℓ∞ norm = max of |centered coefficients|
- ℓ2 norm = sqrt(sum of squared centered coefficients)

(Note: canonical embedding norm is a different thing used in security proofs.
Protocol-level norm checks in LatticeFold/SALSAA operate on coefficient norms directly.)

## 1c. Norm growth under multiplication

Randomly generate pairs of "small norm" polynomials (coefficients ∈ {-1, 0, 1}).
Multiply them and observe the ℓ∞ norm of the product.
Repeat ~10 times. How large does the product norm get compared to the inputs?
