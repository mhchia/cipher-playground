# Lattice Cryptography Exercises — Round 1

Parameters for all exercises: q = 17, d = 4 (ring: Z_17[X]/(X^4 + 1))

## Exercise 1: Ring Arithmetic in Z_q[X]/(X^d + 1)

### 1a. Polynomial multiplication in the quotient ring

Define a = 1 + 2X + 3X^2 + 4X^3 and b = 5 + 6X + 7X^2 + 8X^3.
Compute a * b in Z_17[X]/(X^4 + 1).
Verify that X^4 ≡ -1 takes effect (i.e., higher-degree terms wrap around with a sign flip).

### 1b. Norm computation with centered representation

Compute the ℓ∞ norm and ℓ2 norm of a.

Use centered representation: map coefficients from {0,...,16} to {-8,...,8}.
Mapping: if coeff > q//2, subtract q. So 16 → -1, 15 → -2, ..., 9 → -8.

- ℓ∞ norm = max of |centered coefficients|
- ℓ2 norm = sqrt(sum of squared centered coefficients)

(Note: canonical embedding norm is a different thing used in security proofs.
Protocol-level norm checks in LatticeFold/SALSAA operate on coefficient norms directly.)

### 1c. Norm growth under multiplication

Randomly generate pairs of "small norm" polynomials (coefficients ∈ {-1, 0, 1}).
Multiply them and observe the ℓ∞ norm of the product.
Repeat ~10 times. How large does the product norm get compared to the inputs?


## Exercise 2: NTT / INTT

Prerequisite: 17 ≡ 1 (mod 8), so a primitive 8th root of unity exists in Z_17.

### 2a. Find the primitive root

Find ω ∈ Z_17 such that ω^8 = 1 and ω^4 ≠ 1.
Hint: start from a multiplicative generator g of Z_17*, then ω = g^((17-1)/8).

### 2b. Implement NTT

Given a(X) with coefficients [a0, a1, a2, a3], compute:
  NTT(a) = [a(ω^1), a(ω^3), a(ω^5), a(ω^7)]

The evaluation points are ODD powers of ω, not even powers.
Question to think about: why odd powers?
Hint: the roots of X^4 + 1 in Z_17 are exactly ω^1, ω^3, ω^5, ω^7.

### 2c. Implement INTT

Implement the inverse NTT. Verify that INTT(NTT(a)) = a.

### 2d. Multiplication via NTT

Verify: NTT(a * b) = NTT(a) ⊙ NTT(b)  (pointwise multiplication)
Compare the result with Exercise 1a.


## Exercise 3: Ajtai Commitment

Parameters: κ = 2 (commitment vector length), m = 4 (witness vector length)
Each entry of the matrix and witness is a ring element in Z_17[X]/(X^4 + 1).

### 3a. Commit

1. Generate a random matrix A ∈ R_q^{κ × m}  (each entry is a random polynomial of degree < 4)
2. Generate a small-norm witness z ∈ R_q^m  (each polynomial has coefficients ∈ {-1, 0, 1})
3. Compute commitment c = A · z ∈ R_q^κ  (matrix-vector product over the ring)

### 3b. Binding experiment

Try to find z' ≠ z such that A · z' = c (i.e., A · (z - z') = 0).

- If you restrict z' to small norm: this should be hard (Module-SIS hardness).
- If you allow z' to have arbitrary norm: is it easier? Try it.
  Hint: if ring elements were invertible, you could solve a linear system.
  But not all ring elements are invertible in R_q — observe what happens.
