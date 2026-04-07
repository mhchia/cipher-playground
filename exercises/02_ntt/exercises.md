# Exercise 2: NTT / INTT

Parameters: q = 17, d = 4 (ring: Z_17[X]/(X^4 + 1))
Prerequisite: 17 ≡ 1 (mod 8), so a primitive 8th root of unity exists in Z_17.

## 2a. Find the primitive root

Find ω ∈ Z_17 such that ω^8 = 1 and ω^4 ≠ 1.
Hint: start from a multiplicative generator g of Z_17*, then ω = g^((17-1)/8).

## 2b. Implement NTT

Given a(X) with coefficients [a0, a1, a2, a3], compute:
  NTT(a) = [a(ω^1), a(ω^3), a(ω^5), a(ω^7)]

The evaluation points are ODD powers of ω, not even powers.
Question to think about: why odd powers?
Hint: the roots of X^4 + 1 in Z_17 are exactly ω^1, ω^3, ω^5, ω^7.

## 2c. Implement INTT

Implement the inverse NTT. Verify that INTT(NTT(a)) = a.

## 2d. Multiplication via NTT

Verify: NTT(a * b) = NTT(a) ⊙ NTT(b)  (pointwise multiplication)
Compare the result with Exercise 1a.
