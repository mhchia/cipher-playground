"""
Low-Degree Extension (LDE) over the hypercube [D]^l using a Lagrange basis.
"""
import itertools
from sage.all import *

from ring import Fq



def tensor_product(_r: vector, _D):
    """
    Calculate \tilde r so LDE[w](r) = <\tilde r, w>
    E.g. D=2, return (eq(r_1, 0), eq(r_1, 1)) ⊗ (eq(r_2, 0), eq(r_2, 1)) ... ⊗ (eq(r_l, 0), eq(r_l, 1))
    """
    l = len(_r)
    #
    # Prepare for eqs
    # [
    #   [eq(x_0, 0), ..., eq(x_0, d-1)],
    #   ...,
    #   [eq(x_{l-1}, 0), ..., eq(x_{l-1}, d-1)],
    # ]
    eqs = []
    for j in range(l):
        # in the end of this loop
        # xi_eq_k = (eq(x_j, 0), ..., eq(x_j, d-1))
        xi_eq_k = []
        for k in range(_D):
            lagrange = prod([
                (_r[j] - k_prime) / (Fq(k) - k_prime)
                for k_prime in range(_D) if k_prime != k
            ])
            xi_eq_k.append(lagrange)
        eqs.append(xi_eq_k)

    # Do tensor product on our own

    # eqs_flat = [eq(x, (0,..., 0)), eq(x, (d-1, ..., d-1)]
    eqs_flat = []
    for idx in itertools.product(range(_D), repeat=l):
        # e.g. idx = (0, ..., 0)
        # eq(x, (0,..., 0)) = eq(x_0, 0) * ... * eq(x_{l-1}, 0)
        eq = prod([eqs[j][idx[j]] for j in range(l)])
        eqs_flat.append(eq)
    return eqs_flat


def lde_poly(w: vector, _D: int):
    """
    Build LDE[w] as a symbolic multivariate polynomial over [D]^l using Lagrange basis.
    Returns (polynomial, xs) where xs are the symbolic variables.
    """
    # Pad w to size=d^l
    w_pad, l = pad_vec_to_d_exp(w, _D)

    P = PolynomialRing(w.base_ring(), [f"x{i}" for i in range(l)])
    xs = P.gens()

    eqs_flat = tensor_product(xs, _D)

    # <w, tensor(r)>
    # res = <w, eqs_flat>
    #     = w_0*eq(x, (0,..., 0)) + ... + w_{d^l-1}*eq(x, (d-1, ..., d-1)
    tilde_f = sum(w_pad[i] * eqs_flat[i] for i in range(_D**l))
    return tilde_f, xs


# Pad w to size d^l to make it on hypercube
def pad_vec_to_d_exp(w: vector, _D: int):
    # First padding w to 8 elements so we have [d]^3
    len_w = len(w)
    l = 0
    while _D**l < len_w:
        l += 1
    w_padded = vector(w.base_ring(), list(w) + [Fq(0)]*(_D**l - len_w))
    return w_padded, l


def test_lde_poly():
    # Test with multiple d values on the same w
    w = vector(Fq, [1, 2, 8, 10, 3, 5])

    for D in [2, 3, 4]:
        w_pad, l = pad_vec_to_d_exp(w, D)
        poly, xs = lde_poly(w, D)

        # Every point in [d]^l should recover the padded value
        for bit_repr in itertools.product(range(D), repeat=l):
            bit_repr_reversed = bit_repr[::-1]
            idx = sum([bit * D**i for i, bit in enumerate(bit_repr_reversed)])
            val = poly.subs({xs[i]: bit_repr[i] for i in range(l)})
            assert w_pad[idx] == val, \
                f"lde_poly mismatch at {bit_repr=}, {D=}: expected {w_pad[idx]}, got {val}"

    # Test: size already equals d^l (no padding needed)
    w_exact = vector(Fq, [1, 2, 3, 4])
    for D in [2, 4]:
        w_pad, l = pad_vec_to_d_exp(w_exact, D)
        poly, xs = lde_poly(w_exact, D)
        assert len(w_pad) == len(w_exact) or w_pad[len(w_exact):] == vector(Fq, [0]*(len(w_pad)-len(w_exact)))
        for bit_repr in itertools.product(range(D), repeat=l):
            bit_repr_reversed = bit_repr[::-1]
            idx = sum([bit * D**i for i, bit in enumerate(bit_repr_reversed)])
            assert w_pad[idx] == poly.subs({xs[i]: bit_repr[i] for i in range(l)})

    # Test: single element
    w_single = vector(Fq, [7])
    for D in [2, 3]:
        poly, xs = lde_poly(w_single, D)
        if len(xs) == 0:
            assert poly == Fq(7)
        else:
            assert poly.subs({xs[0]: 0}) == Fq(7)
