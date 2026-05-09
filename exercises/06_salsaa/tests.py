"""
Unit tests for SALSAA modules.

Run: sage -python tests.py
"""
import itertools

from sage.all import *

from ring import q, Fq, d, x, Rq, conjugate, to_centered, _gen_random_low_norm_poly
from lde import lde_poly, pad_vec_to_d_exp, tensor_product
from relations import LinInstance, LinWitness, LinRelation


# ============================================================
# ring.py
# ============================================================

def test_to_centered():
    assert to_centered(Fq(0)) == 0
    assert to_centered(Fq(1)) == 1
    assert to_centered(Fq(q // 2)) == q // 2          # boundary stays positive
    assert to_centered(Fq(q // 2 + 1)) == -(q // 2)   # next step flips sign
    assert to_centered(Fq(q - 1)) == -1
    print("  test_to_centered: OK")


def test_conjugate():
    # σ(constant) = constant
    assert conjugate(Rq(5)) == Rq(5)
    assert conjugate(Rq(0)) == Rq(0)

    # σ(x) = -x^{d-1}
    assert conjugate(x) == -x**(d - 1)

    # σ²(a) = a (involution)
    a = x**2 + 3*x + 7
    assert conjugate(conjugate(a)) == a

    # σ(a · b) = σ(a) · σ(b)  (ring homomorphism)
    a = x + 1
    b = x**2 + 2
    assert conjugate(a * b) == conjugate(a) * conjugate(b)
    print("  test_conjugate: OK")


# ============================================================
# lde.py
# ============================================================

def test_lde_poly():
    """LDE recovers values on hypercube."""
    w = vector(Fq, [1, 2, 8, 10, 3, 5])

    for D in [2, 3, 4]:
        w_pad, l = pad_vec_to_d_exp(w, D)
        poly, xs = lde_poly(w, D)

        for bit_repr in itertools.product(range(D), repeat=l):
            bit_repr_reversed = bit_repr[::-1]
            idx = sum([bit * D**i for i, bit in enumerate(bit_repr_reversed)])
            val = poly.subs({xs[i]: bit_repr[i] for i in range(l)})
            assert w_pad[idx] == val, \
                f"lde_poly mismatch at {bit_repr=}, {D=}: expected {w_pad[idx]}, got {val}"

    w_exact = vector(Fq, [1, 2, 3, 4])
    for D in [2, 4]:
        w_pad, l = pad_vec_to_d_exp(w_exact, D)
        poly, xs = lde_poly(w_exact, D)
        for bit_repr in itertools.product(range(D), repeat=l):
            bit_repr_reversed = bit_repr[::-1]
            idx = sum([bit * D**i for i, bit in enumerate(bit_repr_reversed)])
            assert w_pad[idx] == poly.subs({xs[i]: bit_repr[i] for i in range(l)})

    w_single = vector(Fq, [7])
    for D in [2, 3]:
        poly, xs = lde_poly(w_single, D)
        if len(xs) == 0:
            assert poly == Fq(7)
        else:
            assert poly.subs({xs[0]: 0}) == Fq(7)
    print("  test_lde_poly: OK")


def test_tensor_product_binary():
    """For binary r ∈ {0,1}^l, tensor_product(r) is the delta basis."""
    D = 2
    for b_0, b_1 in itertools.product([0, 1], repeat=2):
        r = vector(Fq, [b_0, b_1])
        result = tensor_product(r, D)
        # itertools.product([0,1], repeat=2) yields (0,0), (0,1), (1,0), (1,1)
        expected_idx = b_0 * 2 + b_1
        for i, val in enumerate(result):
            if i == expected_idx:
                assert val == 1, f"r={(b_0,b_1)} pos {i}: expected 1, got {val}"
            else:
                assert val == 0, f"r={(b_0,b_1)} pos {i}: expected 0, got {val}"
    print("  test_tensor_product_binary: OK")


def test_tensor_product_lde_consistency():
    """<w, tensor_product(r)> == lde_poly(w) evaluated at r — for arbitrary r."""
    D = 2
    w = vector(Fq, [1, 5, 9, 13])     # length 4 = D^l, l=2
    r = vector(Fq, [3, 7])            # arbitrary

    # Method 1: dot product with tensor
    tensor = tensor_product(r, D)
    val_tensor = sum(w[i] * tensor[i] for i in range(len(w)))

    # Method 2: lde_poly evaluation
    poly, xs = lde_poly(w, D)
    val_poly = poly.subs({xs[0]: r[0], xs[1]: r[1]})

    assert val_tensor == val_poly, \
        f"tensor_product / lde_poly mismatch: tensor={val_tensor}, poly={val_poly}"
    print("  test_tensor_product_lde_consistency: OK")


# ============================================================
# relations.py
# ============================================================

def _make_small_norm_W():
    """Build a deterministic small-norm W for tests."""
    return matrix(Rq, [
        [Rq(1),       Rq(0)],
        [Rq(0),       Rq(1)],
        [x,           Rq(-1)],
        [Rq(-1) * x,  Rq(1)],
    ])


def test_lin_relation_happy():
    H = identity_matrix(Rq, 2)
    F_com = matrix(Rq, [
        [Rq(1) + x,    Rq(2),         Rq(3),    Rq(0)],
        [Rq(0),        Rq(5)*x,       Rq(7),    Rq(11)],
    ])
    W = _make_small_norm_W()
    Y = H * F_com * W

    rel = LinRelation(
        instance=LinInstance(H=H, F_com=F_com, F_eval=None, Y=Y, v_square=16),
        witness=LinWitness(W),
    )
    assert rel.hat_n == 2
    assert rel.m == 4
    print("  test_lin_relation_happy: OK")


def test_lin_relation_wrong_y_raises():
    H = identity_matrix(Rq, 2)
    F_com = matrix(Rq, [
        [Rq(1), Rq(2), Rq(3), Rq(4)],
        [Rq(5), Rq(6), Rq(7), Rq(8)],
    ])
    W = _make_small_norm_W()
    Y_correct = H * F_com * W
    Y_wrong = Y_correct + matrix(Rq, [[Rq(1), Rq(0)], [Rq(0), Rq(0)]])

    raised = False
    try:
        LinRelation(
            instance=LinInstance(H=H, F_com=F_com, F_eval=None, Y=Y_wrong, v_square=16),
            witness=LinWitness(W),
        )
    except AssertionError:
        raised = True
    assert raised, "LinRelation should have raised when H*F*W != Y"
    print("  test_lin_relation_wrong_y_raises: OK")


def test_lin_relation_norm_too_big_raises():
    H = identity_matrix(Rq, 2)
    F_com = matrix(Rq, [
        [Rq(1), Rq(0), Rq(0), Rq(0)],
        [Rq(0), Rq(1), Rq(0), Rq(0)],
    ])
    # First entry has coeff 5 (centered = 5), so ‖·‖² = 25 — exceeds v_square=16
    W_big = matrix(Rq, [
        [Rq(5) * x,  Rq(0)],
        [Rq(0),      Rq(0)],
        [Rq(0),      Rq(0)],
        [Rq(0),      Rq(0)],
    ])
    Y = H * F_com * W_big

    raised = False
    try:
        LinRelation(
            instance=LinInstance(H=H, F_com=F_com, F_eval=None, Y=Y, v_square=16),
            witness=LinWitness(W_big),
        )
    except AssertionError:
        raised = True
    assert raised, "LinRelation should have raised when col_norm² > v_square"
    print("  test_lin_relation_norm_too_big_raises: OK")


def test_lin_instance_mismatched_widths_raises():
    """F_com.ncols() != F_eval.ncols() should fail dimension check."""
    H = identity_matrix(Rq, 2)
    F_com = matrix(Rq, [[Rq(1), Rq(2), Rq(3), Rq(4)],
                        [Rq(5), Rq(6), Rq(7), Rq(8)]])
    F_eval = matrix(Rq, [[Rq(1), Rq(2)]])     # 1 × 2  (wrong width!)
    Y = matrix(Rq, [[Rq(0), Rq(0)], [Rq(0), Rq(0)]])

    raised = False
    try:
        LinInstance(H=H, F_com=F_com, F_eval=F_eval, Y=Y, v_square=16)
    except AssertionError:
        raised = True
    assert raised, "LinInstance should have raised for mismatched F_com / F_eval widths"
    print("  test_lin_instance_mismatched_widths_raises: OK")


def test_lin_relation_dimensions_asymmetric():
    """Use n̂ ≠ r and m ≠ r to expose any row/col confusion in dim properties."""
    H = identity_matrix(Rq, 3)                                          # n̂ = 3
    F_com = matrix(Rq, 3, 4, lambda i, j: Rq(i * 4 + j + 1))            # 3 × 4 (n̄=3, m=4)
    W = matrix(Rq, [
        [Rq(1), Rq(0)],
        [Rq(0), Rq(1)],
        [Rq(1), Rq(1)],
        [Rq(0), Rq(0)],
    ])                                                                   # m=4, r=2
    Y = H * F_com * W                                                    # 3 × 2

    rel = LinRelation(
        instance=LinInstance(H=H, F_com=F_com, F_eval=None, Y=Y, v_square=16),
        witness=LinWitness(W),
    )
    assert rel.hat_n == 3
    assert rel.n == 3
    assert rel.m == 4
    assert rel.r == 2     # catches Y.nrows() vs Y.ncols() confusion
    print("  test_lin_relation_dimensions_asymmetric: OK")


def test_with_extra_eval_dimensions():
    """with_extra_eval grows F_eval / Y / H by the right amount."""
    H = identity_matrix(Rq, 2)
    F_com = matrix(Rq, [
        [Rq(1) + x,    Rq(2),         Rq(3),    Rq(0)],
        [Rq(0),        Rq(5)*x,       Rq(7),    Rq(11)],
    ])
    W = _make_small_norm_W()
    Y = H * F_com * W
    inst = LinInstance(H=H, F_com=F_com, F_eval=None, Y=Y, v_square=16)

    # 2 new eval rows; new_Y_rows must be consistent with H_new being identity-extended
    new_F_rows = matrix(Rq, [
        [Rq(1), Rq(0), Rq(0), Rq(0)],
        [Rq(0), Rq(1), Rq(0), Rq(0)],
    ])
    new_Y_rows = new_F_rows * W

    new_inst = inst.with_extra_eval(new_F_rows=new_F_rows, new_Y_rows=new_Y_rows)

    # F_eval grew from None → 2 rows
    assert new_inst.F_eval.nrows() == 2
    assert new_inst.F_eval.ncols() == 4

    # Y grew by 2 rows
    assert new_inst.Y.nrows() == inst.Y.nrows() + 2

    # H grew by 2 in BOTH dims (square grows together)
    assert new_inst.H.nrows() == inst.H.nrows() + 2
    assert new_inst.H.ncols() == inst.H.ncols() + 2

    # F_com unchanged
    assert new_inst.F_com == inst.F_com

    # Should also be wrap-able into a valid LinRelation (HFW=Y still holds)
    new_rel = LinRelation(instance=new_inst, witness=LinWitness(W))
    assert new_rel.hat_n == 4
    print("  test_with_extra_eval_dimensions: OK")


# ============================================================
# Runner
# ============================================================

def main():
    print("ring.py")
    test_to_centered()
    test_conjugate()

    print("lde.py")
    test_lde_poly()
    test_tensor_product_binary()
    test_tensor_product_lde_consistency()

    print("relations.py")
    test_lin_relation_happy()
    test_lin_relation_wrong_y_raises()
    test_lin_relation_norm_too_big_raises()
    test_lin_instance_mismatched_widths_raises()
    test_lin_relation_dimensions_asymmetric()
    test_with_extra_eval_dimensions()

    print("\nAll tests passed.")


if __name__ == "__main__":
    main()
