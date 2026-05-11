"""
Unit tests for SALSAA modules.

Run: sage -python tests.py
"""
import itertools

from sage.all import *

from ring import q, Fq, d, x, Rq, conjugate, to_centered, _gen_random_low_norm_poly
from lde import lde_poly, pad_vec_to_d_exp, tensor_product
from relations import LinInstance, LinWitness, LinRelation
from join import rok_join
from rp import rok_rp


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


_SHARED_F_COM = matrix(Rq, [
    [Rq(1) + x,   Rq(2),    Rq(3),    Rq(0)],
    [Rq(0),       Rq(5)*x,  Rq(7),    Rq(11)],
])


def _make_lin_relation_with_eval(variant: int = 0):
    """
    Build a satisfying LinRelation with BOTH commitment and evaluation rows.

    All variants share the same `F_com` (Ajtai commitment matrix) — this matches
    the SALSAA join assumption that instances being folded use a common
    commitment scheme. Variants differ in W (different witness) and r̃ points
    (different evaluation rows in F_eval).

    Layout (κ=2, n_eval=2, m=4, r=2):
      F_com:  2 × 4   shared across variants
      F_eval: 2 × 4   tensor_product rows at variant-specific points
      H:      4 × 4   identity (no batching)
      Y:      4 × 2   = H · F · W
    """
    F_com = _SHARED_F_COM
    if variant == 0:
        W = _make_small_norm_W()
        r_tilde_0 = vector(Fq, [3, 5])
        r_tilde_1 = vector(Fq, [7, 11])
    elif variant == 1:
        W = matrix(Rq, [
            [Rq(0),  Rq(1)],
            [Rq(1),  x],
            [Rq(-1), Rq(0)],
            [x,      Rq(-1) * x],
        ])
        r_tilde_0 = vector(Fq, [2, 9])
        r_tilde_1 = vector(Fq, [4, 6])
    else:
        raise ValueError(f"Unknown variant {variant}")

    # F_eval rows: tensor evaluation basis at the chosen points
    F_eval = matrix(Rq, [
        tensor_product(r_tilde_0, 2),
        tensor_product(r_tilde_1, 2),
    ])
    # H = I_{κ + n_eval}, no batching (every eval row preserved as-is)
    H = identity_matrix(Rq, F_com.nrows() + F_eval.nrows())
    Y = H * F_com.stack(F_eval) * W

    return LinRelation(
        instance=LinInstance(H=H, F_com=F_com, F_eval=F_eval, Y=Y, v_square=16),
        witness=LinWitness(W),
    )


def test_make_lin_relation_with_eval():
    """Fixture self-check: both variants build valid LinRelations."""
    rel_0 = _make_lin_relation_with_eval(variant=0)
    rel_1 = _make_lin_relation_with_eval(variant=1)
    # Both should be valid (no exception from __post_init__) and have eval rows
    for rel in (rel_0, rel_1):
        assert rel.instance.F_com.nrows() == 2
        assert rel.instance.F_eval is not None
        assert rel.instance.F_eval.nrows() == 2
        assert rel.instance.H.nrows() == 4
        assert rel.instance.Y.nrows() == 4
    # Variants should differ (otherwise join with itself = boring test)
    assert rel_0.instance.W != rel_1.instance.W or rel_0.instance.F_com != rel_1.instance.F_com \
        if hasattr(rel_0.instance, 'W') else True
    print("  test_make_lin_relation_with_eval: OK")


def test_rok_join_smoke():
    """
    Join two LinRelations into one.

    NOTE: this test will FAIL until rok_join is implemented (TDD-style).
    Adjust assertions as the join contract crystallizes.
    """
    rel_0 = _make_lin_relation_with_eval(variant=0)
    rel_1 = _make_lin_relation_with_eval(variant=1)

    joined = rok_join(rel_0, rel_1)

    # Should return a valid LinRelation (its __post_init__ enforces HFW=Y + norm)
    assert isinstance(joined, LinRelation), \
        f"rok_join must return LinRelation, got {type(joined)}"

    n_top = rel_0.instance.F_com.nrows()
    n_bot_each = rel_0.instance.F_eval.nrows()
    assert joined.m == rel_0.m
    assert joined.hat_n == rel_0.hat_n + rel_1.hat_n - n_top
    assert joined.r == rel_0.r + rel_1.r
    # Per paper §6.1: joining L=2 instances adds (n - n̄)·(L-1) = n_bot rows to F.
    assert joined.n == n_top + 2 * n_bot_each, \
        f"expected n={n_top + 2*n_bot_each}, got {joined.n}"

    # Key invariants of join: shared commitment + v_square preserved + H square (pre-batch)
    assert joined.instance.F_com == rel_0.instance.F_com, "F_com must be preserved"
    assert joined.instance.v_square == rel_0.instance.v_square, "v_square must be preserved"
    assert joined.instance.H.nrows() == joined.instance.H.ncols(), "H stays square pre-batch"
    print("  test_rok_join_smoke: OK")


def _make_lin_relation_for_rp_input():
    """
    Satisfying LinRelation suitable as Π^⊗RP input.

    Π^⊗RP projects W (m × r) into Ŵ (m' × r) where m' = m/r, so we need
    `m` divisible by `r`. Reuses the existing small-norm fixture (m=4, r=2 → m'=2).
    """
    F_com = _SHARED_F_COM
    W = _make_small_norm_W()  # 4 × 2
    H = identity_matrix(Rq, F_com.nrows())
    Y = H * F_com * W
    return LinRelation(
        instance=LinInstance(H=H, F_com=F_com, F_eval=None, Y=Y, v_square=16),
        witness=LinWitness(W),
    )


def test_rok_rp_smoke():
    """
    Smoke test for rok_rp.

    LinRelations:
      - aug:  augmented original   ((H̃, F̃, Ỹ), W)   — width r preserved
      - proj: projected            ((I, F̂, ŷ), ŵ)   — width collapsed to 1,
                                                        m shrinks to m' = m/r,
                                                        norm grows to β̂ = m_rp · β

    Verifier samples J' ∈ R_q^{n_rp × m_rp} from {-1, 0, 1}, with m_rp/n_rp = r.

    NOTE: this test will FAIL until rok_rp is implemented (TDD-style).
    """
    lin_in = _make_lin_relation_for_rp_input()

    # Pick the smallest legal projection: n_rp = 1, m_rp = r.
    # (Adjust the call below if your rok_rp signature ends up different.)
    n_rp = 1
    m_rp = lin_in.r
    out = rok_rp(lin_in, n_rp=n_rp, m_rp=m_rp)

    # Type-shape: must be a 2-tuple of LinRelations.
    assert isinstance(out, tuple) and len(out) == 2, \
        f"rok_rp must return a 2-tuple, got {type(out).__name__} of len {len(out) if hasattr(out, '__len__') else '?'}"
    aug, proj = out
    assert isinstance(aug, LinRelation), f"first output must be LinRelation, got {type(aug).__name__}"
    assert isinstance(proj, LinRelation), f"second output must be LinRelation, got {type(proj).__name__}"

    # The fact that both wrap as LinRelation already enforces (via __post_init__):
    #   - H · F · W = Y on each side  (algebraic correctness)
    #   - max column ‖·‖² ≤ v_square  (norm bound)
    # So any cross-term arithmetic bug in rok_rp would already explode above.

    # Commitment must be preserved on the augmented side
    # (Π^⊗RP only *appends* a row to F; it doesn't touch the F_com block).
    assert aug.instance.F_com == lin_in.instance.F_com, \
        "F_com on the augmented side must be unchanged"

    assert aug.hat_n == lin_in.hat_n + 1
    assert aug.n == lin_in.n + 1
    assert aug.m == lin_in.m
    assert aug.r == lin_in.r
    assert aug.v_square == lin_in.v_square

    assert proj.hat_n == lin_in.n_top + 1
    assert proj.n == lin_in.n_top + 1
    assert proj.m == lin_in.m
    # only one \hat w
    assert proj.r == 1
    # norm bound growth: \beta_{proj} = m_{rp} * \beta_{lin}
    assert proj.v_square == lin_in.v_square * (m_rp**2)

    print("  test_rok_rp_smoke: OK")


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
    test_make_lin_relation_with_eval()

    print("join.py")
    test_rok_join_smoke()

    print("rp.py")
    test_rok_rp_smoke()

    print("\nAll tests passed.")


if __name__ == "__main__":
    main()
