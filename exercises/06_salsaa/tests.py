"""
Unit tests for SALSAA modules.

Run: sage -python tests.py
"""
import itertools

from sage.all import *

from ring import q, Fq, d, x, Rq, conjugate, to_centered, _gen_random_low_norm_poly
from lde import lde_poly, pad_vec_to_d_exp, tensor_product
from relations import LinInstance, LinWitness, LinRelation
from rok import (
    rok_join,
    rok_rp,
    rok_fold,
    rok_batch,
    rok_decompose,
    get_l,
    balanced_b_ary_decompose_Fq,
    compose_Fq,
    decompose_W,
)
from salsaa import fold as salsaa_fold, gen_random_W, gen_random_F, gen_H
from ring import n_hat as ring_n_hat, n as ring_n, m as ring_m, r as ring_r, beta as ring_beta


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


# ============================================================
# decompose.py
# ============================================================

def test_get_l():
    """ℓ = number of base-b digits needed for the range [-β, β]."""
    # β=1, b=2: 2β+1=3, need 2 binary digits
    assert get_l(1, 2) == 2
    # β=4, b=2: 2β+1=9, need 4 binary digits
    assert get_l(4, 2) == 4
    # β=7, b=3: 2β+1=15, need 3 ternary digits
    assert get_l(7, 3) == 3

    # β=0 must fire assertion (zero range is degenerate)
    raised = False
    try:
        get_l(0, 2)
    except AssertionError:
        raised = True
    assert raised, "get_l should reject β=0"
    print("  test_get_l: OK")


def test_decompose_Fq_explicit():
    """Concrete digit lists for documented cases."""
    # 7 = 1 + 2 + 4 = 0b0111 → [1, 1, 1, 0]
    assert balanced_b_ary_decompose_Fq(Fq(7), 2, 4) == [Fq(1), Fq(1), Fq(1), Fq(0)]
    # Sign carries through: -7 → [-1, -1, -1, 0]
    assert balanced_b_ary_decompose_Fq(Fq(-7), 2, 4) == [Fq(-1), Fq(-1), Fq(-1), Fq(0)]
    # Zero → all zeros
    assert balanced_b_ary_decompose_Fq(Fq(0), 2, 4) == [Fq(0)] * 4
    # Balanced ternary stress test: 5 with b=3 exercises the carry step.
    # Non-balanced would give [2, 1] (digit 2 ∉ {-1, 0, 1}); balanced uses
    # carry to push 2 → -1 with +3 added to the next position:
    #   5 = (-1)·1 + (-1)·3 + 1·9 → [-1, -1, 1]
    assert balanced_b_ary_decompose_Fq(Fq(5), 3, 3) == [Fq(-1), Fq(-1), Fq(1)]
    # Sign symmetry for the same case
    assert balanced_b_ary_decompose_Fq(Fq(-5), 3, 3) == [Fq(1), Fq(1), Fq(-1)]
    print("  test_decompose_Fq_explicit: OK")


def test_compose_Fq_explicit():
    """Reverse of balanced_b_ary_decompose_Fq on the same documented cases."""
    assert compose_Fq([Fq(1), Fq(1), Fq(1), Fq(0)], 2) == Fq(7)
    assert compose_Fq([Fq(-1), Fq(-1), Fq(-1), Fq(0)], 2) == Fq(-7)
    assert compose_Fq([Fq(0)] * 4, 2) == Fq(0)
    # Balanced ternary recompose: [-1, -1, 1] · (1, 3, 9) = -1 - 3 + 9 = 5
    assert compose_Fq([Fq(-1), Fq(-1), Fq(1)], 3) == Fq(5)
    assert compose_Fq([Fq(1), Fq(1), Fq(-1)], 3) == Fq(-5)
    print("  test_compose_Fq_explicit: OK")


def test_decompose_Fq_roundtrip():
    """compose(decompose(f)) == f for all f ∈ [-β, β] and several (b, β)."""
    for b in [2, 3]:
        for beta in [1, 4, 7]:
            l = get_l(beta, b)
            for f_int in range(-beta, beta + 1):
                f = Fq(f_int)
                coeffs = balanced_b_ary_decompose_Fq(f, b, l)
                assert len(coeffs) == l, \
                    f"balanced_b_ary_decompose_Fq must return ℓ={l} digits, got {len(coeffs)}"
                f_back = compose_Fq(coeffs, b)
                assert f_back == f, \
                    f"roundtrip mismatch: f={f_int}, b={b}, β={beta}, " \
                    f"coeffs={coeffs}, got {f_back}"
    print("  test_decompose_Fq_roundtrip: OK")


def test_decompose_W_roundtrip():
    """W = Σ_k b^k · V_k where V = decompose_W(W, b, ℓ)."""
    W = matrix(Rq, [
        [Rq(1) + 2*x,    Rq(3)],
        [Rq(0),          Rq(-1) - x**2],
    ])
    b = 2
    beta = 4  # max coeff magnitude in W is 3; β=4 gives a safe ℓ
    l = get_l(beta, b)
    V = decompose_W(W, b, l)

    # Shape
    assert len(V) == l, f"expected ℓ={l} matrices, got {len(V)}"
    for k, V_k in enumerate(V):
        assert V_k.nrows() == W.nrows(), f"V_{k} row count mismatch"
        assert V_k.ncols() == W.ncols(), f"V_{k} col count mismatch"

    # Round-trip: Σ_k b^k · V_k must reassemble to W
    W_back = sum([(b ** k) * V[k] for k in range(l)])
    assert W_back == W, f"roundtrip failed:\n  W={W}\n  W_back={W_back}"
    print("  test_decompose_W_roundtrip: OK")


def test_decompose_W_norm_bound():
    """Every coefficient of every V_k lives in {-(b-1), ..., b-1}."""
    W = matrix(Rq, [
        [Rq(7),       Rq(-3) * x],
        [5 * x**2,    Rq(0)],
    ])
    b = 2
    beta = 7
    l = get_l(beta, b)
    V = decompose_W(W, b, l)

    for k, V_k in enumerate(V):
        for poly in V_k.list():
            for c in poly.list():
                cv = to_centered(c)
                assert abs(cv) <= b - 1, \
                    f"V_{k} has coeff with |cv|={abs(cv)} > b-1={b-1}"
    print("  test_decompose_W_norm_bound: OK")


def _make_lins_for_salsaa_fold(L: int = 2):
    """
    Build L satisfying LinRelations to feed into the top-level `salsaa.fold(lins)`.

    Mirrors the pattern from salsaa.main() but generates L instances instead of 1:
      H = I_{n_hat}                    (identity — no batching)
      F_com shared across all L        (per SALSAA join's commitment-scheme invariant)
      W per-instance, random low-norm  (each from `_gen_random_low_norm_poly`)
      F_eval = None                    (matches main(): commitment-only, no eval rows yet)
      Y = H · F_com · W                (computed per instance)

    Norm budget (matches main()):
      beta_square = d                  (each ring coef ∈ {-1, 0, 1} ⇒ ‖w_ij‖² ≤ d)
      v_square    = m · beta_square    (m ring entries per column)
    """
    H = gen_H(ring_n_hat, ring_n)
    F_com = gen_random_F(ring_n, ring_m)   # shared across all L
    beta_square = d
    v_square = ring_m * beta_square

    lins = []
    for _ in range(L):
        W = gen_random_W(ring_m, ring_r)
        Y = H * F_com * W
        lins.append(LinRelation(
            instance=LinInstance(H=H, F_com=F_com, F_eval=None, Y=Y, v_square=v_square),
            witness=LinWitness(W),
        ))
    return lins


def test_salsaa_fold_smoke():
    """
    Smoke test for `salsaa.fold(lins)` — the top-level round driver.

    `fold(lins)` is expected to chain together one round of the SALSAA pipeline:
      join → norm → ⊗RP → fold → join(2) → batch → b-decomp
    feeding L fresh LinRelations in and producing one accumulator LinRelation out.

    NOTE: TDD-style — will FAIL until salsaa.fold is implemented.
    Assertions are deliberately loose (just structural invariants); tighten
    once the chain shape and parameter choices crystallize.
    """
    L = 2
    lins = _make_lins_for_salsaa_fold(L)

    # Sanity-check the fixture itself before feeding into fold
    assert len(lins) == L
    for lin in lins:
        assert isinstance(lin, LinRelation)
    # All instances must share F_com (join's invariant)
    for lin in lins[1:]:
        assert lin.instance.F_com == lins[0].instance.F_com, \
            "All input LinRelations must share F_com"

    # Run one round
    out = salsaa_fold(lins)

    # Output type
    assert isinstance(out, LinRelation), \
        f"salsaa.fold must return LinRelation, got {type(out).__name__}"

    # m (witness rows / commitment width) is preserved across the chain
    assert out.m == lins[0].m, \
        f"m changed: {lins[0].m} → {out.m}"

    # F_com (Ajtai commitment matrix) is preserved across the chain
    # — every sub-protocol takes a row-side or column-side action but never
    # rewrites the commitment matrix itself.
    assert out.instance.F_com == lins[0].instance.F_com, \
        "F_com must be preserved end-to-end across the chain"

    # The LinRelation __post_init__ has already validated:
    #   - H · F · W = Y  on the output (chain produces a valid relation)
    #   - col_norm² ≤ v_square  (norm bounds were tracked correctly through the chain)
    # If any sub-protocol's bound is too tight or any algebra is wrong, the
    # offending step's __post_init__ explodes here with a precise message.

    print("  test_salsaa_fold_smoke: OK")


def test_rok_decompose_smoke():
    """
    Smoke test for rok_decompose (Π^b-decomp).

    Decomposes the witness into ℓ small-norm chunks:
        W = Σ_{k=0..ℓ-1} b^k · V_k         with V_k entries in {-⌊b/2⌋, ..., ⌊b/2⌋}
        ℓ = ⌈log_b(2β + 1)⌉

    Effect on shape:
        Ŵ = (V_0 | V_1 | ... | V_{ℓ-1})    r grows: r → ℓ·r
        Ŷ = (Z_0 | Z_1 | ... | Z_{ℓ-1})    where Z_k = H·F·V_k
        H, F, m, n, hat_n preserved
        v_square shrinks (per-entry coefficients tighten)

    NOTE: TDD-style — will FAIL until rok_decompose is implemented.
    Assumes signature `rok_decompose(lin, b)` mirroring the param pattern
    of `rok_fold(lin, r_out)` and `rok_batch(lin, n_target_eval_rows)`.
    """
    lin_in = _make_lin_relation_with_eval(variant=0)
    b = 2
    out = rok_decompose(lin_in, b)

    # Type
    assert isinstance(out, LinRelation), \
        f"rok_decompose must return LinRelation, got {type(out).__name__}"

    # m, n, hat_n unchanged (decomp only touches W column count and Y)
    assert out.m == lin_in.m, f"m changed: {lin_in.m} → {out.m}"
    assert out.n == lin_in.n, f"n changed: {lin_in.n} → {out.n}"
    assert out.hat_n == lin_in.hat_n, f"hat_n changed: {lin_in.hat_n} → {out.hat_n}"

    # F (both halves) and H unchanged — decomp doesn't touch the statement side
    assert out.instance.F_com == lin_in.instance.F_com, "F_com must be preserved"
    if lin_in.instance.F_eval is None:
        assert out.instance.F_eval is None
    else:
        assert out.instance.F_eval == lin_in.instance.F_eval, "F_eval must be preserved"
    assert out.instance.H == lin_in.instance.H, "H must be preserved"

    # Witness width grew by an integer factor (= ℓ > 1 for non-trivial β/b).
    # For the test fixture (v_square = 16 ⇒ β = 4) and b = 2, ℓ = ⌈log_2(9)⌉ = 4.
    assert out.r > lin_in.r, \
        f"r must strictly grow (decomp must widen witness for β > ⌊b/2⌋): {lin_in.r} → {out.r}"
    assert out.r % lin_in.r == 0, \
        f"r must grow by integer factor ℓ: {lin_in.r} → {out.r}"
    ell = out.r // lin_in.r
    assert ell > 1, f"ℓ must be > 1 for non-trivial decomposition, got {ell}"

    # v_square should not increase (per-entry bound tightens with decomposition).
    # The LinRelation __post_init__ has already validated that the *actual* column
    # norms of Ŵ fit within out.v_square, so this is just a sanity bound.
    assert out.v_square <= lin_in.v_square, \
        f"v_square must not grow: {lin_in.v_square} → {out.v_square}"

    # The LinRelation __post_init__ has already validated:
    #   - H · F · Ŵ = Ŷ                       (decomp arithmetic consistent)
    #   - max_k ‖V_k‖² ≤ out.v_square           (each level fits the tighter bound)

    print("  test_rok_decompose_smoke: OK")


def test_rok_batch_smoke():
    """
    Smoke test for rok_batch (Π^batch).

    Collapses the evaluation rows of H from `n̂ - n̄` down to `n_target_eval_rows`
    by random linear combination (Vandermonde-like in c), applied symmetrically
    to Y. The top `n̄` (commitment) rows of H are preserved as-is. F (both F_com
    and F_eval), W, m, r, v_square all preserved — batch only acts on H and Y.

    NOTE: TDD-style — will FAIL until rok_batch is implemented.
    """
    lin_in = _make_lin_relation_with_eval(variant=0)
    n_target_eval_rows = 1
    out = rok_batch(lin_in, n_target_eval_rows)

    # Type
    assert isinstance(out, LinRelation), \
        f"rok_batch must return LinRelation, got {type(out).__name__}"

    # hat_n must equal n̄ + n_target_eval_rows (commitment rows kept + batched eval rows)
    expected_hat_n = lin_in.n_top + n_target_eval_rows
    assert out.hat_n == expected_hat_n, \
        f"hat_n: expected {expected_hat_n} (= n_top {lin_in.n_top} + target {n_target_eval_rows}), got {out.hat_n}"

    # n, m, r unchanged (batch doesn't touch witness side or F)
    assert out.n == lin_in.n, f"n changed: {lin_in.n} → {out.n}"
    assert out.m == lin_in.m, f"m changed: {lin_in.m} → {out.m}"
    assert out.r == lin_in.r, f"r changed: {lin_in.r} → {out.r}"

    # F (both halves) unchanged
    assert out.instance.F_com == lin_in.instance.F_com, \
        "F_com must be preserved"
    if lin_in.instance.F_eval is None:
        assert out.instance.F_eval is None
    else:
        assert out.instance.F_eval == lin_in.instance.F_eval, \
            "F_eval must be preserved"

    # W unchanged — batch restates the same witness with fewer rows in H/Y
    assert out.witness.W == lin_in.witness.W, "W must be preserved"

    # v_square preserved (W is the same)
    assert out.v_square == lin_in.v_square, \
        f"v_square should be preserved (W unchanged): {lin_in.v_square} → {out.v_square}"

    # The LinRelation __post_init__ has already validated:
    #   - H' · F · W = Y'   (so the random combination was applied consistently)
    #   - column norm² ≤ v_square (still holds since W unchanged)
    # i.e., if rok_batch combines H/Y inconsistently, this test explodes in
    # __post_init__ with "relation doesn't hold".

    print("  test_rok_batch_smoke: OK")


def test_rok_fold_smoke():
    """
    Smoke test for rok_fold (Π^fold, paper §3.4.1).

    Reduces witness width r → 1 by random column linear combination:
        Ŵ := W · c    Ŷ := Y · c    with random c ∈ R_q^r
    H, F preserved. Output ((H, F, Ŷ), Ŵ) ∈ Σ^lin.

    NOTE: TDD-style — will FAIL until rok_fold is implemented.
    """
    lin_in = _make_lin_relation_with_eval(variant=0)
    r_out = 1
    out = rok_fold(lin_in, r_out)

    # Type + structural invariants
    assert isinstance(out, LinRelation), \
        f"rok_fold must return LinRelation, got {type(out).__name__}"

    # Witness width collapsed to r_out
    assert out.r == r_out, f"folded r must be {r_out}, got {out.r}"

    # Other dims unchanged
    assert out.hat_n == lin_in.hat_n, f"hat_n changed: {lin_in.hat_n} → {out.hat_n}"
    assert out.n == lin_in.n,         f"n changed: {lin_in.n} → {out.n}"
    assert out.m == lin_in.m,         f"m changed: {lin_in.m} → {out.m}"

    # F unchanged (fold only touches W and Y by linear combination)
    assert out.instance.F_com == lin_in.instance.F_com, \
        "F_com must be preserved"
    if lin_in.instance.F_eval is None:
        assert out.instance.F_eval is None
    else:
        assert out.instance.F_eval == lin_in.instance.F_eval, \
            "F_eval must be preserved"

    # H unchanged
    assert out.instance.H == lin_in.instance.H, "H must be preserved"

    # The LinRelation __post_init__ has already validated:
    #   - H · F · Ŵ = Ŷ
    #   - column norm² ≤ v_square (so v_square is at least sufficient)
    # i.e., if rok_fold sends a wrong c-application, this test would explode in
    # __post_init__ with a precise "relation doesn't hold" message.

    print("  test_rok_fold_smoke: OK")


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

    print("fold.py")
    test_rok_fold_smoke()

    print("batch.py")
    test_rok_batch_smoke()

    print("decompose.py")
    test_get_l()
    test_decompose_Fq_explicit()
    test_compose_Fq_explicit()
    test_decompose_Fq_roundtrip()
    test_decompose_W_roundtrip()
    test_decompose_W_norm_bound()
    test_rok_decompose_smoke()

    print("salsaa.py")
    test_salsaa_fold_smoke()

    print("\nAll tests passed.")


if __name__ == "__main__":
    main()
