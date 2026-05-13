from sage.all import *

from relations import LinRelation, LinInstance, LinWitness

from ring import Rq, x, to_centered, Fq, d


def get_l(beta: int, b: int) -> int:
    # get l = \ceil(\log {2 beta + 1})
    assert beta > 0
    v = 2 * beta + 1
    l = 0
    while v > 0:
        v = v // b
        l += 1
    return l


def balanced_b_ary_decompose_Fq(f: Fq, b: int, l: int) -> list[Fq]:
    """
    E.g. b = 2. f = 7
        -> [ 1,  1,  1,  0, ...]
    E.g. b = 2. f = -7
        -> [-1, -1, -1,  0, ...]
    E.g. b = 3. f = 5
        -> [-1, -1,  1,  0, ...]
    """
    assert b > 0
    assert l > 0
    f_ct = to_centered(f)
    sign = -1 if f_ct < 0 else 1
    f_abs = abs(f_ct)

    coeffs = []
    for _ in range(l):
        r = f_abs % b
        # We want each digit \in [-b/2, b/2]
        if r > b // 2:
            # r > b/2..., need to make it back in range
            # Make this digit becomes negative so it's within [-b/2, 0]
            r = r - b
            # Since we deduct r by b, add it back to f_abs as "carry"
            f_abs += b
        coeffs.append(Fq(sign * r))
        f_abs = f_abs // b
    return coeffs


def compose_Fq(coeffs: list[Fq], b: int) -> Fq:
    return sum([Fq(c) * (b ** i) for i, c in enumerate(coeffs)])


def decompose_W(W: matrix, b: int, l: int) -> list[matrix]:
    # Store V_i matrices in list for now. We'll transform it back
    # to matrices in the end
    V = [list() for _ in range(l)]
    for r in W.list():
        # r = 4 + 5x + 3x^2
        cur_bit_polys = [Rq(0) for _ in range(l)]
        for exp, c in enumerate(r.list()):
            # c*x^{exp}
            c_decomposed = balanced_b_ary_decompose_Fq(c, b, l)
            # c*x^i = (d_0*b^0 + d_1*b^1 + ...) * x^i
            #       = d_0*b^0*x^i + d_1*b^1*x^i + ...
            #            V_0            V_1     ...
            for b_exp, d in enumerate(c_decomposed):
                # V[b_exp].append()
                cur_bit_polys[b_exp] += d * (x ** exp)
        for k in range(l):
            V[k].append(cur_bit_polys[k])
    # Transform V_i back to matrices
    return [matrix(Rq, W.nrows(), W.ncols(), v) for v in V]


def rok_decompose(lin: LinRelation, b: int) -> LinRelation:
    # beta = \sqrt(v_square)
    beta = isqrt(lin.v_square)
    l = get_l(beta, b)

    #
    # Prover
    #
    H, F, W = lin.instance.H, lin.instance.F, lin.witness.W
    Vs = decompose_W(W, b, l)
    Zs = [H * F * V_k for V_k in Vs]

    # V_tilde = [V_0 || ... || V_{l-1}]
    V_tilde = Vs[0]
    for V_k in Vs[1:]:
        V_tilde = V_tilde.augment(V_k)

    # Send Z_k, k \in [l] to Verifier

    #
    # Verifier
    #
    Y = lin.instance.Y
    # Check Y ?= \sum_{i=0}^{l-1} b^i Z_i
    Y_recovered = sum([(b**i) * Z_k for i, Z_k in enumerate(Zs)])
    assert Y == Y_recovered

    #
    # Both
    #
    # Z_tilde = [Z_0 || ... || Z_{l-1}]
    Z_tilde = Zs[0]
    for Z_k in Zs[1:]:
        Z_tilde = Z_tilde.augment(Z_k)

    new_v_square = lin.m * d * ((b // 2)**2)
    # Verfy H*F*Z_tilde = Y_tilde
    return LinRelation(
        instance=LinInstance(
            H=H,
            F_com=lin.instance.F_com,
            F_eval=lin.instance.F_eval,
            Y=Z_tilde,
            v_square=new_v_square,
        ),
        witness=LinWitness(W=V_tilde),
    )
