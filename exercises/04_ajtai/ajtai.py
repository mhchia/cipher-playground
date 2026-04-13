from sage.all import *

# ============================================================
# Ajtai Commitment without hiding, based on M-SIS
# ============================================================


beta = 1

# toy
q_toy = 17
d_toy = 4
# falcon
q_falcon = 12289
d_falcon = 1024

# κ = commitment vector length (rows of A)
# m = witness vector length (columns of A)
kappa = 2
m = 4


def _to_centered(z: int, Rq):
    _q = Rq.base_ring().order()
    # lift to Z
    _z = ZZ(z)
    if _z <= _q // 2:
        return _z
    else:
        return _z - _q

def _l_inf_norm_Rq(r, Rq) -> int:
    _R = Rq.modulus().parent()
    lifted = _R(r.lift())
    coeffs = lifted.list()
    return max([abs(_to_centered(c, Rq)) for c in coeffs])


def _l_inf_norm_vec(w: vector, Rq) -> int:
    if len(w) == 0:
        return 0
    # return the max coeff among all Rqs in w
    return max([_l_inf_norm_Rq(r, Rq) for r in w])


def _gen_random_low_norm_poly(Rq, beta: int):
    import random
    x = Rq.gen()
    d = Rq.modulus().degree()
    # avoid 0
    poly = 0
    while poly == 0:
        poly = sum([random.randint(-beta, beta) * x**i for i in range(d)])
    return poly



def _gen_random_low_norm_witness(Rq, m, beta: int) -> vector:
    w = vector(Rq, [_gen_random_low_norm_poly(Rq, beta) for _ in range(m)])
    # Sanity check
    w_norm = _l_inf_norm_vec(w, Rq)
    assert w_norm <= beta, f"|w|_\infty is not low-norm: {w_norm=}, {beta=}"
    return w


def setup(kappa: int, m: int, Rq):
    # Vector spaces
    MS = MatrixSpace(Rq, kappa, m)
    return MS.random_element()


def commit(A: vector, z: vector) -> vector:
    """
    A: k x m, random matrix from setup
    z: m x 1, witness

    return c = Az: k x 1
    """
    return A * z


def test_ajtai(q: int, d: int):
    Fq = GF(q)
    R = PolynomialRing(Fq, 'X')
    X = R.gen()
    Rq = R.quotient(X**d + 1, 'x')
    x = Rq.gen()
    print(f"Testing ajtai with {q=}, {d=}")
    A = setup(kappa, m, Rq)
    print(f"{A=}")
    # witness
    z1 = _gen_random_low_norm_witness(Rq, m, beta)
    print(f"{z1=}")
    # Test: commit
    c1 = commit(A, z1)
    print(f"{c1=}")

    # Test: addition homomorphism
    z2 = _gen_random_low_norm_witness(Rq, m, beta)
    c2 = commit(A, z2)

    z12 = z1 + z2
    c12 = commit (A, z12)
    assert c12 == c1 + c2, f"homomorphism doesn't hold: {c1=}, {c2=}, {c12=}"

    # Test: norm growth w.h.p
    norm_z12 = _l_inf_norm_vec(z12, Rq)
    print(f"{norm_z12=}")


def main():
    test_ajtai(q_toy, d_toy)
    # test_ajtai(q_falcon, d_falcon)


if __name__ == '__main__':
    main()
