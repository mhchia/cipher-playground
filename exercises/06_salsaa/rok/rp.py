from sage.all import *

from relations import LinRelation, LinInstance, LinWitness

from ring import Rq

def _sample_J_entry():
    import random
    # χ(0)=1/2, χ(1)=χ(-1)=1/4
    return Rq(random.choice([0, 0, 1, -1]))


def sample_J(n_rp: int, m_rp: int) -> matrix:
    _m = []
    for _ in range(n_rp):
        _m.append([_sample_J_entry() for _ in range(m_rp)])
    return matrix(Rq, _m)


def vec(W: matrix):
    """
    column-major
    E.g. W = [[1,2,3],
              [4,5,6]]
         vec(W) = [1,4,2,5,3,6]
    """
    return vector(W.transpose().list())


def rok_rp(lin: LinRelation, n_rp: int, m_rp: int) -> tuple[LinRelation, LinRelation]:
    """
    Prover proves W satisfying FW=Y with |w| <= beta using Johnson-Lindenstrauss
    random projection.
    Returns
    - lin_orig:  original `lin` plus JL
    - lin_w_hat: \hat w commitment and its low norm
    """

    n_hat, m, n, n_top, r, F_com, F_eval, H, Y = lin.hat_n, lin.m, lin.n, lin.n_top, lin.r, lin.instance.F_com, lin.instance.F_eval, lin.instance.H, lin.instance.Y
    W = lin.witness.W

    assert m_rp == n_rp * r, f"need m_rp = n_rp·r, got m_rp={m_rp}, n_rp={n_rp}, r={r}"

    #
    # Verifier
    #
    # 1. Samples J and send it to Prover
    m_prime = m // r
    J = sample_J(n_rp, m_rp)

    #
    # Prover
    #

    # 1. Calculate \hat J
    size_I = m // m_rp
    I = identity_matrix(Rq, size_I)
    J_hat = I.tensor_product(J)
    assert J_hat.nrows() == size_I * n_rp
    assert J_hat.ncols() == size_I * m_rp

    # 2. Calculate \hat W \in R^{m_prime \time r}
    W_hat = J_hat * W

    # 3. flatten \hat W to \hat w \in R^{m}
    w_hat = vec(W_hat)
    assert len(w_hat) == m

    # 4. Commit \hat w to save proof size
    z_bar = F_com * w_hat

    # Sends `z_bar` to Verifier

    #
    # Verifier
    #
    # 4. Sample c from R_q
    c = Rq.random_element()
    # Sends c to Prover

    #
    # Prover
    #
    # 5-8. Calculate c_0, c_1, and \vec c_0, \vec c_1,
    c_0, c_1 = c ** m_prime, c
    c_0_vec = vector(Rq, [c_0 ** i for i in range(r)])
    c_1_vec = vector(Rq, [c_1 ** i for i in range(m_prime)])
    # \vec c = \vec c_0 ⊗ \vec c_1
    c_vec = c_0_vec.tensor_product(c_1_vec).list()
    assert len(c_vec) == m

    # 9. Calculate \vec r
    r_vec = c_1_vec * W_hat
    # Send `r_vec` to Verifier

    #
    # Both Prover and Verifier
    #
    # 10.1. H_tilde [F              ] W = [Y]
    #               [c_1_vec * J_hat]     [r]
    H_tilde = block_matrix(Rq, [
        [H,                     zero_matrix(Rq, n_hat, 1)],
        [zero_matrix(Rq, 1, n), identity_matrix(Rq, 1)],
    ])
    assert H_tilde.nrows() == n_hat + 1
    assert H_tilde.ncols() == n + 1

    F_eval_tilde = F_eval.stack(c_1_vec * J_hat)

    Y_tilde = Y.stack(r_vec)
    assert Y_tilde.nrows() == n_hat + 1
    assert Y_tilde.ncols() == r

    lin_orig = LinRelation(
        LinInstance(H=H_tilde, F_com=F_com, F_eval=F_eval_tilde, Y=Y_tilde, v_square=lin.v_square),
        LinWitness(W=W),
    )

    # 10.1. I [F_top] w_hat = [z_bar      ]
    #         [c_vec]         [c_0_vec * r]

    F_eval_hat = matrix(Rq, [c_vec])
    y_hat = column_matrix(Rq, list(z_bar) + [c_0_vec * r_vec])
    assert y_hat.nrows() == n_top + 1
    assert y_hat.ncols() == 1

    # \hat \beta = m_rp * \beta
    v_square_new = m_rp**2 * lin.v_square

    lin_w_hat = LinRelation(
        LinInstance(
            H=identity_matrix(Rq, n_top+1),
            F_com=F_com,
            F_eval=F_eval_hat,
            Y=y_hat,
            v_square=v_square_new,
        ),
        LinWitness(W=matrix(Rq, w_hat).transpose()),
    )

    return lin_orig, lin_w_hat
