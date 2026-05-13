from sage.all import *

from relations import LinRelation, LinInstance, LinWitness

from ring import Rq


# TODO: use batch+ in p.34 to reduce n further? It requires another sumcheck.

def rok_batch(lin: LinRelation, n_target_eval_rows: int) -> LinRelation:
    """
    Batch evaluation statements into smaller statements.
    E.g. \tilde f(r)      = s
         \tilde f(\bar r) = \bar s
    -> tilde f(r) + c \tilde f(\bar r) = s + c \bar s
    """

    H = lin.instance.H
    Y = lin.instance.Y
    n_top = lin.n_top
    n_orig_rows = lin.hat_n - n_top

    #
    # Verifier
    #
    c = Rq.random_element()
    c_matrix = matrix(
        Rq,
        n_target_eval_rows,
        n_orig_rows,
        [c ** i for i in range(n_target_eval_rows * n_orig_rows)]
    )
    # Send `c_matrix` to Prover

    #
    # Both
    #
    # Calculate H_tilde = [H_top            ]
    #                     [c_matrix * H_bot]
    H_top = H[:n_top]
    H_bot = H[n_top:]
    assert H_bot.nrows() == n_orig_rows
    new_H_bot = c_matrix * H_bot
    H_tilde = H_top.stack(new_H_bot)

    # Calculate Y_tilde = [Y_top            ]
    #                     [c_matrix * Y_bot]
    Y_top = Y[:n_top]
    assert H_top.nrows() == Y_top.nrows()
    Y_bot = Y[n_top:]
    assert Y_bot.nrows() == n_orig_rows
    Y_tilde = Y_top.stack(c_matrix * Y_bot)
    assert H_tilde.nrows() == Y_tilde.nrows()

    return LinRelation(
        LinInstance(
            H=H_tilde,
            F_com=lin.instance.F_com,
            F_eval=lin.instance.F_eval,
            Y=Y_tilde,
            beta=lin.beta,
        ),
        LinWitness(W=lin.witness.W)
    )
