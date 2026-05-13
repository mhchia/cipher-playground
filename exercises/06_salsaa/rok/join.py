from sage.all import *

from relations import LinRelation, LinInstance, LinWitness
from ring import Rq


def rok_join(lin_0: LinRelation, lin_1: LinRelation) -> LinRelation:
    """

    """
    assert lin_0.m == lin_1.m

    def unpack(lin: LinRelation):
        inst = lin.instance
        return inst.H, inst.F_com, inst.F_eval, inst.Y, lin.witness.W
    H_0, F_0_top, F_0_bot, Y_0, W_0 = unpack(lin_0)
    H_1, F_1_top, F_1_bot, Y_1, W_1 = unpack(lin_1)

    # Check commitment matches
    assert F_0_top == F_1_top
    n_top = F_0_top.nrows()

    # NOTE: H can be non-square matrix
    H_0_bot = H_0[n_top:, n_top:]
    H_1_bot = H_1[n_top:, n_top:]

    #
    # Prover
    #
    Y_01_bot = H_0_bot * F_0_bot * W_1
    Y_10_bot = H_1_bot * F_1_bot * W_0

    # Send Y_01_bot and Y_10_bot to Verifier

    #
    # Both Prover and Verifier
    #

    # H_new = [I                    ]
    #         [     H_0_bot         ]
    #         [              H_1_bot]
    H_new = block_diagonal_matrix([identity_matrix(Rq, n_top), H_0_bot, H_1_bot])

    # F = [F_top  ]
    #     [F_0_bot]
    #     [F_1_bot]
    F_new_bot = F_0_bot.stack(F_1_bot)

    # Y = [
    #     [Y_0[:n_top], Y_1[:n_top]],
    #     [Y_0[n_top:], Y_01_bot],
    #     [Y_10_bot, Y_1[n_top:]],
    # ]
    Y_new_0 = Y_0[:n_top].augment(Y_1[:n_top])
    Y_new_1 = Y_0[n_top:].augment(Y_01_bot)
    Y_new_2 = Y_10_bot.augment(Y_1[n_top:])
    Y_new = Y_new_0.stack(Y_new_1).stack(Y_new_2)

    #
    # Check relation holds
    #

    # W = [W_0, W_1]
    W_new = W_0.augment(W_1)

    # \beta = max(\beta_1, \beta_2)
    new_beta = max(lin_0.beta, lin_1.beta)
    lin_joined = LinRelation(
        LinInstance(H_new, F_0_top, F_new_bot, Y_new, new_beta),
        LinWitness(W_new),
    )

    # lin_joined
    # .m = lin_0.m  (unchagned)
    # .hat_n = hat_n_0 + hat_n_1 - n_top
    # .r = r_0 + r_1
    return lin_joined
