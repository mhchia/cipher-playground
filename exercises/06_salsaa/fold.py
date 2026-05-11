from sage.all import *

from relations import LinRelation, LinInstance, LinWitness

from ring import Rq


def sample_C(r_in: int, r_out):
    # challenge set: larger challenge set over subtractive set
    # rok_rp computes and proves about randomised projections of the witness.
    # This allows us to use a much larger challenge set in the folding step Π fold
    # instead of a subtractive set (used in [KLNO24]), which ultimately removes (roughly) one λ factor from the proof size of [KLNO24].

    # SALSAA paper p.36
    # We take a different approach and sample challenges so that the coefficients
    # are sampled uniformly from the ternary set {−1, 0, 1}. The set of challenges
    # with ternary coefficients is not strong samplings sets per definition,
    # but the probability of sampling elements so that the inverse of two is non-invertible is small.
    import random

    return matrix(
        Rq,
        r_in,
        r_out,
        # p=1/3, different from the one in rok_rp
        lambda i, j: random.choice([-1, 0, 1])
    )


def rok_fold(lin: LinRelation, r_out: int) -> LinRelation:
    r_in = lin.r

    # TODO: confirm it's correct.
    # gamma
    # Since c_ij is +-1 or 0,
    # |c_ij * w|^2 = |c_ij w_1|^2 + ... + |c_ij w_r|^2
    #              <= |w|^2
    expansion_factor = 1

    #
    # Verifier
    #
    C = sample_C(r_in, r_out)
    Y = lin.instance.Y
    Y_tilde = Y * C

    #
    # Prover
    #
    H, F_com, F_eval, W = lin.instance.H, lin.instance.F_com, lin.instance.F_eval, lin.witness.W
    # \tilde W = W * C
    W_tilde = W * C

    # Derive new bound
    # \beta = r_in * \gamma * \beta
    v_square = (r_in * expansion_factor) ** 2 * lin.v_square

    return LinRelation(
        instance=LinInstance(
            H=H,
            F_com=F_com,
            F_eval=F_eval,
            Y=Y_tilde,
            v_square=v_square,
        ),
        witness=LinWitness(W_tilde),
    )
