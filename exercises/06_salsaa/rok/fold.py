
from sage.all import *

from relations import LinRelation, LinInstance, LinWitness

from ring import Rq, x, d


def _sample_small_Rq():
    import random
    return sum([
        # p=1/3, different from the one in rok_rp
        random.choice([-1, 0, 1]) * (x**i)
        for i in range(d)
    ])


def sample_C(r_in: int, r_out: int) -> matrix:
    # challenge set: larger challenge set over subtractive set
    # rok_rp computes and proves about randomised projections of the witness.
    # This allows us to use a much larger challenge set in the folding step Π fold
    # instead of a subtractive set (used in [KLNO24]), which ultimately removes (roughly) one λ factor from the proof size of [KLNO24].

    # SALSAA paper p.36
    # We take a different approach and sample challenges so that the coefficients
    # are sampled uniformly from the ternary set {−1, 0, 1}. The set of challenges
    # with ternary coefficients is not strong samplings sets per definition,
    # but the probability of sampling elements so that the inverse of two is non-invertible is small.
    return matrix(
        Rq,
        r_in,
        r_out,
        [
            _sample_small_Rq()
            for _ in range(r_in * r_out)
        ],
    )


def rok_fold(lin: LinRelation, r_out: int) -> LinRelation:
    r_in = lin.r

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

    # Since every entry of C, C_ij, is sampled by every coefficient
    # being ternary, {-1,0,1}.
    # Since C_ij is low-norm, e.g. 1+x+...+x^{d-1}
    # |C_ij * W_i| <= d * |W_i| = d * \beta
    # |W| = |\sum_{i=1}^r C_ij * W_i| <= r * d * \beta
    new_beta = r_in * d * lin.beta

    return LinRelation(
        instance=LinInstance(
            H=H,
            F_com=F_com,
            F_eval=F_eval,
            Y=Y_tilde,
            beta=new_beta,
        ),
        witness=LinWitness(W_tilde),
    )
