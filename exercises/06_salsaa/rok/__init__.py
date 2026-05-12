"""
SALSAA Reduction-of-Knowledge (RoK) sub-protocols.

Each module implements one step from the SALSAA paper's RoK chain:

  Π^join     → join.py     — combine relations sharing F_com
  Π^norm     → norm.py     — rok_norm + rok_bar_sum (sumcheck-driven norm proof)
  Π^⊗RP      → rp.py       — JL tensor random projection
  Π^fold     → fold.py     — collapse witness width r → 1 via RLC
  Π^batch    → batch.py    — collapse H eval rows via Vandermonde RLC
  Π^b-decomp → decompose.py — balanced b-ary witness decomposition
"""

from .batch import rok_batch
from .decompose import (
    rok_decompose,
    get_l,
    balanced_b_ary_decompose_Fq,
    compose_Fq,
    decompose_W,
)
from .fold import rok_fold
from .join import rok_join
from .norm import rok_norm, rok_bar_sum
from .rp import rok_rp
