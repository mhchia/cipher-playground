# cipher-playground

Learning exercises for lattice-based cryptography, ZK proof systems, and folding schemes.

## Setup

This project uses **SageMath 10.7** (bundles Python 3.13).

### Install sage
Make sure you're installed sage.
```bash
# Install SageMath (macOS)
brew install --cask sage
sage --version
```

### Run a script

```bash
# Run a script
sage exercises/04_ajtai/ajtai.py

# Open Jupyter
sage -sh -c "jupyter lab"

# Install the SageMath Jupyter kernel (one-time)
sage -sh -c "python -m ipykernel install --user --name sagemath --display-name 'SageMath 10.7'"
```

## Structure

```
notebook/              Jupyter notebooks
  algebra/             Visual Algebra
exercises/             .py implementations — structured learning progression
  00_sage_basics/      SageMath API cheatsheet
  01_ring_arith/       Polynomial ring arithmetic, norm
  02_ntt/              NTT / INTT
  03_field_ext/        Field extensions (Fp²)
  04_ajtai/            Ajtai commitment (Module-SIS)
  05_sumcheck/         Multivariate sumcheck protocol
  06_salsaa/           SALSAA lattice folding scheme (full RoK chain)
  10_ecc/              Elliptic curve cryptography
  11_pairings/         High-level usage of pairing, G_1, G_2, GT
  12_kzg/              KZG polynomial commitment
  13_groth16/          R1CS, QAP, Groth16
```

## Exercise Progression

### Lattice track (01–09)
- [x] 01 — Polynomial arithmetic in Z_q[X]/(X^d + 1), norm
- [x] 02 — NTT / INTT (standard + negacyclic)
- [x] 03 — Field extensions (Fp²)
- [x] 04 — Ajtai commitment (M-SIS), SIS experiment
- [x] 05 — Sumcheck (toy multivariate sumcheck protocol)
- [x] 06 — SALSAA — full lattice folding RoK chain (`join → norm → ⊗RP → fold → join → batch → b-decomp`)

### SNARK/STARK track (10–16)
- [x] 10 — ECC — group law, scalar multiplication
- [x] 11 — Pairings — High-level usage of pairing, G_1, G_2, GT
- [x] 12 — KZG
- [x] 13 — Groth16 — simplified: no ZK, no public input verification
- [ ] 14 — PLONK — gates, permutation argument
- [ ] 15 — FRI — Reed-Solomon proximity testing
- [ ] 16 — zkVM — instruction traces, execution proofs
