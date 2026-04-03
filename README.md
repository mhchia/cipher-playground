# cipher-playground

Learning exercises for lattice-based cryptography, ZK proof systems, and folding schemes.

## Setup

This project uses **SageMath 10.7** (bundles Python 3.13).

```bash
# Run a script
sage -python exercises/03_ajtai/ajtai.py

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
  01_ring_arith/       Polynomial ring arithmetic, NTT
  02_field_ext/        Field extensions (Fp²)
  10_ecc/              Elliptic curve cryptography
  11_pairings/         High-level usage of pairing, G_1, G_2, GT
  12_kzg/              KZG polynomial commitment
  13_groth16/          R1CS, QAP, Groth16
```

## Exercise Progression

### Lattice track (01–09)
- [x] 01 — Polynomial arithmetic in Z_q[X]/(X^d + 1)
- [ ] 01 — NTT / INTT (in progress)
- [x] 02 — Field extensions (Fp²)
- [ ] 03 — Ajtai commitment (in progress)
- [ ] 04 — Embedding comparison (NTT vs coefficient)
- [ ] 05 — Sumcheck (toy version)
- [ ] 06 — Simple folding step
- [ ] 07 — b-ary decomposition
- [ ] 08 — CCS satisfaction check

### SNARK/STARK track (10–16)
- [x] 10 — ECC — group law, scalar multiplication
- [x] 11 — Pairings — High-level usage of pairing, G_1, G_2, GT
- [x] 12 — KZG
- [x] 13 — Groth16 — simplified: no ZK, no public input verification
- [ ] 14 — PLONK — gates, permutation argument
- [ ] 15 — FRI — Reed-Solomon proximity testing
- [ ] 16 — zkVM — instruction traces, execution proofs
