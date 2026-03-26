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
notebook/              Jupyter notebooks — exploration, coursework, homework
  algebra/             Visual Algebra, MoonMath, polynomial rings
exercises/             .py implementations — structured learning progression
  00_sage_basics/      SageMath API cheatsheet
  01_ring_arith/       Polynomial ring arithmetic, NTT
  03_ajtai/            Ajtai commitment (Module-SIS)
  10_ecc/              Elliptic curve cryptography
  11_pairings/         Bilinear pairings, BLS signatures
  12_poly_commit/      KZG polynomial commitment
  13_groth16/          R1CS, QAP, Groth16
  14_plonk/            Gates, permutation argument, PLONK
  15_fri/              FRI, Reed-Solomon proximity testing
  16_zkvm/             Instruction traces, execution proofs
```

## Exercise Progression

### Lattice track (01–09)
1. Polynomial arithmetic in Z_q[X]/(X^d + 1)
2. NTT / INTT
3. Ajtai commitment
4. Embedding comparison (NTT vs coefficient)
5. Sumcheck (toy version)
6. Simple folding step
7. b-ary decomposition
8. CCS satisfaction check

### SNARK/STARK track (10–16)
9. ECC — group law, scalar multiplication
10. Pairings — bilinear maps, BLS signatures
11. Polynomial commitments — KZG
12. Groth16 — R1CS, QAP
13. PLONK — gates, permutation argument
14. FRI — Reed-Solomon proximity testing
15. zkVM — instruction traces, execution proofs
