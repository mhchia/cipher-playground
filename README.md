# cipher-playground

Learning exercises for lattice-based cryptography, ZK proofs, and folding schemes.

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
notebook/           Jupyter notebooks — exploration, coursework, homework
  algebra/          Visual Algebra, MoonMath, polynomial rings
exercises/          .py implementations — structured learning progression
  00_sage_basics/   SageMath API cheatsheet
  03_ajtai/         Ajtai commitment (Module-SIS)
  kzg.py            KZG polynomial commitment (standalone)
```

### Exercise progression (from CLAUDE.md)

1. Polynomial arithmetic in Z_q[X]/(X^d + 1)
2. NTT / INTT
3. Ajtai commitment
4. Embedding comparison (NTT vs coefficient)
5. Sumcheck (toy version)
6. Simple folding step
7. b-ary decomposition
8. CCS satisfaction check
