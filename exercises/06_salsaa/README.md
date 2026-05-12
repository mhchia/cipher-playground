# SALSAA

Toy Sage implementation of the SALSAA paper's folding rok chain.

```
Π^join → Π^norm → Π^⊗RP → Π^fold → Π^join → Π^batch → Π^b-decomp
```

## Layout

```
06_salsaa/
├── rok/             # RoK sub-protocols (one per file)
│   ├── join.py      # Π^join — combine relations sharing F_com
│   ├── norm.py      # Π^norm — rok_norm + rok_bar_sum (sumcheck-driven)
│   ├── rp.py        # Π^⊗RP  — JL tensor random projection (W → ŵ)
│   ├── fold.py      # Π^fold — collapse W width r → 1 via RLC
│   ├── batch.py     # Π^batch — collapse H eval rows via Vandermonde RLC
│   └── decompose.py # Π^b-decomp — balanced b-ary witness decomposition
│
├── salsaa.py        # top-level driver: assembles the full chain (salsaa.fold)
├── tests.py         # unit + smoke tests for every module
```

## How to run

Requires [SageMath](https://www.sagemath.org/) ≥ 10.

```bash
sage -python exercises/06_salsaa/tests.py
```


## References

- [SALSAA-note-by-Yingfei-1.md](https://github.com/coset-io/baby-lattice-folding/blob/main/notes-paper/4-SALSAA-I-2026-03-01.pdf)
- [SALSAA-note-by-Yingfei-2.md](https://github.com/coset-io/baby-lattice-folding/blob/main/notes-paper/4-SALSAA-II-2026-03-08.pdf)
- [SALSAA paper](https://eprint.iacr.org/2025/2124) (main protocol reference)
