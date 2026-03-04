# quicktools

A collection of independent computational chemistry scripts.

## Tools

### `nanoflake/` — Graphene Nanoflake Generator

Generates graphene nanoflakes (hexagonal and triangular) with configurable edge types and optional hydrogen passivation. See [`nanoflake/README.md`](nanoflake/README.md) for full usage.

**Dependencies:** Python 3, ASE, NumPy

```bash
conda env create -f nanoflake/environment.yml
conda activate graphene
python nanoflake/Codes/graph_nanoflake.py -n 4
python nanoflake/Codes/triangular.py -n 4 --saturate -v
```

### `constrained_dft/` — Constrained DFT Charge-Transfer Analysis

Computes charge transfer between two molecular fragments using NWChem DFT single-point calculations and electron density cube file differencing.

**Dependencies:** Python 3, ASE, NumPy, NWChem (on PATH)

```bash
cd constrained_dft
python charge_transfer.py frag1.xyz frag2.xyz \
    --frag1-charge 0 --frag1-mult 2 \
    --total-charge 0 --total-mult 1
```
