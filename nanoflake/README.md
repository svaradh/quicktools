# Graphene Nanoflake Generator

Python tools for generating graphene nanoflakes with customizable shapes and edge types. Outputs atomic structure files in XYZ format readable by ASE, VESTA, VMD, and Avogadro.

**Dependencies:** Python 3, ASE, NumPy

## Setup

```bash
conda env create -f environment.yml
conda activate graphene
```

Or manually:

```bash
conda create -n graphene python=3.11 -y
pip install ase numpy
```

---

## `graph_nanoflake.py` — Hexagonal and triangular flakes

```bash
python Codes/graph_nanoflake.py [options]
```

| Flag | Description | Default |
|------|-------------|---------|
| `-n`, `--size N` | Size parameter — atoms along each edge | `4` |
| `-s`, `--shape` | `hexagonal` or `triangular` | `hexagonal` |
| `-e`, `--edge-type` | `zigzag`, `armchair`, or `alternating` | `zigzag` |
| `--no-saturate` | Disable hydrogen passivation | — |
| `-o`, `--orientation` | Plane: `xy`, `xz`, `yz` | `xy` |
| `-v`, `--visualize` | Open ASE GUI viewer | — |
| `--output FILE` | Output file path | `graphene_flake.xyz` |

```bash
python Codes/graph_nanoflake.py -n 4                          # hexagonal, zigzag
python Codes/graph_nanoflake.py -n 4 -e alternating -v        # alternating edges, visualize
python Codes/graph_nanoflake.py -n 3 --no-saturate            # no H passivation
```

### Python API

```python
from Codes.graph_nanoflake import create_graphene_flake

flake = create_graphene_flake(n=4, shape='hexagonal', edge_type='zigzag',
                               saturated=True, orientation='xy', visualize=False)
```

---

## `triangular.py` — Triangular zigzag flakes

Generates triangular graphene nanoflakes with pure zigzag edges using a
point-in-triangle lattice cut. Validates edge character via sublattice analysis.

> **Note:** Valid zigzag triangular flakes only exist for `n` not divisible by 3
> (n = 4, 5, 7, 8, 10, 11, …). The script reports a warning for invalid sizes.

```bash
python Codes/triangular.py [options]
```

| Flag | Description | Default |
|------|-------------|---------|
| `-n`, `--size N` | Size parameter: side length = n × 2.46 Å | `6` |
| `--truncate-corners` | Remove 1-coordinated corner C atoms | — |
| `--saturate` | Passivate edge carbons with hydrogen | — |
| `-v`, `--visualize` | Open ASE GUI viewer | — |
| `--output FILE` | Output file path | `triangular_flake.xyz` |

```bash
python Codes/triangular.py -n 4                               # basic triangular flake
python Codes/triangular.py -n 4 --truncate-corners            # clip the three corner atoms
python Codes/triangular.py -n 4 --saturate -v                 # H-passivated, visualize
python Codes/triangular.py -n 4 --truncate-corners --saturate # truncate then passivate
```

---

## Physical Parameters

| Parameter | Value |
|-----------|-------|
| C–C bond length | 1.42 Å |
| C–H bond length | 1.09 Å |
| Lattice constant | 2.46 Å |
| Bond detection cutoff | 1.6 Å |
