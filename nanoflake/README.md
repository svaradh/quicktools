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

## `graph_nanoflake.py` — Unified entry point (zigzag only)

Dispatches to `hexagonal.py` or `triangular.py` based on `-s`. Only zigzag
edges are produced.

```bash
python Codes/graph_nanoflake.py [options]
```

| Flag | Description | Default |
|------|-------------|---------|
| `-s`, `--shape` | `hexagonal` or `triangular` | `hexagonal` |
| `-n`, `--size N` | Rings (hexagonal) or side-length units (triangular) — see below | `3` |
| `--saturate` | Passivate edge carbons with hydrogen | — |
| `--truncate-corners` | Triangular only: remove 1-coordinated corner C atoms | — |
| `--output FILE` | Output file path | `hexagonal_flake.xyz` / `triangular_flake.xyz` |
| `-v`, `--visualize` | Open ASE GUI viewer | — |

**`-n` meaning depends on shape:**
- **Hexagonal:** number of hexagon rings (n=1 → coronene, n=2 → circumcoronene)
- **Triangular:** side length in units of 2.46 Å; valid zigzag sizes are n not divisible by 3

```bash
python Codes/graph_nanoflake.py -s hexagonal -n 3 --saturate
python Codes/graph_nanoflake.py -s triangular -n 4 --truncate-corners --saturate -v
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

## `hexagonal.py` — Hexagonal zigzag flakes by hexagon tiling

Generates hexagonal graphene nanoflakes by tiling flat-top benzene-ring hexagons
outward from a central hexagon. Because the flake is built from hexagon vertices
rather than cut from a sheet, zigzag edges are guaranteed by construction for any
ring count. Named PAHs: n=0 → benzene, n=1 → coronene, n=2 → circumcoronene.

```bash
python Codes/hexagonal.py [options]
```

| Flag | Description | Default |
|------|-------------|---------|
| `-n`, `--rings N` | Number of hexagon rings around the central hexagon | `3` |
| `--saturate` | Passivate edge carbons with hydrogen | — |
| `--output FILE` | Output file path | `hexagonal_flake.xyz` |
| `-v`, `--visualize` | Open ASE GUI viewer | — |

```bash
python Codes/hexagonal.py -n 2                    # circumcoronene (54 C)
python Codes/hexagonal.py -n 3 --saturate         # H-passivated, n=3
python Codes/hexagonal.py -n 1 --saturate -v      # coronene, visualize
```

### Python API

```python
from Codes.hexagonal import generate_hexagonal_flake, saturate

flake = generate_hexagonal_flake(n=3)        # 96 C atoms
flake = saturate(flake)                      # adds 24 H atoms
```

### Algorithm

1. Place hexagon centres on a triangular lattice with primitive vectors
   `v1 = (3a/2, a√3/2)` and `v2 = (0, a√3)` (edge-sharing neighbour directions).
2. Select centres satisfying the axial cube-coordinate condition
   `|i| ≤ n`, `|j| ≤ n`, `|i+j| ≤ n` — this gives `3n² + 3n + 1` hexagons.
3. Compute 6 vertices per hexagon at angles 0°, 60°, 120°, … (flat-top).
4. Deduplicate via a rounded-coordinate dict (shared vertices are exact in
   floating-point with these lattice vectors).
5. Optionally passivate: each edge C (coordination 2) gets 1 H opposite its
   C–C bisector; each corner C (coordination 1) gets 2 H at ±120°.

---

## Physical Parameters

| Parameter | Value |
|-----------|-------|
| C–C bond length | 1.42 Å |
| C–H bond length | 1.09 Å |
| Lattice constant | 2.46 Å |
| Bond detection cutoff | 1.6 Å |
