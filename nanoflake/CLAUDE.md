# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Python tools for generating hexagonal and triangular graphene nanoflakes with zigzag edges. Only zigzag edges are supported. Outputs atomic structure files in XYZ format.

## Environment Setup

```bash
# Create and activate conda environment
~/miniconda3/bin/conda create -n graphene python=3.11 -y
~/miniconda3/envs/graphene/bin/pip install ase numpy
```

**Dependencies:** Python 3, ASE (Atomic Simulation Environment), NumPy

## Project Structure

```
Codes/
  graph_nanoflake.py    # Unified entry point — dispatches to hexagonal.py / triangular.py
  hexagonal.py          # Hexagonal zigzag flakes (hexagon-tiling approach)
  triangular.py         # Triangular zigzag flakes (lattice-cut approach)
```

## Running the Code

Only zigzag edges are produced.

```bash
PYTHON=~/miniconda3/envs/graphene/bin/python

# Unified entry point (recommended)
$PYTHON Codes/graph_nanoflake.py -s hexagonal -n 3 --saturate
$PYTHON Codes/graph_nanoflake.py -s triangular -n 4 --truncate-corners --saturate

# Backends can also be run directly
$PYTHON Codes/hexagonal.py -n 3 --saturate
$PYTHON Codes/triangular.py -n 4 --saturate
```

## Architecture

### `graph_nanoflake.py` (unified entry point)

Thin CLI dispatcher. Imports from `hexagonal.py` and `triangular.py` and routes
based on `-s`. No geometry logic of its own.

**Python import:**
```python
from Codes.hexagonal import generate_hexagonal_flake, saturate
flake = saturate(generate_hexagonal_flake(n=3))

from Codes.triangular import generate_triangular_flake, saturate, truncate_corners
flake = saturate(truncate_corners(generate_triangular_flake(n=4)))
```

### `hexagonal.py`

**Key functions:** `generate_hexagonal_flake(n)`, `saturate(flake)`

**Algorithm:**
1. Hexagon centres placed on triangular lattice with primitive vectors
   `v1 = (3a/2, a√3/2)`, `v2 = (0, a√3)` — the actual edge-sharing neighbour
   directions for flat-top hexagons, so shared vertices are identical in
   floating-point and deduplication via a rounded-coordinate dict is exact.
2. Axial cube-coordinate condition `|i|≤n, |j|≤n, |i+j|≤n` selects the
   hexagonal patch (`3n²+3n+1` hexagons).
3. 6 vertices per hexagon at 0°, 60°, 120°, … (flat-top).
4. Zigzag edges are structural — they follow from the hexagon tiling geometry.

### `triangular.py`

**Key functions:** `generate_triangular_flake(n)`, `check_zigzag_edges(flake)`,
`truncate_corners(flake)`, `saturate(flake)`

Lattice-cut approach: builds a graphene sheet via ASE, cuts a triangular region
using a cross-product point-in-triangle test, and validates zigzag edges via
sublattice analysis. Valid zigzag sizes: n not divisible by 3 (n = 4, 5, 7, 8, …).

## Physical constants

- C-C bond length: 1.42 Å
- C-H bond length: 1.09 Å
- Lattice constant: 2.46 Å (√3 × C-C)
- Bond detection cutoff: 1.6 Å

## Output

The XYZ file contains atomic coordinates readable by:
- ASE (`ase.io.read`)
- VMD, VESTA, Avogadro, and other molecular viewers
- Quantum chemistry packages (VASP, Gaussian, etc. via ASE conversion)
