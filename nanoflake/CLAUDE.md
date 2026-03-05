# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Graphene is a Python tool for generating graphene nanoflakes with customizable shapes and edge types. Supports hexagonal and triangular flakes with zigzag, armchair, or alternating edges. Outputs atomic structure files in XYZ format.

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
  graph_nanoflake.py    # Hexagonal and triangular flakes (masking approach)
  triangular.py         # Triangular zigzag flakes (lattice-cut approach)
  hexagonal.py          # Hexagonal zigzag flakes (hexagon-tiling approach)
```

## Running the Code

```bash
PYTHON=~/miniconda3/envs/graphene/bin/python

# Hexagonal flake by tiling (guaranteed zigzag, preferred)
$PYTHON Codes/hexagonal.py -n 3
$PYTHON Codes/hexagonal.py -n 3 --saturate --output coronene.xyz

# Triangular zigzag flake
$PYTHON Codes/triangular.py -n 4 --saturate

# Hexagonal/triangular via masking (legacy)
$PYTHON Codes/graph_nanoflake.py -n 4 -s hexagonal -e zigzag
```

## Architecture

### `hexagonal.py` (preferred for hexagonal flakes)

**Key functions:** `generate_hexagonal_flake(n)`, `saturate(flake)`

**Algorithm:**
1. Hexagon centres placed on triangular lattice:
   `v1 = (3a/2, a√3/2)`, `v2 = (0, a√3)` — these are the actual edge-sharing
   neighbour directions for flat-top hexagons, so shared vertices are identical
   in floating-point and deduplication via a rounded-coordinate dict is exact.
2. Axial cube-coordinate condition `|i|≤n, |j|≤n, |i+j|≤n` selects the
   hexagonal patch (`3n²+3n+1` hexagons).
3. 6 vertices per hexagon at 0°, 60°, 120°, … (flat-top).
4. Zigzag edges are structural — they follow from the hexagon tiling, not from
   a cut alignment.

### `graph_nanoflake.py` (legacy, hexagonal + triangular)

**Main function:** `create_graphene_flake(n, shape, edge_type, saturated, vacuum, orientation, visualize)`

Generates a large honeycomb supercell and masks it with a shape-specific
boundary polygon. Edge type (zigzag/armchair/alternating) controls the polygon
orientation. Works for both hexagonal and triangular shapes.

**Python import:**
```python
from Codes.graph_nanoflake import create_graphene_flake
flake = create_graphene_flake(n=4, shape='hexagonal', edge_type='zigzag',
                               saturated=True, orientation='xy')

from Codes.hexagonal import generate_hexagonal_flake, saturate
flake = saturate(generate_hexagonal_flake(n=3))
```

### `triangular.py`

Lattice-cut approach for triangular zigzag flakes. Valid sizes: n not divisible
by 3 (n = 4, 5, 7, 8, …). Reports a warning for invalid sizes.

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
