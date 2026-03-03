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
  graph_nanoflake.py    # Main nanoflake generator
```

## Running the Code

```bash
# Using the conda environment
~/miniconda3/envs/graphene/bin/python Codes/graph_nanoflake.py [options]
```

**CLI options:**
```
-n, --size N            Size parameter - atoms along edge (default: 4)
-s, --shape             hexagonal, triangular (default: hexagonal)
-e, --edge-type         zigzag, armchair, alternating (default: zigzag)
--no-saturate           Disable hydrogen saturation
-o, --orientation       Plane: xy, xz, yz (default: xy)
-v, --visualize         Open ASE GUI viewer
--output FILE           Output file (default: graphene_flake.xyz)
```

**Examples:**
```bash
# Hexagonal flake with zigzag edges (default)
python Codes/graph_nanoflake.py -n 4

# Triangular flake with armchair edges
python Codes/graph_nanoflake.py -n 5 -s triangular -e armchair

# Hexagonal with alternating edges, visualize
python Codes/graph_nanoflake.py -n 4 -e alternating -v

# No hydrogen passivation, custom output
python Codes/graph_nanoflake.py -n 3 --no-saturate --output my_flake.xyz
```

**Python import:**
```python
from Codes.graph_nanoflake import create_graphene_flake
flake = create_graphene_flake(n=4, shape='hexagonal', edge_type='zigzag',
                               saturated=True, orientation='xy', visualize=True)
```

## Architecture

**Main module:** `Codes/graph_nanoflake.py`

**Main function:** `create_graphene_flake(n, shape, edge_type, saturated, vacuum, orientation, visualize)`

- `n`: Size parameter (atoms along each edge)
- `shape`: `'hexagonal'` or `'triangular'`
- `edge_type`: `'zigzag'`, `'armchair'`, or `'alternating'` (hexagonal only)
- `saturated`: When True, passivates edge carbons with hydrogen
- `vacuum`: Spacing around structure in Angstroms (default: 5.0)
- `orientation`: Plane orientation - `'xy'`, `'xz'`, or `'yz'`
- `visualize`: When True, opens structure in ASE GUI viewer

**Algorithm flow:**
1. Generates honeycomb lattice positions from unit cell vectors
2. Applies shape-specific boundary mask (hexagon or triangle)
3. Mask geometry depends on edge type (zigzag/armchair/alternating)
4. Removes duplicate atoms at boundaries
5. Optionally saturates edge carbons with H atoms
6. Rotates to desired plane orientation
7. Returns ASE `Atoms` object

**Physical constants:**
- C-C bond length: 1.42 Å
- C-H bond length: 1.09 Å
- Unit cell width: 2.46 Å (√3 × C-C bond)
- Bond detection cutoff: 1.5 Å

**Known issues:**
- Armchair triangular flakes may need further refinement

## Output

The XYZ file contains atomic coordinates readable by:
- ASE (`ase.io.read`)
- VMD, VESTA, Avogadro, and other molecular viewers
- Quantum chemistry packages (VASP, Gaussian, etc. via ASE conversion)
