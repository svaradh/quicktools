# Graphene Nanoflake Generator

A Python tool for generating graphene nanoflakes with customizable shapes and edge types. Supports hexagonal and triangular flakes with zigzag, armchair, or alternating edges. Outputs atomic structure files in XYZ format.

## Setup

```bash
# Create and activate conda environment
conda env create -f environment.yml
conda activate graphene
```

Or manually:

```bash
conda create -n graphene python=3.11 -y
pip install ase numpy
```

**Dependencies:** Python 3, ASE (Atomic Simulation Environment), NumPy

## Usage

```bash
python Codes/graph_nanoflake.py [options]
```

### Options

| Flag | Description | Default |
|------|-------------|---------|
| `-n`, `--size N` | Size parameter — atoms along each edge | `4` |
| `-s`, `--shape` | `hexagonal` or `triangular` | `hexagonal` |
| `-e`, `--edge-type` | `zigzag`, `armchair`, or `alternating` | `zigzag` |
| `--no-saturate` | Disable hydrogen passivation | — |
| `-o`, `--orientation` | Plane: `xy`, `xz`, `yz` | `xy` |
| `-v`, `--visualize` | Open ASE GUI viewer | — |
| `--output FILE` | Output file path | `graphene_flake.xyz` |

### Examples

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

### Python API

```python
from Codes.graph_nanoflake import create_graphene_flake

flake = create_graphene_flake(
    n=4,
    shape='hexagonal',
    edge_type='zigzag',
    saturated=True,
    orientation='xy',
    visualize=True
)
```

## Output

The XYZ file contains atomic coordinates readable by:
- ASE (`ase.io.read`)
- VMD, VESTA, Avogadro, and other molecular viewers
- Quantum chemistry packages (VASP, Gaussian, etc. via ASE conversion)

## Physical Parameters

| Parameter | Value |
|-----------|-------|
| C–C bond length | 1.42 Å |
| C–H bond length | 1.09 Å |
| Unit cell width | 2.46 Å |
| Bond detection cutoff | 1.5 Å |

