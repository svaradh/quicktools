# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Repository Overview

**quicktools** is a collection of independent computational chemistry scripts. Each subdirectory is a self-contained tool:

- `nanoflake/` — Graphene nanoflake generator (ASE + NumPy). Has its own `CLAUDE.md` with full details.
- `constrained_dft/` — DFT charge-transfer analysis between two molecular fragments (ASE + NWChem).

## Environment Setup

### nanoflake
```bash
conda env create -f nanoflake/environment.yml   # creates 'graphene' env
conda activate graphene
# or manually: pip install ase numpy
```

### constrained_dft
Requires NWChem installed and on PATH. Python dependencies: ASE, NumPy.
```bash
pip install ase numpy
```

## Running the Tools

### nanoflake
```bash
# From repo root, using graphene conda env
~/miniconda3/envs/graphene/bin/python nanoflake/Codes/graph_nanoflake.py [options]

# Key options: -n SIZE, -s hexagonal|triangular, -e zigzag|armchair|alternating,
#              --no-saturate, -o xy|xz|yz, -v (visualize), --output FILE
```
See `nanoflake/CLAUDE.md` for full CLI reference and Python API.

### constrained_dft
```bash
cd constrained_dft
python charge_transfer.py frag1.xyz frag2.xyz \
    --frag1-charge CHARGE --frag1-mult MULT \
    --total-charge CHARGE --total-mult MULT \
    [-np NPROCS] [--basis BASIS] [--xc XC] [--vacuum VACUUM] [--scf-iterations N]

# Example (acetate radical + hydrogen atom):
python charge_transfer.py acetate.xyz proton.xyz \
    --frag1-charge 0 --frag1-mult 2 \
    --total-charge 0 --total-mult 1 --frag2-mult 2
```

## Architecture

### constrained_dft/charge_transfer.py

Single-file script with no importable API. Flow:
1. Parses CLI args; derives `frag2_charge = total_charge - frag1_charge`
2. Reads both XYZ files via ASE, concatenates into a combined system, and calls `combined.center(vacuum=...)` — this shared cell ensures all three DFT calculations run on an identical grid
3. Runs three NWChem DFT single-point calculations (combined, frag1, frag2) via a custom `NWChemWithDplot` subclass that appends a `DPLOT` block to each NWChem input file to dump electron density cube files
4. Computes formation energy: `dE = E_combined - E_frag1 - E_frag2`
5. Reads the three cube files and writes `nwchem_scratch/density_difference.cube` = ρ(combined) − ρ(frag1) − ρ(frag2)

NWChem scratch files land in `nwchem_scratch/` (created automatically).

### nanoflake/Codes/graph_nanoflake.py

See `nanoflake/CLAUDE.md` for the full architecture description.

## Output Files

- `nanoflake`: XYZ file (default: `graphene_flake.xyz`) readable by ASE, VESTA, VMD, Avogadro
- `constrained_dft`: `nwchem_scratch/density_difference.cube` for visualization in VESTA/VMD; formation energy printed to stdout in eV and kcal/mol
