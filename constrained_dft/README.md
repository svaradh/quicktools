# quicktools

A collection of computational chemistry tools built with ASE and NWChem.

## Tools

### charge_transfer.py

DFT study of charge transfer between two molecular fragments using ASE + NWChem.

Given two fragments as XYZ files, the script:
1. Computes the energy of the combined system
2. Computes the energy of each isolated fragment
3. Calculates the formation energy of the combined system relative to isolated fragments
4. Generates electron density difference cube files (rho_combined - rho_frag1 - rho_frag2)

All three calculations use an identical grid (same cell and origin) so the density difference is physically meaningful.

**Usage:**

```bash
python charge_transfer.py frag1.xyz frag2.xyz \
    --frag1-charge CHARGE --frag1-mult MULT \
    --total-charge CHARGE --total-mult MULT \
    [--frag2-mult MULT] [-np NPROCS] [--basis BASIS] [--xc XC] [--vacuum VACUUM] \
    [--scf-iterations ITERATIONS]
```

**Required arguments:**

| Argument | Description |
|---|---|
| `frag1.xyz` | XYZ file for fragment 1 |
| `frag2.xyz` | XYZ file for fragment 2 |
| `--frag1-charge` | Charge on fragment 1 |
| `--frag1-mult` | Spin multiplicity of fragment 1 |
| `--total-charge` | Total charge of the combined system |
| `--total-mult` | Total spin multiplicity of the combined system |

**Optional arguments:**

| Argument | Default | Description |
|---|---|---|
| `--frag2-mult` | same as frag1 | Spin multiplicity of fragment 2 |
| `-np` | 1 | Number of MPI processors for NWChem |
| `--basis` | 6-311++G** | Basis set |
| `--xc` | B3LYP | Exchange-correlation functional |
| `--vacuum` | 6.0 | Vacuum padding in angstroms |
| `--scf-iterations` | 500 | Maximum number of SCF iterations |

Fragment 2 charge is derived automatically as `total_charge - frag1_charge`.

**Example** (acetate radical + hydrogen atom):

```bash
python charge_transfer.py acetate.xyz proton.xyz \
    --frag1-charge 0 --frag1-mult 2 \
    --total-charge 0 --total-mult 1 \
    --frag2-mult 2
```

**Dependencies:** Python 3, NumPy, ASE, NWChem
