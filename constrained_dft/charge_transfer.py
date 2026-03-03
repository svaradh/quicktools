#!/usr/bin/env python3
"""
DFT study of charge transfer between two molecular fragments using ASE + NWChem.

Reads two fragments from XYZ files, computes:
  1. Energy of the combined system (fragment 1 + fragment 2)
  2. Energy of isolated fragment 1 (with specified charge/multiplicity)
  3. Energy of isolated fragment 2 (charge derived from total - fragment 1)
  4. Formation energy of the combined system relative to isolated fragments
  5. Electron density difference: rho(combined) - rho(frag1) - rho(frag2)
"""

import os
import argparse
import numpy as np
from ase import Atoms
from ase.io import read as ase_read
from ase.calculators.nwchem import NWChem
from ase.units import Ha

# ---------------------------------------------------------------------------
# Command-line arguments
# ---------------------------------------------------------------------------
parser = argparse.ArgumentParser(
    description='DFT charge-transfer study between two fragments (ASE + NWChem).'
)
parser.add_argument(
    'frag1_xyz', type=str,
    help='Path to XYZ file for fragment 1'
)
parser.add_argument(
    'frag2_xyz', type=str,
    help='Path to XYZ file for fragment 2'
)
parser.add_argument(
    '--frag1-charge', type=int, required=True,
    help='Charge to constrain on fragment 1'
)
parser.add_argument(
    '--frag1-mult', type=int, required=True,
    help='Spin multiplicity of fragment 1'
)
parser.add_argument(
    '--total-charge', type=int, required=True,
    help='Total charge of the combined system'
)
parser.add_argument(
    '--total-mult', type=int, required=True,
    help='Total spin multiplicity of the combined system'
)
parser.add_argument(
    '--frag2-mult', type=int, default=None,
    help='Spin multiplicity of fragment 2 (default: same as fragment 1)'
)
parser.add_argument(
    '-np', '--nprocs', type=int, default=1,
    help='Number of MPI processors for NWChem (default: 1)'
)
parser.add_argument(
    '--basis', type=str, default='6-311++G**',
    help='Basis set (default: 6-311++G**)'
)
parser.add_argument(
    '--xc', type=str, default='B3LYP',
    help='Exchange-correlation functional (default: B3LYP)'
)
parser.add_argument(
    '--vacuum', type=float, default=6.0,
    help='Vacuum padding in angstroms (default: 6.0)'
)
parser.add_argument(
    '--scf-iterations', type=int, default=500,
    help='Maximum number of SCF iterations (default: 500)'
)
args = parser.parse_args()

# Derive fragment 2 charge and multiplicity
frag2_charge = args.total_charge - args.frag1_charge
frag2_mult = args.frag2_mult if args.frag2_mult is not None else args.frag1_mult

NWCHEM_COMMAND = f'mpirun -np {args.nprocs} nwchem PREFIX.nwi > PREFIX.nwo'
print(f"NWChem command: {NWCHEM_COMMAND}")

# ---------------------------------------------------------------------------
# 1. Read molecular geometries from XYZ files
# ---------------------------------------------------------------------------

frag1 = ase_read(args.frag1_xyz)
frag2 = ase_read(args.frag2_xyz)

# Build combined system by concatenating the two fragments
combined = frag1 + frag2

# Center the combined system with vacuum; then place isolated fragments
# in the same cell so that all three calculations share an identical grid
# for a consistent density difference.
combined.center(vacuum=args.vacuum)

# The centering operation shifts all atoms; recover the per-fragment
# positions from the (now-centered) combined system.
n1 = len(frag1)
frag1_isolated = Atoms(
    symbols=combined.symbols[:n1],
    positions=combined.positions[:n1],
    cell=combined.cell.copy(),
    pbc=combined.pbc,
)
frag2_isolated = Atoms(
    symbols=combined.symbols[n1:],
    positions=combined.positions[n1:],
    cell=combined.cell.copy(),
    pbc=combined.pbc,
)

print(f"\nFragment 1: {args.frag1_xyz} ({len(frag1)} atoms)")
print(f"  Charge: {args.frag1_charge}  Multiplicity: {args.frag1_mult}")
print(f"Fragment 2: {args.frag2_xyz} ({len(frag2)} atoms)")
print(f"  Charge: {frag2_charge}  Multiplicity: {frag2_mult}")
print(f"Combined system: {len(combined)} atoms")
print(f"  Total charge: {args.total_charge}  Multiplicity: {args.total_mult}")

# ---------------------------------------------------------------------------
# 2. Common NWChem DFT settings
# ---------------------------------------------------------------------------

basis = args.basis
xc = args.xc
grid_setting = 'fine'

common_dft = {
    'xc': xc,
    'grid': grid_setting,
    'convergence': 'energy 1e-7',
    'iterations': args.scf_iterations,
}

scratch_dir = os.path.abspath('nwchem_scratch')
os.makedirs(scratch_dir, exist_ok=True)


def cube_filename(label):
    return os.path.join(scratch_dir, f'{label}_density.cube')


def make_dplot_block(label, atoms):
    """Generate NWChem DPLOT block for writing electron density cube file."""
    cfile = cube_filename(label)
    cell = atoms.cell.array
    origin = np.zeros(3)
    end = cell.diagonal()
    return (
        f'  TITLE "{label} electron density"\n'
        f'  LimitXYZ\n'
        f'    {origin[0]:.4f} {end[0]:.4f} 100\n'
        f'    {origin[1]:.4f} {end[1]:.4f} 100\n'
        f'    {origin[2]:.4f} {end[2]:.4f} 100\n'
        f'  spin total\n'
        f'  gaussian\n'
        f'  output {cfile}\n'
    )


class NWChemWithDplot(NWChem):
    """NWChem calculator that appends a DPLOT block after the main task."""

    def __init__(self, dplot_extra='', **kwargs):
        self._dplot_extra = dplot_extra
        super().__init__(**kwargs)

    def write_input(self, atoms, properties=None, system_changes=None):
        super().write_input(atoms, properties, system_changes)
        if self._dplot_extra:
            infile = os.path.join(self.directory, self.input_filename())
            with open(infile, 'a') as f:
                f.write(self._dplot_extra)


def get_nwchem_calculator(label, charge, mult, atoms, extra_dft=None):
    """Create an NWChem calculator with the given settings."""
    dft_block = dict(common_dft)
    if mult > 1:
        dft_block['odft'] = None
        dft_block['mult'] = mult
    if extra_dft:
        dft_block.update(extra_dft)

    dplot_block = make_dplot_block(label, atoms)
    dplot_extra = f'\n\ndplot\n{dplot_block}end\n\ntask dplot\n'

    calc = NWChemWithDplot(
        dplot_extra=dplot_extra,
        label=os.path.join(scratch_dir, label),
        command=NWCHEM_COMMAND,
        basis=basis,
        charge=charge,
        task='energy',
        dft=dft_block,
    )
    return calc


# ---------------------------------------------------------------------------
# 3. Run calculations
# ---------------------------------------------------------------------------

print("\n" + "=" * 70)
print("DFT CALCULATIONS: CHARGE TRANSFER BETWEEN FRAGMENTS")
print(f"Method: {xc}/{basis}")
print("=" * 70)

# --- 3a. Combined system ---
print(f"\n>>> Calculation 1: Combined system")
print(f"    Total charge: {args.total_charge} | Multiplicity: {args.total_mult}")
combined.calc = get_nwchem_calculator(
    'combined', charge=args.total_charge, mult=args.total_mult, atoms=combined
)
E_combined = combined.get_potential_energy()
print(f"    Energy: {E_combined:.6f} eV  ({E_combined / Ha:.6f} Ha)")

# --- 3b. Isolated fragment 1 ---
print(f"\n>>> Calculation 2: Fragment 1 ({os.path.basename(args.frag1_xyz)})")
print(f"    Charge: {args.frag1_charge} | Multiplicity: {args.frag1_mult}")
frag1_isolated.calc = get_nwchem_calculator(
    'frag1', charge=args.frag1_charge, mult=args.frag1_mult, atoms=frag1_isolated
)
E_frag1 = frag1_isolated.get_potential_energy()
print(f"    Energy: {E_frag1:.6f} eV  ({E_frag1 / Ha:.6f} Ha)")

# --- 3c. Isolated fragment 2 ---
print(f"\n>>> Calculation 3: Fragment 2 ({os.path.basename(args.frag2_xyz)})")
print(f"    Charge: {frag2_charge} | Multiplicity: {frag2_mult}")
frag2_isolated.calc = get_nwchem_calculator(
    'frag2', charge=frag2_charge, mult=frag2_mult, atoms=frag2_isolated
)
E_frag2 = frag2_isolated.get_potential_energy()
print(f"    Energy: {E_frag2:.6f} eV  ({E_frag2 / Ha:.6f} Ha)")

# ---------------------------------------------------------------------------
# 4. Formation energy
# ---------------------------------------------------------------------------

E_formation = E_combined - (E_frag1 + E_frag2)

print("\n" + "=" * 70)
print("RESULTS")
print("=" * 70)
print(f"\nE(combined)   = {E_combined:.6f} eV")
print(f"E(fragment 1) = {E_frag1:.6f} eV")
print(f"E(fragment 2) = {E_frag2:.6f} eV")
print(f"\nFormation energy of combined system:")
print(f"  dE = E(combined) - E(frag1) - E(frag2)")
print(f"  dE = {E_formation:.6f} eV  ({E_formation * 23.0609:.4f} kcal/mol)")

# ---------------------------------------------------------------------------
# 5. Electron density difference from cube files
# ---------------------------------------------------------------------------

def read_cube(filename):
    """Read a Gaussian cube file. Returns (origin, axes, data, atoms_info)."""
    with open(filename, 'r') as f:
        f.readline()
        f.readline()
        parts = f.readline().split()
        natoms = int(parts[0])
        origin = np.array([float(x) for x in parts[1:4]])
        nx_line = f.readline().split()
        nx = int(nx_line[0])
        ax1 = np.array([float(x) for x in nx_line[1:4]])
        ny_line = f.readline().split()
        ny = int(ny_line[0])
        ax2 = np.array([float(x) for x in ny_line[1:4]])
        nz_line = f.readline().split()
        nz = int(nz_line[0])
        ax3 = np.array([float(x) for x in nz_line[1:4]])
        for _ in range(abs(natoms)):
            f.readline()
        data = []
        for line in f:
            data.extend(float(x) for x in line.split())
        data = np.array(data).reshape((nx, ny, nz))
    return origin, (ax1, ax2, ax3), (nx, ny, nz), data


def write_cube(filename, origin, axes, grid_shape, data, comment="density difference"):
    """Write a Gaussian cube file (no atom info, just volumetric data)."""
    nx, ny, nz = grid_shape
    ax1, ax2, ax3 = axes
    with open(filename, 'w') as f:
        f.write(f"{comment}\n")
        f.write("Generated by charge_transfer.py\n")
        f.write(f"    0 {origin[0]:12.6f} {origin[1]:12.6f} {origin[2]:12.6f}\n")
        f.write(f"{nx:5d} {ax1[0]:12.6f} {ax1[1]:12.6f} {ax1[2]:12.6f}\n")
        f.write(f"{ny:5d} {ax2[0]:12.6f} {ax2[1]:12.6f} {ax2[2]:12.6f}\n")
        f.write(f"{nz:5d} {ax3[0]:12.6f} {ax3[1]:12.6f} {ax3[2]:12.6f}\n")
        count = 0
        for ix in range(nx):
            for iy in range(ny):
                for iz in range(nz):
                    f.write(f" {data[ix, iy, iz]:13.5E}")
                    count += 1
                    if count % 6 == 0:
                        f.write("\n")
        if count % 6 != 0:
            f.write("\n")


print("\n--- Electron Density Difference ---")

cube_comb = cube_filename('combined')
cube_f1 = cube_filename('frag1')
cube_f2 = cube_filename('frag2')

if all(os.path.exists(f) for f in [cube_comb, cube_f1, cube_f2]):
    origin_c, axes_c, shape_c, rho_combined = read_cube(cube_comb)
    origin_1, axes_1, shape_1, rho_frag1 = read_cube(cube_f1)
    origin_2, axes_2, shape_2, rho_frag2 = read_cube(cube_f2)

    if shape_c == shape_1 == shape_2:
        delta_rho = rho_combined - rho_frag1 - rho_frag2

        voxel_vol = np.abs(np.dot(axes_c[0], np.cross(axes_c[1], axes_c[2])))

        total_excess = np.sum(delta_rho) * voxel_vol
        max_accum = np.max(delta_rho)
        min_deple = np.min(delta_rho)

        print(f"  Grid dimensions: {shape_c}")
        print(f"  Integrated delta_rho (should be ~0): {total_excess:.6f} e")
        print(f"  Max density accumulation:  {max_accum:.6e} e/bohr^3")
        print(f"  Max density depletion:     {min_deple:.6e} e/bohr^3")

        diff_cube = os.path.join(scratch_dir, 'density_difference.cube')
        write_cube(diff_cube, origin_c, axes_c, shape_c, delta_rho,
                   comment="Density difference: combined - frag1 - frag2")
        print(f"\n  Density difference cube file written to:\n    {diff_cube}")
        print("  (Visualize with VESTA, VMD, or similar)")
    else:
        print("  WARNING: Cube file grids do not match. Interpolation needed.")
        print(f"    Combined grid:   {shape_c}")
        print(f"    Fragment 1 grid: {shape_1}")
        print(f"    Fragment 2 grid: {shape_2}")
        print("  Skipping density difference. Use identical grid specs for all systems.")
else:
    missing = [f for f in [cube_comb, cube_f1, cube_f2] if not os.path.exists(f)]
    print("  Cube files not found (DPLOT may not have run):")
    for f in missing:
        print(f"    {f}")
    print("  Density difference analysis skipped.")

print("\n" + "=" * 70)
print("DONE")
print("=" * 70)
