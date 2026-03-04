"""
Generate triangular graphene nanoflakes with zigzag edges.

Builds a graphene sheet via ASE, cuts a triangular region using a cross-product
point-in-triangle test, and validates that the resulting edges are pure zigzag
via sublattice analysis. Only sizes where n % 3 != 0 yield valid zigzag edges.
"""
import argparse
import warnings
import numpy as np
from ase.build import graphene


def generate_triangular_flake(n: int) -> "ase.Atoms":
    """
    Generate a triangular graphene nanoflake.

    Parameters
    ----------
    n : int
        Controls side length: side_length = n * 2.46 Å

    Returns
    -------
    ase.Atoms
        The triangular flake (no hydrogen passivation).
    """
    a = 2.46
    side_length = n * a
    sheet_cells = 3 * n + 4  # generous margin to avoid boundary artefacts

    ase_sheet = graphene(formula="C2", a=a, size=(sheet_cells, sheet_cells, 1))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        ase_sheet.center(vacuum=5, axis=2)

    xyz = ase_sheet.get_positions()
    # Use the true cell centre (not atom mean) so the triangle aligns correctly
    # regardless of the oblique cell geometry.
    cell = ase_sheet.get_cell()
    center = np.zeros(3)
    center[:2] = ((cell[0] + cell[1]) / 2)[:2]
    center[2] = np.mean(xyz[:, 2])
    h = (np.sqrt(3) / 2) * side_length

    v1 = center + np.array([0, 2 / 3 * h, 0])
    v2 = center + np.array([-side_length / 2, -1 / 3 * h, 0])
    v3 = center + np.array([side_length / 2, -1 / 3 * h, 0])

    def sign(p, a, b):
        return (p[:, 0] - b[0]) * (a[1] - b[1]) - (a[0] - b[0]) * (p[:, 1] - b[1])

    d1 = sign(xyz, v1, v2)
    d2 = sign(xyz, v2, v3)
    d3 = sign(xyz, v3, v1)
    has_neg = (d1 < 0) | (d2 < 0) | (d3 < 0)
    has_pos = (d1 > 0) | (d2 > 0) | (d3 > 0)
    inside = ~(has_neg & has_pos)

    flake = ase_sheet[inside]
    flake.set_pbc(False)
    flake.center(vacuum=5.0)

    return flake


def check_zigzag_edges(flake, bond_cutoff=1.6):
    """
    Validate that all edge atoms lie on the same sublattice.

    In a honeycomb lattice, zigzag edges have all edge atoms on one sublattice
    (A or B). If edge atoms appear on both sublattices, the edges are armchair
    or mixed, meaning the triangle size is not commensurate with the lattice.

    Returns
    -------
    is_zigzag : bool
    message : str
    """
    pos = flake.get_positions()
    n = len(flake)

    neighbors = [[] for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            if np.linalg.norm(pos[i] - pos[j]) < bond_cutoff:
                neighbors[i].append(j)
                neighbors[j].append(i)

    coord = np.array([len(nb) for nb in neighbors])

    # Assign sublattice A=0 / B=1 via BFS
    sublattice = np.full(n, -1, dtype=int)
    sublattice[0] = 0
    queue = [0]
    while queue:
        atom = queue.pop(0)
        for nb in neighbors[atom]:
            if sublattice[nb] == -1:
                sublattice[nb] = 1 - sublattice[atom]
                queue.append(nb)

    edge_atoms = np.where(coord == 2)[0]
    if len(edge_atoms) == 0:
        return False, "No 2-coordinated edge atoms found"

    sub = sublattice[edge_atoms]
    n_A = np.sum(sub == 0)
    n_B = np.sum(sub == 1)

    if n_A == 0 or n_B == 0:
        return True, f"All {len(edge_atoms)} edge atoms on sublattice {'A' if n_B==0 else 'B'} — zigzag edges confirmed"
    else:
        return False, f"Edge atoms split across both sublattices (A:{n_A}, B:{n_B}) — not a valid zigzag size"


def saturate(flake, bond_cutoff=1.6, ch_bond=1.09):
    """
    Passivate under-coordinated edge carbons with hydrogen atoms.

    sp2 carbon rules (120° bond angles, in-plane):
      coord=2 → 1 H, placed opposite the bisector of the two C-C bonds
      coord=1 → 2 H, placed at ±120° from the single C-C bond
    """
    from ase import Atoms

    pos = flake.get_positions()
    n = len(flake)

    neighbors = [[] for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            if np.linalg.norm(pos[i] - pos[j]) < bond_cutoff:
                neighbors[i].append(j)
                neighbors[j].append(i)

    h_pos = []
    for i in range(n):
        p = pos[i]
        nbs = neighbors[i]
        coord = len(nbs)
        if coord == 3:
            continue

        unit = [(pos[j] - p) / np.linalg.norm(pos[j] - p) for j in nbs]

        if coord == 2:
            # bisector direction, then negate
            h_dir = -(unit[0] + unit[1])
            h_dir /= np.linalg.norm(h_dir)
            h_pos.append(p + ch_bond * h_dir)

        elif coord == 1:
            # rotate the bond vector by ±120° in the xy plane
            v = unit[0]
            for sign in (+1, -1):
                angle = sign * 2 * np.pi / 3
                c, s = np.cos(angle), np.sin(angle)
                h_dir = np.array([c * v[0] - s * v[1],
                                  s * v[0] + c * v[1],
                                  0.0])
                h_pos.append(p + ch_bond * h_dir)

    if not h_pos:
        return flake
    result = flake.copy()
    result.extend(Atoms("H" * len(h_pos), positions=h_pos))
    return result


def truncate_corners(flake, bond_cutoff=1.6):
    """Remove 1-coordinated corner atoms, leaving 2-coordinated C as the new tips."""
    pos = flake.get_positions()
    n = len(flake)
    coord = np.zeros(n, dtype=int)
    for i in range(n):
        for j in range(i + 1, n):
            if np.linalg.norm(pos[i] - pos[j]) < bond_cutoff:
                coord[i] += 1
                coord[j] += 1
    return flake[coord != 1]


def main():
    parser = argparse.ArgumentParser(
        description="Generate triangular graphene nanoflakes with zigzag edges"
    )
    parser.add_argument(
        "-n", "--size", type=int, default=6,
        help="Size parameter: side_length = n * 2.46 Å (default: 6)",
    )
    parser.add_argument(
        "-v", "--visualize", action="store_true",
        help="Open structure in ASE GUI viewer",
    )
    parser.add_argument(
        "--output", type=str, default="triangular_flake.xyz",
        help="Output XYZ file (default: triangular_flake.xyz)",
    )
    parser.add_argument(
        "--truncate-corners", action="store_true",
        help="Remove 1-coordinated corner atoms, leaving 2-coordinated C as new tips",
    )
    parser.add_argument(
        "--saturate", action="store_true",
        help="Passivate under-coordinated edge carbons with hydrogen atoms",
    )
    args = parser.parse_args()

    flake = generate_triangular_flake(args.size)

    is_zigzag, msg = check_zigzag_edges(flake)
    print(f"Edge validation: {msg}")
    if not is_zigzag:
        print("WARNING: structure does not have pure zigzag edges — try a different size")

    if args.truncate_corners:
        flake = truncate_corners(flake)
        print("Truncated 3 corner atoms (1-coordinated C removed)")

    if args.saturate:
        n_c = len(flake)
        flake = saturate(flake)
        print(f"Saturated: added {len(flake) - n_c} H atoms")

    flake.write(args.output)
    print(f"Wrote {len(flake)} atoms to {args.output}")

    if args.visualize:
        from ase.visualize import view
        view(flake)


if __name__ == "__main__":
    main()
