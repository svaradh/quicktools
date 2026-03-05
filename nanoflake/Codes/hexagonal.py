"""
Generate hexagonal graphene nanoflakes with zigzag edges by tiling hexagons.

Builds the flake by placing flat-top benzene-ring hexagons outward from a
central hexagon using axial (cube) coordinates. Adjacent hexagon centres are
connected at the correct lattice positions so shared vertices are exact in
floating-point, and deduplication via a rounded-coordinate dict is reliable.
"""
import argparse
import numpy as np
from ase import Atoms


def generate_hexagonal_flake(n: int) -> "ase.Atoms":
    """
    Generate a hexagonal graphene nanoflake.

    Parameters
    ----------
    n : int
        Number of hexagon rings around the central hexagon.
        n=0 →  6 C (benzene)
        n=1 → 24 C (coronene)
        n=2 → 54 C (circumcoronene)

    Returns
    -------
    ase.Atoms
        Hexagonal graphene flake, no hydrogen passivation, pbc=False.
    """
    a = 1.42  # C-C bond length = hexagon edge length, Å

    # Primitive vectors for the triangular lattice of hexagon centres.
    # For flat-top hexagons, edge-sharing neighbours are at 30°, 90°, 150°, …
    # so the two independent primitive directions are 30° and 90°:
    #   v1 = a√3 · (cos30°, sin30°) = (3a/2,  a√3/2)
    #   v2 = a√3 · (cos90°, sin90°) = (0,     a√3  )
    v1 = np.array([1.5 * a,           a * np.sqrt(3) / 2.0])
    v2 = np.array([0.0,               a * np.sqrt(3)])

    # Collect unique carbon positions.
    # Axial cube-coordinate condition |i|≤n, |j|≤n, |i+j|≤n selects the
    # hexagonal patch of N rings.
    atom_dict = {}
    for i in range(-n, n + 1):
        for j in range(-n, n + 1):
            if abs(i + j) > n:
                continue
            cx, cy = i * v1 + j * v2
            for k in range(6):
                angle = k * np.pi / 3.0
                x = cx + a * np.cos(angle)
                y = cy + a * np.sin(angle)
                key = (round(x, 6), round(y, 6))
                atom_dict[key] = (x, y)

    positions_3d = np.array([[x, y, 0.0] for x, y in atom_dict.values()])
    flake = Atoms("C" * len(positions_3d), positions=positions_3d)
    flake.set_pbc(False)
    flake.center(vacuum=5.0)
    return flake


def saturate(flake, bond_cutoff=1.6, ch_bond=1.09):
    """
    Passivate under-coordinated edge carbons with hydrogen atoms.

    sp2 carbon rules (120° bond angles, in-plane):
      coord=2 → 1 H, placed opposite the bisector of the two C-C bonds
      coord=1 → 2 H, placed at ±120° from the single C-C bond
    """
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
            h_dir = -(unit[0] + unit[1])
            h_dir /= np.linalg.norm(h_dir)
            h_pos.append(p + ch_bond * h_dir)

        elif coord == 1:
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


def main():
    parser = argparse.ArgumentParser(
        description="Generate hexagonal graphene nanoflakes"
    )
    parser.add_argument(
        "-n", "--rings", type=int, default=3,
        help="Number of hexagon rings around the central hexagon (default: 3). "
             "n=0 → benzene, n=1 → coronene, n=2 → circumcoronene",
    )
    parser.add_argument(
        "--saturate", action="store_true",
        help="Passivate under-coordinated edge carbons with hydrogen atoms",
    )
    parser.add_argument(
        "--output", type=str, default="hexagonal_flake.xyz",
        help="Output XYZ file (default: hexagonal_flake.xyz)",
    )
    parser.add_argument(
        "-v", "--visualize", action="store_true",
        help="Open structure in ASE GUI viewer",
    )
    args = parser.parse_args()

    flake = generate_hexagonal_flake(args.rings)
    print(f"Generated {len(flake)} C atoms (n={args.rings})")

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
