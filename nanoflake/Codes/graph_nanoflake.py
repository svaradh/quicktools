"""
Unified graphene nanoflake generator — zigzag edges only.

Dispatches to the appropriate backend based on shape:
  hexagonal  →  hexagonal.py   (hexagon-tiling, guaranteed zigzag by construction)
  triangular →  triangular.py  (lattice-cut, zigzag validated via sublattice check)

Only zigzag edges are produced. For armchair or alternating edges see the
legacy implementation in git history.
"""
import argparse
import sys
import os

sys.path.insert(0, os.path.dirname(__file__))

from hexagonal import generate_hexagonal_flake, saturate as saturate_hex
from triangular import (
    generate_triangular_flake,
    check_zigzag_edges,
    truncate_corners,
    saturate as saturate_tri,
)


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Generate graphene nanoflakes with zigzag edges. "
            "Use -s to select hexagonal (tiling) or triangular (lattice-cut) shape."
        )
    )
    parser.add_argument(
        "-s", "--shape", choices=["hexagonal", "triangular"], default="hexagonal",
        help="Flake shape (default: hexagonal)",
    )
    parser.add_argument(
        "-n", "--size", type=int, default=3,
        help=(
            "Size parameter. "
            "Hexagonal: number of rings around the central hexagon "
            "(n=0→benzene, n=1→coronene, n=2→circumcoronene; default: 3). "
            "Triangular: side length in units of 2.46 Å "
            "(valid zigzag sizes: n not divisible by 3; default: 3)."
        ),
    )
    parser.add_argument(
        "--saturate", action="store_true",
        help="Passivate under-coordinated edge carbons with hydrogen atoms",
    )
    parser.add_argument(
        "--truncate-corners", action="store_true",
        help="Triangular only: remove 1-coordinated corner C atoms",
    )
    parser.add_argument(
        "--output", type=str, default=None,
        help="Output XYZ file (default: hexagonal_flake.xyz or triangular_flake.xyz)",
    )
    parser.add_argument(
        "-v", "--visualize", action="store_true",
        help="Open structure in ASE GUI viewer",
    )
    args = parser.parse_args()

    if args.shape == "hexagonal":
        if args.truncate_corners:
            parser.error("--truncate-corners is only valid for triangular flakes")

        flake = generate_hexagonal_flake(args.size)
        print(f"Generated hexagonal flake: {len(flake)} C atoms (n={args.size})")

        if args.saturate:
            n_c = len(flake)
            flake = saturate_hex(flake)
            print(f"Saturated: added {len(flake) - n_c} H atoms")

        output = args.output or "hexagonal_flake.xyz"

    else:  # triangular
        flake = generate_triangular_flake(args.size)
        print(f"Generated triangular flake: {len(flake)} C atoms (n={args.size})")

        is_zigzag, msg = check_zigzag_edges(flake)
        print(f"Edge validation: {msg}")
        if not is_zigzag:
            print("WARNING: not a valid zigzag size — try n not divisible by 3")

        if args.truncate_corners:
            before = len(flake)
            flake = truncate_corners(flake)
            print(f"Truncated {before - len(flake)} corner atoms")

        if args.saturate:
            n_c = len(flake)
            flake = saturate_tri(flake)
            print(f"Saturated: added {len(flake) - n_c} H atoms")

        output = args.output or "triangular_flake.xyz"

    flake.write(output)
    print(f"Wrote {len(flake)} atoms to {output}")

    if args.visualize:
        from ase.visualize import view
        view(flake)


if __name__ == "__main__":
    main()
