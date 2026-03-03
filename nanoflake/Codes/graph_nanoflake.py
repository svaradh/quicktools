from ase import Atoms
import numpy as np

def create_graphene_flake(n, shape='hexagonal', edge_type='zigzag', saturated=True,
                          vacuum=5.0, orientation='xy', visualize=False):
    """
    Create a graphene nanoflake with specified shape and edge type.

    n: size parameter (atoms along each edge)
    shape: 'hexagonal' or 'triangular'
    edge_type: 'zigzag', 'armchair', or 'alternating' (hexagonal only)
    saturated: if True, passivate edges with H atoms
    vacuum: vacuum spacing in Angstroms
    orientation: plane orientation - 'xy', 'xz', or 'yz'
    visualize: if True, open structure in ASE GUI viewer
    """
    a_cc = 1.42  # C-C bond length

    # Honeycomb lattice vectors
    a1 = np.array([np.sqrt(3) * a_cc, 0])
    a2 = np.array([np.sqrt(3) * a_cc / 2, 1.5 * a_cc])

    # Two atoms per unit cell
    basis = [np.array([0, 0]), np.array([0, a_cc])]

    # Generate large supercell
    size = 2 * n + 4
    all_positions = []
    for i in range(-size, size + 1):
        for j in range(-size, size + 1):
            cell_origin = i * a1 + j * a2
            for b in basis:
                all_positions.append(cell_origin + b)

    all_positions = np.array(all_positions)

    # Center positions
    center = np.mean(all_positions, axis=0)
    all_positions -= center

    # Build flake based on shape and edge type
    if shape == 'triangular':
        # Direct construction for triangular flakes
        unique_positions = _build_triangular_flake(n, a_cc, edge_type)
    else:
        # Use masking for hexagonal flakes
        if edge_type == 'zigzag':
            mask = _hexagon_mask_zigzag(all_positions, n, a_cc)
        elif edge_type == 'armchair':
            mask = _hexagon_mask_armchair(all_positions, n, a_cc)
        elif edge_type == 'alternating':
            mask = _hexagon_mask_alternating(all_positions, n, a_cc)
        else:
            raise ValueError(f"Unknown edge_type: {edge_type}")

        positions = all_positions[mask]
        unique_positions = _remove_duplicates(positions)

    # Create 3D positions
    positions_3d = np.array([[p[0], p[1], 0] for p in unique_positions])

    flake = Atoms('C' * len(positions_3d), positions=positions_3d)
    flake.center(vacuum=vacuum)

    if saturated:
        flake = _saturate_edges(flake)

    # Rotate to desired orientation
    if orientation == 'xz':
        positions = flake.positions.copy()
        flake.positions[:, 1] = -positions[:, 2]
        flake.positions[:, 2] = positions[:, 1]
        flake.center(vacuum=vacuum)
    elif orientation == 'yz':
        positions = flake.positions.copy()
        flake.positions[:, 0] = -positions[:, 2]
        flake.positions[:, 2] = positions[:, 0]
        flake.center(vacuum=vacuum)

    if visualize:
        from ase.visualize import view
        view(flake)

    return flake


def _saturate_edges(flake):
    """Add hydrogen atoms to passivate edge carbons."""
    from ase import Atom

    cutoff = 1.5  # Å
    h_positions = []

    for i in range(len(flake)):
        pos_i = flake.positions[i]
        # Find neighbors
        neighbors = []
        for j in range(len(flake)):
            if i != j:
                d = np.linalg.norm(pos_i - flake.positions[j])
                if d < cutoff:
                    neighbors.append(flake.positions[j])

        n_bonds = len(neighbors)
        if n_bonds < 3:
            # Calculate existing bond angles
            existing_angles = []
            for neighbor in neighbors:
                diff = neighbor[:2] - pos_i[:2]
                angle = np.arctan2(diff[1], diff[0])
                existing_angles.append(angle)

            # Add H atoms in missing directions (~120° apart)
            for _ in range(3 - n_bonds):
                best_angle = _find_best_angle(existing_angles)
                if best_angle is not None:
                    direction = np.array([np.cos(best_angle), np.sin(best_angle), 0])
                    h_positions.append(pos_i + 1.09 * direction)
                    existing_angles.append(best_angle)

    if h_positions:
        h_atoms = Atoms('H' * len(h_positions), positions=h_positions)
        flake = flake + h_atoms

    return flake


def _find_best_angle(existing_angles):
    """Find angle furthest from existing bond angles."""
    if not existing_angles:
        return 0.0
    best_angle = None
    best_min_dist = -1
    for test_angle in np.linspace(-np.pi, np.pi, 360):
        min_dist = min(abs(np.arctan2(np.sin(test_angle - ea), np.cos(test_angle - ea)))
                      for ea in existing_angles)
        if min_dist > best_min_dist:
            best_min_dist = min_dist
            best_angle = test_angle
    return best_angle


def _hexagon_mask_zigzag(positions, n, a_cc):
    """Zigzag-edged hexagon mask."""
    # Zigzag edge length
    R = n * np.sqrt(3) * a_cc
    mask = []
    for p in positions:
        x, y = abs(p[0]), abs(p[1])
        # Hexagon with flat top/bottom (zigzag edges)
        in_hex = (y <= R * np.sqrt(3) / 2) and (y <= np.sqrt(3) * (R - x))
        mask.append(in_hex and x <= R)
    return np.array(mask)


def _hexagon_mask_armchair(positions, n, a_cc):
    """Armchair-edged hexagon mask."""
    # Armchair: pointy top (rotated 30 degrees from zigzag)
    R = n * 1.5 * a_cc
    mask = []
    for p in positions:
        x, y = abs(p[0]), abs(p[1])
        # Hexagon with pointy top (armchair edges)
        in_hex = (x <= R * np.sqrt(3) / 2) and (x <= np.sqrt(3) * (R - y))
        mask.append(in_hex and y <= R)
    return np.array(mask)


def _hexagon_mask_alternating(positions, n, a_cc):
    """Alternating zigzag/armchair edges (dodecagonal-like)."""
    # Use 12-sided polygon approximation
    R = n * a_cc * 1.5
    mask = []
    for p in positions:
        r = np.sqrt(p[0]**2 + p[1]**2)
        if r < 0.1:
            mask.append(True)
            continue
        theta = np.arctan2(p[1], p[0])
        # 12-fold symmetry boundary
        theta_mod = abs(theta) % (np.pi / 6)
        R_boundary = R / np.cos(theta_mod - np.pi / 12)
        mask.append(r <= R_boundary * 1.05)
    return np.array(mask)


def _build_triangular_flake(n, a_cc, edge_type):
    """Build triangular graphene flake with proper edge termination.

    For zigzag: all 3 edges are zigzag, edge atoms on same sublattice
    For armchair: all 3 edges are armchair
    n = number of edge atoms per side
    """
    positions = []

    if edge_type == 'zigzag':
        # Zigzag triangle: edge atoms form sublattice A in triangular arrangement
        # B atoms fill interior where they have 3 A neighbors

        # Sublattice A: triangular lattice with n atoms per edge
        # Use axial coordinates for triangular lattice
        a_tri = a_cc * np.sqrt(3)  # triangular lattice constant

        # Basis vectors for triangular lattice of A atoms
        v1 = np.array([a_tri, 0])
        v2 = np.array([a_tri / 2, a_tri * np.sqrt(3) / 2])

        # A atoms: indices (i, j) where i >= 0, j >= 0, i + j <= n - 1
        A_positions = {}
        for i in range(n):
            for j in range(n - i):
                pos = i * v1 + j * v2
                A_positions[(i, j)] = pos
                positions.append(pos.tolist())

        # B atoms: placed at center of "up" triangles formed by 3 A atoms
        # A B atom at (i, j) connects A(i,j), A(i+1,j), A(i,j+1)
        # Only add if all 3 A neighbors exist
        b_offset = (v1 + v2) / 3  # center of triangle

        for i in range(n - 1):
            for j in range(n - 1 - i):
                # Check all 3 A neighbors exist
                if (i, j) in A_positions and (i+1, j) in A_positions and (i, j+1) in A_positions:
                    b_pos = A_positions[(i, j)] + b_offset
                    positions.append(b_pos.tolist())

    elif edge_type == 'armchair':
        # Armchair triangle: edges have alternating A-B atoms
        # Build using dimer rows along armchair direction

        # For armchair, we build row by row where each row is an armchair line
        # Row spacing (perpendicular to armchair direction)
        row_spacing = a_cc * 3 / 2

        # Within-row atom spacing
        atom_spacing = a_cc * np.sqrt(3)

        # Vertical offset between A and B in same dimer
        dimer_offset = a_cc / 2

        for row in range(n):
            n_dimers = n - row  # number of A-B pairs in this row

            # Row y position
            y_base = row * row_spacing * np.sqrt(3) / 2
            x_offset = row * row_spacing / 2

            for i in range(n_dimers):
                x = x_offset + i * atom_spacing

                # A atom (lower of dimer)
                positions.append([x, y_base])

                # B atom (upper of dimer, offset in x and y)
                if row < n - 1 or i < n_dimers - 1:  # Skip B at top vertex
                    positions.append([x + atom_spacing / 2, y_base + a_cc])

    else:
        raise ValueError(f"Unknown edge_type: {edge_type}")

    return positions


def _remove_duplicates(positions, tol=0.1):
    """Remove duplicate positions within tolerance."""
    unique = []
    for p in positions:
        is_dup = False
        for u in unique:
            if np.linalg.norm(p - u) < tol:
                is_dup = True
                break
        if not is_dup:
            unique.append(p)
    return unique


# Keep old function name for backwards compatibility
def create_wulff_graphene_flake(n, saturated=True, vacuum=5.0, orientation='xy', visualize=False):
    """Backwards compatible wrapper - creates hexagonal flake with zigzag edges."""
    return create_graphene_flake(n, shape='hexagonal', edge_type='zigzag',
                                  saturated=saturated, vacuum=vacuum,
                                  orientation=orientation, visualize=visualize)

    if saturated:
        from ase import Atom

        # Identify under-coordinated C atoms and add H along missing bond directions
        cutoff = 1.5  # Å
        h_positions = []

        # Expected bond directions in sp2 graphene (120° apart)
        bond_angles = [0, 2*np.pi/3, 4*np.pi/3]

        for i in range(len(flake)):
            pos_i = flake.positions[i]
            # Find neighbors
            neighbors = []
            for j in range(len(flake)):
                if i != j:
                    d = np.linalg.norm(pos_i - flake.positions[j])
                    if d < cutoff:
                        neighbors.append(flake.positions[j])

            n_bonds = len(neighbors)
            if n_bonds < 3:
                # Find missing bond directions
                if n_bonds == 0:
                    # Isolated atom - add 3 H atoms
                    for angle in bond_angles:
                        direction = np.array([np.cos(angle), np.sin(angle), 0])
                        h_positions.append(pos_i + 1.09 * direction)
                else:
                    # Calculate existing bond angles
                    existing_angles = []
                    for neighbor in neighbors:
                        diff = neighbor[:2] - pos_i[:2]
                        angle = np.arctan2(diff[1], diff[0])
                        existing_angles.append(angle)

                    # Find missing directions (should be ~120° from existing bonds)
                    for _ in range(3 - n_bonds):
                        # Find angle furthest from all existing angles
                        best_angle = None
                        best_min_dist = -1
                        for test_angle in np.linspace(-np.pi, np.pi, 360):
                            min_dist = min(abs(np.arctan2(np.sin(test_angle - ea), np.cos(test_angle - ea)))
                                          for ea in existing_angles)
                            if min_dist > best_min_dist:
                                best_min_dist = min_dist
                                best_angle = test_angle
                        if best_angle is not None:
                            direction = np.array([np.cos(best_angle), np.sin(best_angle), 0])
                            h_positions.append(pos_i + 1.09 * direction)
                            existing_angles.append(best_angle)

        if h_positions:
            h_atoms = Atoms('H' * len(h_positions), positions=h_positions)
            flake += h_atoms

    # Rotate flake to desired orientation
    if orientation == 'xz':
        # Rotate 90° around x-axis: y -> z, z -> -y
        positions = flake.positions.copy()
        flake.positions[:, 1] = -positions[:, 2]
        flake.positions[:, 2] = positions[:, 1]
        flake.center(vacuum=vacuum)
    elif orientation == 'yz':
        # Rotate 90° around y-axis: x -> z, z -> -x
        positions = flake.positions.copy()
        flake.positions[:, 0] = -positions[:, 2]
        flake.positions[:, 2] = positions[:, 0]
        flake.center(vacuum=vacuum)
    # 'xy' is default, no rotation needed

    if visualize:
        from ase.visualize import view
        view(flake)

    return flake

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='Generate graphene nanoflakes')
    parser.add_argument('-n', '--size', type=int, default=4,
                        help='Size parameter - atoms along each edge (default: 4)')
    parser.add_argument('-s', '--shape', choices=['hexagonal', 'triangular'], default='hexagonal',
                        help='Flake shape (default: hexagonal)')
    parser.add_argument('-e', '--edge-type', choices=['zigzag', 'armchair', 'alternating'],
                        default='zigzag', help='Edge type (default: zigzag). '
                        'Note: alternating only available for hexagonal shape')
    parser.add_argument('--no-saturate', action='store_true',
                        help='Disable hydrogen saturation of edges')
    parser.add_argument('-o', '--orientation', choices=['xy', 'xz', 'yz'], default='xy',
                        help='Flake plane orientation (default: xy)')
    parser.add_argument('-v', '--visualize', action='store_true',
                        help='Open structure in ASE GUI viewer')
    parser.add_argument('--output', type=str, default='graphene_flake.xyz',
                        help='Output filename (default: graphene_flake.xyz)')

    args = parser.parse_args()

    # Validate edge type for triangular shape
    if args.shape == 'triangular' and args.edge_type == 'alternating':
        parser.error("Triangular flakes don't support 'alternating' edge type")

    flake = create_graphene_flake(
        n=args.size,
        shape=args.shape,
        edge_type=args.edge_type,
        saturated=not args.no_saturate,
        orientation=args.orientation,
        visualize=args.visualize
    )
    flake.write(args.output)
    print(f'Wrote {len(flake)} atoms to {args.output}')   
