"""
extract_coords.py
-----------------
Extract atomic coordinates from a molecular dynamics trajectory using MDAnalysis.

For each frame (or a slice thereof), the script selects atoms matching a given
MDAnalysis selection string, collects their residue names, residue IDs, atom names,
and positions along the requested spatial axes, and stacks everything into a 3D array
of shape (n_frames, n_atoms, n_cols) where n_cols = metadata columns + coordinate columns.

The result is saved as a compressed NumPy archive (.npz).

Output array column layout (example for coords_type='xz'):
    col 0 : resname  (str)
    col 1 : resid    (int)
    col 2 : atom name (str)
    col 3 : x position (float)
    col 4 : z position (float)

Usage
-----
    python extract_coords.py \\
        --top   ../top.pdb \\
        --traj  ../traj.xtc \\
        --sel   "protein" \\
        --start 0 \\
        --end   100 \\
        --step  1 \\
        --coords xz \\
        --out   coords.npz

    # Minimal usage (all frames, xyz, default selection):
    python extract_coords.py --top top.pdb --traj traj.xtc

Dependencies
------------
    MDAnalysis >= 2.0
    NumPy >= 1.20
"""

import argparse
import sys

import MDAnalysis as mda
import numpy as np


# ---------------------------------------------------------------------------
# Mapping from shorthand axis label to column indices in positions array
# ---------------------------------------------------------------------------
COORDS_INDICES = {
    'x':   0,
    'y':   1,
    'z':   2,
    'xy':  [0, 1],
    'xz':  [0, 2],
    'yz':  [1, 2],
    'xyz': [0, 1, 2],
}


def parse_args():
    parser = argparse.ArgumentParser(
        description="Extract MD trajectory coordinates into a compressed .npz array.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--top",    required=True,          help="Topology file path (e.g. .pdb, .gro, .psf)")
    parser.add_argument("--traj",   required=True,          help="Trajectory file path (e.g. .xtc, .dcd, .trr)")
    parser.add_argument("--sel",    default="protein",     help="MDAnalysis atom selection string")
    parser.add_argument("--start",  type=int,   default=0,  help="First frame index (inclusive)")
    parser.add_argument("--end",    type=int,   default=None, help="Last frame index (exclusive). Default: all frames")
    parser.add_argument("--step",   type=int,   default=1,  help="Step between frames")
    parser.add_argument("--coords", default="xyz",          choices=list(COORDS_INDICES.keys()),
                        help="Spatial axes to extract")
    parser.add_argument("--out",    default="coords.npz",   help="Output .npz file path")
    return parser.parse_args()


def extract(top_path, traj_path, selection_str, start_frame, end_frame, step, coords_type):
    """
    Load universe, iterate over the requested frame slice, and return a 3D array.

    Returns
    -------
    coords : np.ndarray, shape (n_frames, n_atoms, n_cols)
        Mixed-type object array containing metadata + positions per atom per frame.
    """
    u = mda.Universe(top_path, traj_path, in_memory=False)
    sel = u.select_atoms(selection_str)

    if sel.n_atoms == 0:
        sys.exit(f"[ERROR] Selection '{selection_str}' matched 0 atoms. Check your selection string.")

    print(f"[INFO] Topology  : {top_path}")
    print(f"[INFO] Trajectory: {traj_path}")
    print(f"[INFO] Selection : '{selection_str}' → {sel.n_atoms} atoms")
    print(f"[INFO] Frames    : [{start_frame} : {end_frame} : {step}]")
    print(f"[INFO] Axes      : {coords_type}")

    col_idx = COORDS_INDICES[coords_type]  # int or list of ints
    frames = []

    for ts in u.trajectory[start_frame:end_frame:step]:
        # Metadata columns — reshaped to (n_atoms, 1) for horizontal stacking
        resname  = sel.atoms.resnames.reshape(-1, 1)          # e.g. [['LIP'], ['LIP'], ...]
        resids   = sel.atoms.resids.reshape(-1, 1)             # e.g. [[1], [2], ...]
        names    = sel.atoms.names.reshape(-1, 1)              # e.g. [['CA'], ['N'], ...]

        # Coordinate columns — .copy() is critical: positions is a live view into
        # MDAnalysis internal buffer and will be overwritten on the next frame
        positions = sel.positions[:, col_idx].copy()           # (n_atoms,) or (n_atoms, n_axes)
        if positions.ndim == 1:
            positions = positions.reshape(-1, 1)               # ensure 2D for concatenation

        # Stack metadata + positions horizontally → (n_atoms, 3 + n_axes)
        frame_data = np.concatenate([resname, resids, names, positions], axis=1)
        frames.append(frame_data)

    # Stack all frames → (n_frames, n_atoms, n_cols)
    coords = np.stack(frames, axis=0)
    print(f"[INFO] Output shape: {coords.shape}  (frames × atoms × cols)")
    return coords


def main():
    args = parse_args()

    coords = extract(
        top_path      = args.top,
        traj_path     = args.traj,
        selection_str = args.sel,
        start_frame   = args.start,
        end_frame     = args.end,
        step          = args.step,
        coords_type   = args.coords,
    )

    # Save as compressed archive; load with np.load(args.out)['coords']
    np.savez_compressed(args.out, coords=coords)
    print(f"[INFO] Saved → {args.out}")


if __name__ == "__main__":
    main()