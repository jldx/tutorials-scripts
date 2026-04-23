################################################################################
# PyMOL Backbone Spine Visualization Plugin
################################################################################
#
# DESCRIPTION:
#   Custom PyMOL script to create smooth, high-quality backbone spline 
#   visualizations. This script generates a continuous tube representation
#   of the protein backbone using B-spline interpolation on CA atoms.
#
# Was used for representing oligourea helices during the postdoc at LBT
#
# PURPOSE:
#   - Create smooth cartoon-like backbone representation
#   - Higher visual quality than default trace representation
#   - Useful for presentations and publications
#   - Customizable tube radius and interpolation density
#   - Flexible coloring options
#
# USAGE:
#   1. Load this script in PyMOL:
#      run backbone_spine.py
#
#   2. After loading PDB structure, call:
#      bb_spine()                    # Uses default parameters
#      bb_spine(sel='chain A and name CA', radius=0.5, color=(1,0,0), name='my_spine')
#      bb_spine(sel='name CA', steps_per_residue=500)  # High-quality spline
#
# DEPENDENCIES:
#   - PyMOL (with cmd and cgo modules)
#   - NumPy (numerical operations)
#   - SciPy (splprep/splev for B-spline interpolation)
#   - Python 3.x
#
# OUTPUT:
#   - CGO (Compiled Graphics Object) representation in PyMOL
#   - Default object name: spiner{radius}_s{steps_per_residue}
#     (Example: spiner0.3_s200)
#
################################################################################

from pymol import cmd, cgo
import numpy as np
from scipy.interpolate import splprep, splev
import string


def bb_spine(sel='name CA', radius=0.3, steps_per_residue=200, color=(0.5, 0.5, 0.5), name=None):
    """
    Draw a super-smooth CGO backbone spline using B-spline interpolation.
    
    PARAMETERS:
    -----------
    sel : str
        PyMOL selection string for atoms to trace
        Default: 'name CA' (alpha carbon backbone trace)
        Examples:
          - 'name CA': all CA atoms (standard backbone)
          - 'chain A and name CA': CA atoms in chain A
          - 'name CA + name C': CA and C atoms for finer backbone
          - 'resi 10-50 and name CA': specific residue range
    
    radius : float
        Tube thickness in Ångströms (default: 0.3 Å)
        Visual effects:
          - 0.1-0.2: thin wireframe effect
          - 0.3-0.5: standard representation
          - 0.7-1.0: thick ribbon-like representation
          - 1.5+: blocky tube effect
    
    steps_per_residue : int
        Interpolation density - number of spline segments per residue
        (default: 200, which gives smooth curves)
        Visual quality vs. performance:
          - 50-100: Fast, reasonable smoothness
          - 200: Default (good balance)
          - 300-500: High quality (slower rendering)
          - 1000+: Very smooth (slow on large proteins)
    
    color : tuple of 3 floats
        RGB color (range 0.0-1.0)
        Default: (0.5, 0.5, 0.5) - gray
        Examples:
          - (1.0, 0.0, 0.0): red
          - (0.0, 0.0, 1.0): blue
          - (0.2, 0.8, 0.2): green
          - (1.0, 0.5, 0.0): orange
    
    name : str or None
        CGO object name in PyMOL (default: None)
        If None, automatically generates name from parameters
        Generated name format: spiner{radius}_s{steps_per_residue}
        Example: spiner0.3_s200
        Custom name: 'my_backbone_spine'
    
    RETURNS:
    --------
    None (creates CGO object visible in PyMOL viewport)
    
    EXAMPLE USAGE:
    ---------------
    # Standard backbone trace
    bb_spine()
    
    # Thick red backbone for chain A
    bb_spine(sel='chain A and name CA', radius=0.5, color=(1,0,0), name='chain_A_spine')
    
    # High-quality blue backbone
    bb_spine(steps_per_residue=500, color=(0,0,1), name='hires_spine')
    
    # Multiple spines with different colors for different chains
    bb_spine(sel='chain A and name CA', color=(1,0,0), name='chainA')
    bb_spine(sel='chain B and name CA', color=(0,0,1), name='chainB')
    """
    
    # ========================================================================
    # INPUT VALIDATION AND CONVERSION
    # ========================================================================
    # Convert string arguments to proper types (PyMOL sometimes passes strings)
    radius = float(radius)
    steps_per_residue = int(steps_per_residue)

    # ========================================================================
    # OBJECT NAMING
    # ========================================================================
    # Generate default name if not provided
    # This allows running the function multiple times without name conflicts
    if not name:
        name = f'spiner{radius}_s{steps_per_residue}'

    # ========================================================================
    # EXTRACT COORDINATES
    # ========================================================================
    # Get PyMOL model object from selection
    # model.atom contains all atoms in the selection
    model = cmd.get_model(sel)
    
    # Extract 3D coordinates as NumPy array
    # Shape: (N_atoms, 3) where each row is [x, y, z]
    coords = np.array([a.coord for a in model.atom])

    # ========================================================================
    # B-SPLINE INTERPOLATION
    # ========================================================================
    # Create smooth spline through CA atom coordinates
    # splprep: Parametric spline representation
    #   coords.T: Transpose to get (3, N_atoms) format (x, y, z coordinates)
    #   s=0: No smoothing (exact interpolation through all points)
    #   k=3: Cubic spline (smooth and reasonably fast)
    # Returns:
    #   tck: Spline coefficients (knots, coefficients, degree)
    #   u: Parameter values for input points
    tck, _ = splprep(coords.T, s=0, k=3)
    
    # Generate high-density interpolated points
    # Calculate total number of interpolation points
    # More points = smoother spline but slower rendering
    n_points = len(coords) * steps_per_residue
    
    # Create evenly-spaced parameter values from 0 to 1
    # These are used to sample the spline at high density
    u_new = np.linspace(0, 1, int(n_points))
    
    # Evaluate spline at high-density points
    # splev returns separate x, y, z coordinate arrays
    x_new, y_new, z_new = splev(u_new, tck)

    # ========================================================================
    # BUILD CGO OBJECT (Compiled Graphics Object)
    # ========================================================================
    # CGO is PyMOL's low-level graphics format
    # It's much faster than other representations for custom geometry
    obj = []
    
    # Loop through consecutive pairs of interpolated points
    # Each pair becomes a cylinder in the CGO
    for i in range(len(x_new) - 1):
        obj.extend([
            # Set color for this cylinder
            cgo.COLOR, *color,
            
            # Draw cylinder between point i and i+1
            cgo.CYLINDER,
            # Start point coordinates
            x_new[i], y_new[i], z_new[i],
            # End point coordinates
            x_new[i+1], y_new[i+1], z_new[i+1],
            # Tube radius (in Ångströms)
            radius,
            # Color at start and end of cylinder (both same as above)
            *color, *color
        ])

    # ========================================================================
    # LOAD AND DISPLAY CGO OBJECT
    # ========================================================================
    # Convert CGO list to PyMOL displayable object
    # Creates new object in viewport with given name
    cmd.load_cgo(obj, name)


# ============================================================================
# REGISTER COMMAND WITH PYMOL
# ============================================================================
# Make bb_spine function available as a PyMOL command
# After loading this script, users can call: bb_spine(...)
cmd.extend("bb_spine", bb_spine)

# ============================================================================
# EXAMPLE USAGE (commented out)
# ============================================================================
# Uncomment lines below to automatically generate spines when script loads

# # Define color (RGB values 0.0-1.0)
# obj_col_code = (0.486, 0.529, 0.557)  # Steel blue-gray
#
# # Select backbone atoms (CA atoms is standard, or add C atoms for finer backbone)
# cmd.do(f'select backbon_spline, (name CA + name C)')
#
# # Create smooth backbone spine
# bb_spine(sel=f'backbon_spline', color=obj_col_code, name=f'spline')

