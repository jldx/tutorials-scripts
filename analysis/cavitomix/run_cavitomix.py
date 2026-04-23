#!/usr/bin/env python3
"""
run_cavitomix.py (Suggested explicit name)
---------------------------------------------------
CavitoMix: Automated Protein Cavity Detection (Headless Mode)

This is a **command-line automation script** designed to run cavity detection
without opening the PyMOL GUI. It is meant for:
  - Batch processing multiple protein structures
  - Integration into automated pipelines and workflows
  - Running on headless/remote systems (no display needed)
  - Scripting and parallel execution
  - High-throughput screening of protein cavities

CavitoMix combines multiple cavity detection methods (geometric and energetic)
to provide robust identification of potential drug binding sites.

Workflow (Headless / Automated):
  1. Initialize CavFind object with analysis parameters
  2. Load protein structure from PDB file
  3. Run cavity detection algorithms (no GUI)
  4. Export cavity results to PDB format
  5. (Optional) Process results programmatically or via PyMOL in batch mode

Usage (Command-line):
    python3 detect_cavities_batch.py
    python3 detect_cavities_batch.py --input protein.pdb --output cavities.pdb

Suggested rename to: `detect_cavities_batch.py`
  - Clearly indicates batch/headless operation
  - Emphasizes no PyMOL GUI involvement
  - Better describes automated workflow intent
  - Distinguishes from interactive PyMOL plugin usage

Requirements:
  - CavitoMix package (cavitomix module in current directory)
  - PDB structure file of the protein
  - (Optional) Custom settings/parameters for cavity detection

Output:
  - PDB file with detected cavities as pseudo-atoms or spheres
  - Cavity statistics and analysis data (ready for post-processing)
  - No PyMOL GUI interaction required
"""

from cavitomix import *
from cavitomix.cavfind import CavFind

# ==============================================================================
# USER CONFIGURATION
# ==============================================================================

# Protein structure file (PDB format)
# Replace with your actual protein PDB file path
protein_pdb_file = "<path_to_protein_structure.pdb>"

# Output file for detected cavities (PDB format)
output_cavities_file = "<output_cavities_with_results.pdb>"

# Optional: Custom CavitoMix settings (leave as None to use defaults)
# You can customize parameters like:
#   - Probe radius for cavity detection
#   - Minimum cavity volume threshold
#   - Clustering parameters
#   - etc.
custom_settings = None

# ==============================================================================
# CAVITY DETECTION WORKFLOW
# ==============================================================================

# Step 1: Initialize CavFind object
# ----
# Creates a CavFind instance with cavity detection parameters.
# Parameters:
#   - name: identifier for this analysis run (for logging)
#   - obj_name: name of the molecular object being analyzed
#   - settings: custom parameter dictionary (None = use defaults from cavfind.py)
# Default settings include:
#   - grid_spacing: 0.7 Å (granularity for 3D grid)
#   - probe_radius: 1.4 Å (solvent probe size, typically water)
#   - ligsite_cutoff: 5 (threshold for LigSite algorithm)
#   - keep_hydrogens: False, keep_hetatms: False (clean protein atoms only)
# Implementation: cavfind.py CavFind.__init__() and default_settings dict

cav = CavFind(
    name="protein_analysis",           # Analysis identifier
    obj_name="protein_structure",       # Object name
    settings=custom_settings            # Custom settings (None = defaults)
)

print("✓ CavFind object initialized")

# Step 2: Load protein structure from PDB file
# ----
# Parses PDB format and filters atoms based on settings:
#   - Removes hydrogens (unless keep_hydrogens=True)
#   - Removes heteroatoms (unless keep_hetatms=True)
#   - Removes water molecules (unless keep_waters=True)
#   - Handles alternate conformations (selects conformation 'A' by default)
# Then extracts atomic properties and stores as numpy arrays:
#   - coords: 3D coordinates of retained atoms (n_atoms × 3)
#   - radii: van der Waals radii (depends on UA/AA setting)
#   - hp: hydrophobicity parameters (from radii.py charges/radii tables)
#   - charges: partial atomic charges
# Implementation: cavfind.py struct_from_pdb() and pdb_structure.py PDBStructure class

with open(protein_pdb_file, "r") as pdb_file:
    cav.struct_from_pdb(pdb_file)

print(f"✓ Protein structure loaded from: {protein_pdb_file}")

# Step 3: Run cavity detection algorithms
# ----
# Executes the complete LigSite cavity detection pipeline:
#   A. Grid Setup (setup_grid from ligsite.py):
#      - Creates 3D grid with spacing=0.7 Å around protein
#      - Grid origin and extent calculated from min/max atomic coordinates
#      - Includes cushion parameter for padding around structure
#   B. Grid Masking (mask_grid from ligsite.py):
#      - Marks grid points occupied by protein atoms (PROTEIN_FLAG = -100)
#      - Marks "soft shell" around atoms for probe interaction (SOFT_FLAG = -99)
#      - Uses atom radii × radius_factor + probe_radius for exclusion zones
#   C. LigSite Algorithm (do_ligsite from ligsite.py):
#      - Analyzes grid to identify cavity-like regions
#      - Computes distance field from protein atoms
#      - Identifies connected components of grid points
#   D. Cavity Detection (find_cavities from ligsite.py):
#      - Clusters grid points into distinct cavities
#      - Filters by ligsite_cutoff threshold and gap parameter
#      - Computes cavity volume and geometry
#   E. Cavity Annotation (if annotate=True):
#      - Calculates Coulomb potential for each cavity
#      - Computes hydrophobicity descriptors
#      - Identifies nearby residues and binding pocket properties
# Implementation: cavfind.py run() method calling ligsite.py functions

cav.run()

print("✓ Cavity detection completed")

# Step 4: Export detected cavities to PDB file
# ----
# Converts cavity objects to PDB format with pseudo-atoms representing
# cavity centers and properties:
#   - Each cavity represented by one or more pseudo-atoms
#   - Coordinates: cavity center/grid points
#   - B-factor field: cavity properties (volume, score, etc.)
#   - Residue name: cavity identifier
# Output can be visualized directly in PyMOL, Chimera, or VMD
# Allows overlaying cavities with protein structure for analysis
# Implementation: cavfind.py write_cavities() method

cav.write_cavities(filename=output_cavities_file)

print(f"✓ Cavity results written to: {output_cavities_file}")

# ==============================================================================
# POST-ANALYSIS (Optional)
# ==============================================================================

# You can further analyze or filter cavities based on computed properties:
# - cav.cavities: list of detected cavity objects
# - cav.get_cavity_properties(): extract volume, shape, residues, etc.
# - cav.filter_by_volume(min_volume, max_volume): volume-based filtering
# - etc.

print("\n✅ CavitoMix analysis completed successfully!")
print(f"\nNext steps:")
print(f"  1. Visualize cavities in PyMOL using the local plugin:")
print(f"     - Launch PyMOL from this directory (cavitomix/)")
print(f"     - The PyMOL plugin will be automatically available")
print(f"     - Load protein: {protein_pdb_file}")
print(f"     - Load cavities: {output_cavities_file}")
print(f"     - Use plugin features for interactive cavity analysis")
print(f"  2. Examine cavity properties and select promising drug targets")
print(f"  3. Perform molecular docking in selected cavities")
                
