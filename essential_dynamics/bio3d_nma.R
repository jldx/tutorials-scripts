#!/usr/bin/env Rscript

################################################################################
# NORMAL MODE ANALYSIS (NMA) USING BIO3D
################################################################################
#
# Purpose:
#   This script performs Normal Mode Analysis (NMA) on protein structures using
#   the bio3d R package. NMA predicts protein motions and collective dynamics by
#   analyzing the harmonic oscillations of the protein structure around its
#   equilibrium conformation.
#
# Workflow:
#   1. Load protein structure from PDB file
#   2. Select specific chains (if multi-chain protein)
#   3. Filter to C-alpha atoms only (reduces computational complexity)
#   4. Perform NMA with multiple force fields:
#      - calpha: Simple elastic network model based on C-alpha distances
#      - anm: Anisotropic Network Model (standard NMA)
#      - pfanm: Parameter-free ANM
#      - reach: Refined elastic network contact-based model
#      - sdenm: Smooth domain elastic network model
#   5. Generate trajectory files for modes 7-10 (lowest frequency motions)
#
# Output:
#   - PDB trajectory files for each force field and mode
#   - Format: {forcefield}_modes_{mode_number}.pdb
#   - Example: calpha_modes_7.pdb, anm_modes_8.pdb, etc.
#
# Prerequisites:
#   - R with bio3d package installed
#   - Input PDB structure file
#   - Single or multi-chain protein structure
#
# Force Fields Explained:
#   - calpha: Simplest model; treats protein as network of C-alpha nodes
#   - anm: Standard approach; considers anisotropic (directional) interactions
#   - pfanm: Automatically optimizes spring constant parameters
#   - reach: Accounts for network connectivity and spatial reach
#   - sdenm: Models domain dynamics with smooth transitions
#
# Author: Essential Dynamics Tutorial
# Last Modified: 2026-04-27
#
################################################################################

# ===== LOAD REQUIRED LIBRARY =====
# bio3d: Bioinformatics tools for structural analysis and dynamics
library(bio3d)


# ===== READ PDB STRUCTURE =====
# Load protein structure from PDB file
# Input file should be a valid PDB structure (download from RCSB if needed)
pdb <- read.pdb("<input_file>")

# Display basic structure information
cat("Structure loaded:\n")
cat("  Total atoms:", nrow(pdb$atom), "\n")
cat("  Unique chains:", paste(unique(pdb$atom$chain), collapse=", "), "\n")
cat("  Residues:", nrow(pdb$seqres), "\n\n")


# ===== CHAIN SELECTION (If Multi-Chain) =====
# Select specific chains for analysis (e.g., chains A and B from a dimer)
# This step is optional but recommended for multi-chain complexes
cat("Selecting chains for analysis...\n")

sel_chains <- atom.select(pdb, 
                          chain = c('A', 'B'),      # Chains to include
                          verbose = TRUE)


# ===== C-ALPHA ATOM SELECTION =====
# Filter structure to include ONLY C-alpha atoms
# This significantly reduces computational complexity while retaining key dynamics
# C-alpha atoms are standard proxies for backbone dynamics in NMA
cat("Filtering to C-alpha atoms only...\n")

sel_calpha <- atom.select(pdb,
                          elety = "CA",             # Elemental type: C-alpha
                          chain = c('A', 'B'),      # Apply to selected chains
                          verbose = TRUE)

# Apply chain and C-alpha selection
pdb_trim <- trim.pdb(pdb, sel_calpha)

cat("  Retained atoms:", nrow(pdb_trim$atom), "\n")
cat("  Retained residues:", nrow(pdb_trim$seqres), "\n\n")


# ===== NORMAL MODE ANALYSIS WITH MULTIPLE FORCE FIELDS =====
# Perform NMA using different force field parameterizations
# Each force field has different strengths for capturing specific types of motion
cat("Starting Normal Mode Analysis with multiple force fields...\n\n")

# Iterate through force field types
for (f in c('calpha', 'anm', 'pfanm', 'reach', 'sdenm')) {
  cat("Processing force field:", f, "\n")
  
  # Perform NMA on C-alpha filtered structure
  # Parameters:
  #   - pdb_trim: Trimmed structure (C-alpha only)
  #   - ff: Force field type
  #   - temp: Temperature in Kelvin (affects thermal scaling)
  #   - mass: Use atomic masses for weighting (generally improves results)
  modes <- nma(pdb_trim,
               ff = f,
               temp = 300,                   # Room temperature (physiological relevance)
               mass = TRUE)                  # Include mass weighting in calculations
  
  # Generate trajectory files for individual modes
  # Modes 7-10 are typically low-frequency collective motions (biologically relevant)
  # Modes 1-6 are rigid body motions (translation/rotation)
  for (m in 7:10) {
    # Create trajectory PDB file for this mode
    # Format: {forcefield}_modes_{mode_number}.pdb
    output_file <- paste(paste(f, 'modes', m, sep="_"), 'pdb', sep = '.')
    
    cat("  Generating trajectory for mode", m, "->", output_file, "\n")
    
    # Generate trajectory with interpolation
    # Creates a multi-frame PDB showing motion along this normal mode
    mktrj.nma(modes,
              mode = m,                      # Mode number to visualize
              file = output_file)            # Output trajectory filename
  }
  
  cat("  Force field", f, "complete.\n\n")
}

cat("NMA analysis complete. Trajectory files generated.\n")
cat("Next steps: Visualize trajectories in PyMOL or similar software.\n")