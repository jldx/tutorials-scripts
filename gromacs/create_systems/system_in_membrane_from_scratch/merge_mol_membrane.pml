# PyMOL Script: Merge Molecule and Membrane
# Description: This script merges a molecule and a membrane, removing overlapping lipids to avoid clashes.

# Retain the original atom order to avoid issues with topology
set retain_order, 1

# Load the centered molecule and membrane
load mol_center.gro, mol
load membrane_center.pdb, membrane

# Select and remove membrane lipids within 4 Å of the molecule to avoid clashes
select neighbors, byres (membrane within 4 of mol)
remove neighbors

# Optionally, select and remove POPC lipids that may still be clashing (e.g., inside pores)

# Save the merged structure
save mol_membrane.pdb
