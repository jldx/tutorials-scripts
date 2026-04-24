#!python3

################################################################################
#
# This script modifies a GROMACS .gro file to replace a specific water molecule with 
# a chloride ion (CLM), while ensuring the system remains electrically neutral. 
# It updates residue IDs, atom indices, and the total atom count accordingly.
#
# Key Steps:
#     1. File Parsing: Reads the input .gro file and splits it into sections (molecule, membrane, ions, water, box).
#     2. Chloride Creation: Replaces the coordinates of a selected water oxygen with a new chloride ion.
#     3. Water Removal: Deletes the selected water molecule (3 lines: O, H1, H2).
#     4. ID/Index Update: Adjusts residue IDs and atom indices for the remaining water molecules to maintain continuity.
#     
# Output: Writes the modified structure to a new .gro file.
#
# Assumptions:
#     - The input file follows the standard GROMACS .gro format.
#     - The water molecule to replace is specified by its index in the water section.
#     - The system must remain neutral, so one chloride is removed after adding the new one.
#
################################################################################


import numpy as np

# Define the new residue name for the chloride ion
new_resname = 'CLM'
# Input and output file names
file = 'conf21625.gro'
fout = 'test.gro'

# GRO file structure:
# line 0: header (Gromacs metadata)
# line 1: number of atoms
# lines 2-56398: molecule C8O
# lines 56399-100082: lipid POPC
# lines 100083-100173: sodium ion (SOD)
# lines 100174-100264: chloride ion (CLA)
# lines 100265-203305: water
# last line: box dimensions

# Read the input file
with open(file) as f:
    lines = f.readlines()

# Extract sections from the file
n_atoms = lines[1]  # Number of atoms
lines_mol = lines[2:56398]  # Molecule C8O
lines_memb = lines[56398:100082]  # Lipid POPC
lines_sod = lines[100082:100173]  # Sodium ion (SOD)
lines_cla = lines[100173:100264]  # Chloride ion (CLA)
lines_wat = lines[100264:-1]  # Water
box = lines[-1]  # Box dimensions

# --- Create a new chloride ion with the coordinates of the water oxygen ---
# Index of the water oxygen to replace (in the water section)
water_to_replace = 49629
# Extract the water line and its coordinates
water_line = lines_wat[water_to_replace]
water_coords = water_line[20:].rstrip()  # Coordinates of the water oxygen

# Extract residue and index from the first chloride ion
resid_new_cla = int(lines_cla[0][0:5])
index_new_cla = int(lines_cla[0][15:20])

# Create a new chloride line with the water coordinates
new_cla_line = f"{resid_new_cla:>5}{new_resname}{new_resname:>7}{index_new_cla:5d}{water_coords}\n"
# Add the new chloride and remove the last one to maintain neutrality
new_lines_cla = [new_cla_line] + lines_cla[1:-1]

# --- Delete the water and update water residue IDs and atom indices ---
# Remove the water molecule (3 lines: O, H1, H2)
del lines_wat[water_to_replace:(water_to_replace + 3)]
n_wat = len(lines_wat) // 3  # Number of remaining water molecules

# Calculate new residue IDs for water
first_resid = int(new_lines_cla[-1][0:5]) + 1
last_resid = first_resid + n_wat + 1
new_wat_resids = np.repeat(np.arange(first_resid, last_resid), 3)

# Calculate new atom indices for water
first_index = int(new_lines_cla[-1][15:20]) + 1
last_index = int(lines_wat[-1][15:20]) - 4  # Adjust for deleted water and chloride
new_wat_indices = np.hstack((np.arange(first_index, 99999 + 1), np.arange(last_index + 1)))

# Rebuild water lines with updated IDs and indices
new_lines_wat = []
for k, line in enumerate(lines_wat):
    resid = new_wat_resids[k]
    index = new_wat_indices[k]
    resname = line[5:10]
    atomname = line[10:15]
    coords = line[20:].rstrip()
    new_lines_wat.append(f"{resid:>5}{resname}{atomname}{index:>5}{coords}\n")

# Update the total number of atoms
new_n_atoms = f'{int(n_atoms) - 4}\n'  # Subtract 4 atoms (1 water + 1 chloride)

# --- Merge all sections and write the output file ---
new_lines = [
    'Updated gro\n',
    new_n_atoms,
    *lines_mol,
    *lines_memb,
    *lines_sod,
    *new_lines_cla,
    *new_lines_wat,
    box
]

# Write the new GRO file
with open(fout, 'w') as out:
    out.writelines(new_lines)
