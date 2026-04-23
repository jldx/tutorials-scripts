#!python3

# This script reorders the atoms in a POPC membrane PDB file to match the expected 
# order defined in the topology. This is necessary to ensure compatibility with 
# GROMACS topology files and avoid issues during simulation setup.

# Usage: python3 reorder_atoms_POPC.py input.pdb output.pdb

import sys

# Define the correct order of atom names as per the topology for POPC (no hydrogens)
correct_order = [
    "N", "C13", "C14", "C15", "C12", "C11", "P", "O13", "O14", "O11", "O12", "C1",
    "C2", "O21", "C21", "O22", "C22", "C3", "O31", "C31", "O32", "C32", "C23",
    "C24", "C25", "C26", "C27", "C28", "C29", "C210", "C211", "C212", "C213",
    "C214", "C215", "C216", "C217", "C218", "C33", "C34", "C35", "C36", "C37",
    "C38", "C39", "C310", "C311", "C312", "C313", "C314", "C315", "C316"
]

# Read the original PDB file, skipping the first two and last two lines (header/footer)
with open(sys.argv[1], "r") as f:
    lines = f.readlines()[2:-2]

n_atoms = len(correct_order)
# Split the lines into residues, assuming each residue has `n_atoms` atoms
residues = [lines[i:(i + n_atoms)] for i in range(0, len(lines), n_atoms)]

atom_num = 1

# Write the reordered atoms to the output file
with open(sys.argv[2], "w") as f:
    for res in residues:
        atom_dict = {}
        # Create a dictionary of atoms for the current residue
        for atom in res:
            atom_name = atom.split()[2]
            atom_dict[atom_name] = atom

        # Write atoms in the correct order
        for atom in correct_order:
            f.write(f"ATOM{atom_num:>7} {atom_dict[atom][12:]}")
            atom_num += 1

        # Add TER record to indicate the end of the residue
        f.write("TER\n")
