#!/usr/bin/env python3

"""
================================================================================
EXTRACT PROTEIN SEQUENCE FROM PDB STRUCTURE
================================================================================

Purpose:
    This script extracts the amino acid sequence from a PDB (Protein Data Bank)
    structure file and exports it in FASTA format. This is useful for sequence
    alignment, database searches, and as input for homology modeling workflows.

Workflow:
    1. Parse PDB structure file using Biopython
    2. Iterate through all residues in the structure
    3. Convert 3-letter amino acid codes to 1-letter codes
    4. Write sequence to FASTA file for downstream analysis

Limitations:
    - This script extracts only standard 20 amino acids
    - Post-translational modifications (PTMs) are NOT included
    - Non-standard residues and heteroatoms are ignored
    - For detailed information on modified residues, see:
      * restyp.lib (standard residue definitions)
      * restyp_accelrys.lib (extended residue library)

Input:
    - PDB structure file (*.pdb)

Output:
    - FASTA sequence file (*.fasta)
    - Format: Header line (>pdb_code) followed by sequence line

Prerequisites:
    - Biopython (BioPython) package installed
    - Input PDB file in working directory or specified path

Usage:
    python3 sequence_from_pdb.py <input_pdb> <output_fasta>
    
    Arguments:
        input_pdb:     Path to input PDB structure file
        output_fasta:  Path to output FASTA sequence file
================================================================================
"""

# ===== IMPORTS =====
from Bio.PDB import PDBParser      # PDB structure parsing and manipulation
from Bio.PDB import PDBIO          # PDB file I/O operations (for potential future use)
import sys                          # Command-line argument handling
import os                           # File path and directory operations


# ===== AMINO ACID MAPPING =====
# Dictionary mapping 3-letter amino acid codes (PDB format) to 1-letter codes (FASTA format)
# Includes all 20 standard amino acids
AA = {
    "ALA": "A",  # Alanine
    "ARG": "R",  # Arginine
    "ASN": "N",  # Asparagine
    "ASP": "D",  # Aspartic acid
    "CYS": "C",  # Cysteine
    "GLN": "Q",  # Glutamine
    "GLU": "E",  # Glutamic acid
    "GLY": "G",  # Glycine
    "HIS": "H",  # Histidine
    "ILE": "I",  # Isoleucine
    "LEU": "L",  # Leucine
    "LYS": "K",  # Lysine
    "MET": "M",  # Methionine
    "PHE": "F",  # Phenylalanine
    "PRO": "P",  # Proline
    "SER": "S",  # Serine
    "THR": "T",  # Threonine
    "TRP": "W",  # Tryptophan
    "TYR": "Y",  # Tyrosine
    "VAL": "V"   # Valine
}


# ===== STRUCTURE PARSING =====
# Initialize PDB parser for reading and processing PDB structure files
parser = PDBParser()

# Parse command-line arguments
if len(sys.argv) != 3:
    print("Usage: python3 sequence_from_pdb.py <input_pdb> <output_fasta>")
    print("\nArguments:")
    print("  input_pdb:     Path to input PDB structure file")
    print("  output_fasta:  Path to output FASTA sequence file")
    sys.exit(1)

# Extract input and output file paths from command-line arguments
input_pdb = sys.argv[1]
output_fasta = sys.argv[2]

# Validate that input PDB file exists
if not os.path.isfile(input_pdb):
    print(f"Error: Input PDB file not found: {input_pdb}")
    sys.exit(1)

# Extract PDB code from input filename (remove path and extension)
pdb_code = os.path.splitext(os.path.basename(input_pdb))[0]

# Parse the PDB structure file and create a structure object
# Parameters:
#   - pdb_code: identifier for the structure (used for naming)
#   - input_pdb: path to the input PDB file (from command-line argument)
struct = parser.get_structure(pdb_code, input_pdb)


# ===== SEQUENCE EXTRACTION =====
# Initialize list to accumulate amino acid codes
seq = []

# Iterate through all residues in the structure
# Note: This includes all chains and models present in the PDB file
for r in struct.get_residues():
    # Extract residue name in 3-letter format (as stored in PDB)
    resname_3_letters = r.get_resname()
    
    # Convert to 1-letter code and append to sequence
    # Note: Non-standard residues will cause KeyError - consider adding error handling
    seq.append(AA[resname_3_letters])


# ===== FASTA OUTPUT =====
# Write extracted sequence to FASTA format file
# FASTA format: Header line (>) followed by sequence data
with open(output_fasta, 'w+') as fout:
    # Write FASTA header: ">" followed by sequence identifier
    fout.write(f'>{pdb_code}\n')
    
    # Write sequence: join all 1-letter codes into continuous string
    # This becomes the sequence record in FASTA format
    sequence_str = ''.join(seq)
    fout.write(f'{sequence_str}\n')


print(f"Sequence extraction complete. Output written to {output_fasta}")
print(f"Sequence length: {len(seq)} residues")