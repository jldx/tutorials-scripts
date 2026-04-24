################################################################################
# Utilities for Pore Structure Cleaning and Atom Mapping
################################################################################
#
# This module provides functions for processing crystallographic PDB structures
# to prepare them for molecular dynamics simulations with force fields.
#
# FUNCTIONS:
#   - reorder_pdb_by_atom_number(): Reorder PDB atoms by force field indices
#   - map_and_remap_pdb_fragments(): Map crystal indices to force field indices
#
################################################################################

import os
import sys
import string
import glob
import pandas as pd
import MDAnalysis as mda

def split_molecular_fragments(input_pdb, out_folder):
    """
    Split a multi-molecule PDB structure into individual molecular fragments.
    
    PARAMETERS:
    -----------
    input_pdb : str
        Path to input PDB file containing one or more molecules.
        File must be readable by MDAnalysis.
        Examples:
          - 'structure.pdb'
          - '/path/to/crystal_structure.pdb'
          - 'complex.pdb' (containing multiple chains/molecules)
    
    out_folder : str
        Path to output directory where individual molecule PDB files
        will be saved. Directory is created if it doesn't exist.
        
    RETURNS:
    --------
    int : Number of molecular fragments found and written
    
    RAISES:
    -------
    FileNotFoundError: If input_pdb file does not exist
    ValueError: If input_pdb cannot be parsed by MDAnalysis
    
    WORKFLOW:
    ---------
    1. Create output directory if needed
    2. Load structure from PDB file
    3. Perform bond guessing to identify connected atoms
    4. Identify all molecular fragments (connected components)
    5. Write each fragment to separate file (tmp_0.pdb, tmp_1.pdb, ...)
    
    NOTES:
    ------
    - Fragment detection uses covalent bonding information
    - MDAnalysis.guess_bonds() estimates bond connectivity from distances
    - Output files are named tmp_0.pdb, tmp_1.pdb, etc.
    - Numbering starts from 0 and continues sequentially
    
    EXAMPLE:
    --------
    n_molecules = split_molecular_fragments('crystal.pdb', 'molecules/')
    print(f"Found and split {n_molecules} molecules")
    """
    
    # ========================================================================
    # VALIDATE INPUT FILE
    # ========================================================================
    if not os.path.exists(input_pdb):
        raise FileNotFoundError(f"Input PDB file not found: {input_pdb}")
    
    # ========================================================================
    # CREATE OUTPUT DIRECTORY
    # ========================================================================
    if not os.path.exists(out_folder):
        os.makedirs(out_folder)
        print(f"Created output directory: {out_folder}")
    
    # ========================================================================
    # LOAD MOLECULAR STRUCTURE
    # ========================================================================
    # MDAnalysis Universe: Main container for molecular data
    # Both input_pdb arguments: first for coordinates, second for topology
    # (when identical, loads same file for both)
    try:
        mol = mda.Universe(input_pdb, input_pdb)
        print(f"Loaded structure from: {input_pdb}")
        print(f"Total atoms in structure: {len(mol.atoms)}")
    except Exception as e:
        raise ValueError(f"Failed to parse PDB file: {e}")
    
    # ========================================================================
    # BOND DETECTION AND FRAGMENT IDENTIFICATION
    # ========================================================================
    # guess_bonds(): Estimate covalent bonds based on interatomic distances
    # Uses van der Waals radii and distance cutoffs to detect bonding
    # Creates connectivity graph of the structure
    mol.atoms.guess_bonds()
    
    # fragments: Identifies connected components in the molecular graph
    # Each fragment is a set of atoms that are connected via bonds
    # Returns list of AtomGroup objects (one per molecular fragment)
    n_mols = len(mol.atoms.fragments)
    print(f"Identified {n_mols} molecular fragment(s)")
    
    # ========================================================================
    # SAVE INDIVIDUAL MOLECULAR FRAGMENTS
    # ========================================================================
    # Iterate through all detected molecular fragments
    for m in range(n_mols):
        # Get the mth fragment (subset of atoms forming one molecule)
        fragment = mol.atoms.fragments[m]
        
        # Define output filename with sequential numbering
        output_file = os.path.join(out_folder, f'tmp_{m}.pdb')
        
        # Write fragment to PDB file
        fragment.write(output_file)
        
        # Print progress information
        n_atoms = len(fragment)
        print(f"  Fragment {m}: {n_atoms} atoms → {output_file}")
    
    return n_mols
    
def read_atom_map(map_file):
    """
    Load and prepare atom mapping dataframe from H5_atoms.map file.
    
    The H5_atoms.map file contains the mapping between crystallographic
    atom indices and force field atom indices/names. This function reads,
    validates, and organizes the mapping for efficient lookup during
    atom remapping operations.
    
    INPUT FILE FORMAT:
    ------------------
    Tab-separated values (TSV) with columns:
        xindex    index    name    resname    resid
        1         1        C1      OUA        1
        2         2        C2      OUA        1
        ...       ...      ...     ...        ...
    
    - xindex (int): Original crystal atom index (1-based)
                    Used as lookup key for crystal PDB files
    - index (int): Target force field atom index
                   Will replace xindex in remapped PDB
    - name (str): Force field atom name (4 chars max)
                  Examples: 'C1', 'N1', 'O1', etc.
    - resname (str): Force field residue name (3 chars max)
                     Example: 'OUA' (oligourea)
    - resid (int): Force field residue ID within molecule
    
    PARAMETERS:
    -----------
    map_file : str
        Path to atom mapping file (e.g., 'H5_atoms.map')
        Should be tab or space-separated
        Lines starting with '#' are treated as comments (ignored)
    
    RETURNS:
    --------
    pandas.DataFrame : Prepared atom mapping table with:
        - Index: xindex values (crystal atom indices)
        - Sorted in ascending xindex order (1, 2, 3, ...)
        - Reset index for efficient pandas.loc[] lookups
        - Ready for use in atom remapping operations
    
    PROCESSING STEPS:
    -----------------
    1. Read TSV file using pandas.read_table()
       - sep='\\s+': Handle both tabs and spaces as delimiters
       - comment='#': Ignore lines starting with '#'
    2. Sort by 'xindex' column (ascending order)
       - Ensures sequential access and easier verification
    3. Reset index (drop old index)
       - Allows direct xindex-based lookup with .loc[]
    
    EXAMPLE:
    --------
    # Load atom mapping
    amap = read_atom_map('H5_atoms.map')
    
    # Look up mapping for atom 1
    mapping_1 = amap.loc[1]
    print(mapping_1['index'])    # → 1 (new index)
    print(mapping_1['name'])     # → 'C1' (atom name)
    print(mapping_1['resname'])  # → 'OUA' (residue name)
    
    # Look up mapping for atom 42
    mapping_42 = amap.loc[42]
    # Used during PDB line remapping in map_and_remap_pdb_fragments()
    
    USAGE IN WORKFLOW:
    ------------------
    The returned dataframe is passed to:
        map_and_remap_pdb_fragments(..., amap=amap, ...)
    
    Inside that function, atom mappings are looked up via:
        xatom_mapped = amap.loc[xindex]
        new_index = xatom_mapped['index']
        new_name = xatom_mapped['name']
    """
    
    return (
        pd.read_table(map_file, sep='\\s+', comment='#')    # Read TSV with flexible delimiters
        .sort_values('xindex')                               # Sort by crystal atom index
        .reset_index(drop=True)                              # Reset index for .loc[] access
    )


def reorder_pdb_by_atom_number(pdb_content):
    """
    Reorder PDB file entries to match force field atom numbering sequence.
    
    This function reorganizes PDB records to ensure atoms appear in the order
    specified by the force field (H5_atoms.map), with ANISOU records properly
    paired with their corresponding HETATM records.
    
    CONTEXT:
    --------
    This function is a critical part of the atom mapping workflow. It's called
    by map_and_remap_pdb_fragments() after atoms have been remapped with new
    indices, to ensure the final PDB file has atoms in the correct force field
    order.
    
    WORKFLOW INTEGRATION:
    ---------------------
    1. map_and_remap_pdb_fragments() reads fitted PDB files
    2. Updates atom indices and names (xindex → force field index)
    3. Calls reorder_pdb_by_atom_number() to sort atoms
    4. Writes final remapped PDB with atoms in force field sequence
    
    INPUT:
    ------
    pdb_content : str
        Complete PDB file content as a single string containing:
        - CRYST1: Crystallography record (should be first)
        - HETATM: Atom records (will be sorted)
        - ANISOU: Anisotropic temperature factor records (paired with HETATM)
        - Other records: TER, END, REMARK, etc. (moved to end)
    
    OUTPUT:
    -------
    str : Reordered PDB content with structure:
        1. CRYST1 record(s) first
        2. HETATM records sorted by atom number
        3. ANISOU records paired with their corresponding HETATM
        4. All other records at the end
    
    PROCESS:
    --------
    The function performs the following operations:
    
    1. Parse input into record types
       - CRYST1 records: Crystallography box parameters
       - HETATM records: Heteroatom (non-protein) coordinates
       - ANISOU records: Anisotropic B-factors for thermal motion
       - Other: TER, END, REMARK, etc.
    
    2. Extract atom serial numbers from HETATM and ANISOU records
       - PDB format: columns 7-11 contain the atom serial number
       - Must handle fixed-width format parsing
    
    3. Sort HETATM and ANISOU records by atom serial number
       - Ensures atoms appear in ascending numerical order
       - Order is critical for force field assignment
    
    4. Match ANISOU records to HETATM records
       - Creates dictionary of ANISOU records keyed by atom number
       - Each HETATM is immediately followed by its ANISOU (if present)
       - Maintains crystallographic data association with atoms
    
    5. Reconstruct PDB output
       - CRYST1 records first (box parameters)
       - For each sorted HETATM: output HETATM then matching ANISOU
       - Other records appended at end
    
    EXAMPLE:
    --------
    # Read a PDB file
    with open('fragment_sc.pdb', 'r') as f:
        pdb_content = f.read()
    
    # Reorder by force field sequence
    pdb_sorted = reorder_pdb_by_atom_number(pdb_content)
    
    # Write reordered PDB
    with open('fragment_mapped.pdb', 'w') as f:
        f.write(pdb_sorted)
    
    TECHNICAL NOTES:
    ----------------
    - HETATM atoms are heteroatoms (not standard protein atoms)
    - ANISOU records provide anisotropic temperature factors
    - PDB format uses fixed-width columns (not tab/space delimited)
    - Atom number is 1-based in PDB files
    - Sorting ensures reproducible force field assignment
    - ANISOU data links thermal motion to specific atoms
    """
    
    # ========================================================================
    # PARSE PDB CONTENT INTO RECORD TYPES
    # ========================================================================
    # Split content into individual lines
    lines = pdb_content.strip().split('\n')
    
    # Initialize lists for each record type
    cryst_lines = []      # Crystallography records (box parameters)
    hetatm_lines = []     # Heteroatom coordinate records
    anisou_lines = []     # Anisotropic B-factor records
    other_lines = []      # TER, END, REMARK, etc.
    
    # Categorize each line by record type
    for line in lines:
        if line.startswith('CRYST1'):
            cryst_lines.append(line)
        elif line.startswith('HETATM'):
            hetatm_lines.append(line)
        elif line.startswith('ANISOU'):
            anisou_lines.append(line)
        else:
            other_lines.append(line)
    
    # ========================================================================
    # HELPER FUNCTION: Extract atom number from PDB record
    # ========================================================================
    def get_atom_number(line):
        """
        Extract atom serial number from PDB HETATM/ANISOU record.
        
        PDB format specification:
        - Columns 7-11 (0-indexed: 6-11): Atom serial number (right-aligned)
        - Format: 5-character field containing integer
        
        Example line:
        "HETATM    1  C1  OUA A   1      10.000  20.000  30.000  1.00  0.00           C"
         123456789...
         ^^^^^^^   ← columns 7-11 (this line: " 1" → 1)
        
        Args:
            line (str): Single PDB record line
            
        Returns:
            int: Atom serial number (1-based)
        """
        return int(line[6:11].strip())
    
    # ========================================================================
    # SORT RECORDS BY ATOM NUMBER
    # ========================================================================
    # Sort HETATM records by atom serial number
    # This ensures atoms appear in ascending order (1, 2, 3, ...)
    hetatm_lines.sort(key=get_atom_number)
    
    # Sort ANISOU records by atom serial number
    # Must match HETATM order for proper pairing
    anisou_lines.sort(key=get_atom_number)
    
    # ========================================================================
    # CREATE ANISOU LOOKUP DICTIONARY
    # ========================================================================
    # Build dictionary: atom_number → ANISOU_record
    # Allows O(1) lookup of ANISOU record for any atom
    anisou_dict = {}
    for anisou_line in anisou_lines:
        atom_num = get_atom_number(anisou_line)
        anisou_dict[atom_num] = anisou_line
    
    # ========================================================================
    # RECONSTRUCT REORDERED PDB OUTPUT
    # ========================================================================
    # Initialize result list to store reordered records in proper sequence
    result = []
    
    # Step 1: Add all CRYST1 records first
    # These contain crystallography box parameters needed before atoms
    result.extend(cryst_lines)
    
    # Step 2: Add HETATM records with paired ANISOU records
    # Loop through sorted HETATM records in order
    for hetatm_line in hetatm_lines:
        atom_num = get_atom_number(hetatm_line)
        
        # Add the HETATM record
        result.append(hetatm_line)
        
        # Add corresponding ANISOU record if it exists
        # Not all atoms may have ANISOU data
        if atom_num in anisou_dict:
            result.append(anisou_dict[atom_num])
    
    # Step 3: Add any remaining records (TER, END, etc.) at the end
    result.extend(other_lines)
    
    # ========================================================================
    # RETURN REORDERED CONTENT
    # ========================================================================
    # Join all lines with newline separator
    return '\n'.join(result)


def map_and_remap_pdb_fragments(fragments, pore, amap, seg='POR1', reorder_fn=reorder_pdb_by_atom_number):
    """
    Map atom indices and remap PDB files for all fragments in a pore.
    
    This function performs the following operations for each fragment:
    1. Reads the fitted PDB file (with complete sidechains)
    2. Maps crystal atom indices to force field indices using amap
    3. Updates atom names and residue names according to the mapping
    4. Assigns unique chain IDs (A-Z, then a-z) to each fragment
    5. Adds segment ID to all atoms
    6. Reorders atoms according to the force field specification
    7. Writes remapped PDB files
    
    PARAMETERS:
    -----------
    fragments : list of str
        List of fragment identifiers (e.g., ['tmp_103', 'tmp_118', ...])
        Expects fitted PDB files at: tmp_{pore}/{fragment}_sc.pdb
    
    pore : str
        Pore identifier (e.g., 'pore_1')
        Used to construct directory paths: tmp_{pore}/
    
    amap : pandas.DataFrame
        Atom mapping table with columns:
        - 'xindex': Original crystal atom index (1-based)
        - 'index': Target force field atom index
        - 'name': Force field atom name
        - 'resname': Force field residue name
        - 'resid': Force field residue ID
    
    seg : str, optional
        Segment ID for all atoms (default: 'POR1')
        Used in PDB SEGID field (positions 73-76)
    
    reorder_fn : callable, optional
        Function to reorder PDB content by atom number
        Default: reorder_pdb_by_atom_number (defined earlier in notebook)
        Should accept PDB content string and return sorted content
    
    RETURNS:
    --------
    dict : Statistics dictionary containing:
        - 'n_fragments': Number of fragments processed
        - 'n_chains': Number of unique chain IDs assigned
        - 'chain_ids': List of assigned chain IDs
        - 'mapped_files': List of output file paths
    
    OUTPUT FILES:
    --------------
    Creates mapped PDB files at:
        tmp_{pore}/{fragment}_mapped.pdb
    
    Each file contains:
    - Remapped atom indices (crystal → force field)
    - Updated atom names and residue names
    - Unique chain ID (A, B, C, ..., a, b, c, ...)
    - Segment ID (e.g., POR1)
    - Atoms reordered by force field specification
    
    WORKFLOW:
    ---------
    For each fragment:
    
    1. Read fitted PDB and extract CRYST1 line
    
    2. Process each PDB line:
       - HETATM/ANISOU lines:
         * Extract original atom index from columns 7-11
         * Look up mapping in amap dataframe
         * Reconstruct line with new index, atom name, residue name, chain ID, segment ID
       - CONECT lines:
         * Extract atom indices and map each one
         * Reconstruct with new indices
       - Other lines: Skip
    
    3. Reorder atoms according to force field (amap index field)
    
    4. Write remapped PDB file with CRYST1 + sorted atoms + END
    
    EXAMPLE:
    --------
    pore = 'pore_1'
    fragments = ['tmp_103', 'tmp_118', 'tmp_139', ...]
    
    stats = map_and_remap_pdb_fragments(
        fragments=fragments,
        pore=pore,
        amap=amap,
        seg='POR1',
        reorder_fn=reorder_pdb_by_atom_number
    )
    
    print(f"Processed {stats['n_fragments']} fragments")
    print(f"Assigned chain IDs: {stats['chain_ids']}")
    """
    
    # ========================================================================
    # SETUP
    # ========================================================================
    n_mol = len(fragments)
    
    # Generate chain IDs: A-Z (26 letters) then a-z (26 letters)
    chain_id_list = (string.ascii_uppercase + string.ascii_lowercase)[:n_mol]
    
    # Track output files and statistics
    mapped_files = []
    
    # ========================================================================
    # PROCESS EACH FRAGMENT
    # ========================================================================
    for f, chain in zip(fragments, chain_id_list):
        
        # Define input and output file paths
        file_sc = f'tmp_{pore}/{f}_sc.pdb'
        file_mapped = f'tmp_{pore}/{f}_mapped.pdb'
        
        # ====================================================================
        # READ FITTED PDB
        # ====================================================================
        with open(file_sc, 'r') as fin:
            pdb_lines = fin.readlines()
        
        # Extract CRYST1 line (first line with crystallography info)
        crystal_line = pdb_lines[0]
        
        new_atom_lines = []
        
        # ====================================================================
        # MAP ATOM INDICES AND UPDATE PDB RECORDS
        # ====================================================================
        for line in pdb_lines[1:]:
            header = line[:6]
            
            # Handle HETATM and ANISOU records
            if header in ['HETATM', 'ANISOU']:
                # Extract original atom index (columns 7-11, 1-based)
                xindex = int(line[6:11].strip())
                
                # Look up the mapping for this atom
                xatom_mapped = amap.loc[xindex]
                new_index = xatom_mapped['index']
                
                # Reconstruct PDB line with mapped information
                # PDB format (fixed-width columns):
                # - Cols 1-6: Record type (HETATM/ANISOU)
                # - Cols 7-11: Atom serial number (5 digits, right-aligned)
                # - Col 12: Space
                # - Cols 13-16: Atom name (4 chars, left-aligned)
                # - Cols 18-20: Residue name (3 chars, left-aligned)
                # - Col 22: Chain identifier (1 char)
                # - Cols 23-26: Residue seq number (4 digits)
                # - Cols 27-54: Coordinates (keep as-is)
                # - Cols 55-60: Occupancy/B-factor (keep as-is)
                # - Cols 77-78: Element symbol (keep as-is)
                # - Cols 73-76: Segment ID (4 chars)
                
                new_line = (
                    line[:6] +                                    # Record type
                    f"{new_index:>5} " +                          # New atom index
                    line[11:12] +                                 # Space
                    f"{xatom_mapped['name']:<4}" +                # Mapped atom name
                    f"{xatom_mapped['resname']:<3} " +            # Mapped residue name
                    chain +                                       # Chain ID (unique per fragment)
                    f"{xatom_mapped['resid']:>4}" +               # Residue ID
                    line[26:72] +                                 # Coordinates + occupancy (unchanged)
                    f"{seg}" +                                    # Segment ID
                    line[76:]                                     # Rest of line (unchanged)
                )
                
                new_atom_lines.append(new_line)
            
            # Handle CONECT records (connectivity information)
            elif header == 'CONECT':
                # CONECT records list bonded atoms: "CONECT atom_idx bonded_atom1 bonded_atom2 ..."
                new_conect_line = 'CONECT'
                
                # Extract all atom indices from the line
                atoms = [int(a) for a in line.split()[1:]]
                
                # Map each bonded atom to new index
                for xindex in atoms:
                    new_index = amap.loc[xindex]['index']
                    new_conect_line += f'{new_index:5d}'
                
                new_atom_lines.append(new_conect_line + '\n')
            
            # Skip other record types (TER, END, REMARK, etc.)
            else:
                continue
        
        # ====================================================================
        # REORDER ATOMS ACCORDING TO FORCE FIELD SPECIFICATION
        # ====================================================================
        # Build temporary PDB content to reorder
        temp_pdb_content = crystal_line + ''.join(new_atom_lines) + 'END\n'
        
        # Use reorder function to sort atoms by force field index
        pdb_sorted = reorder_fn(temp_pdb_content)
        
        # ====================================================================
        # WRITE REMAPPED PDB FILE
        # ====================================================================
        with open(file_mapped, 'w') as fout:
            fout.write(pdb_sorted)
        
        # Track output file and statistics
        mapped_files.append(file_mapped)
        print(f"Mapped {f:12} → Chain {chain} → {file_mapped}")
    
    # ========================================================================
    # RETURN STATISTICS
    # ========================================================================
    return {
        'n_fragments': n_mol,
        'n_chains': len(chain_id_list),
        'chain_ids': list(chain_id_list),
        'mapped_files': mapped_files,
        'segment_id': seg
    }


def generate_fitting_script(pore, reference_fragment, fragments_lst):
    """
    Generate PyMOL script for fitting all fragments to a reference structure.
    
    This function creates a PyMOL macro file (.pml) that automates the process
    of structural alignment. It uses a complete reference fragment as the
    template and fits all other fragments to it, capturing complete sidechain
    information in the process.
    
    CONTEXT:
    --------
    In crystallographic structures, some molecules may be incomplete or have
    missing atomic coordinates. By fitting incomplete fragments to a complete
    reference structure (aligned on backbone atoms), we can effectively
    "complete" the incomplete structures with proper sidechain coordinates
    from the reference.
    
    This is a critical step before atom mapping:
    1. All fragments must have the same number of atoms (129 for oligourea)
    2. Fitting adds missing atoms from reference to incomplete fragments
    3. Fitted structures are then remapped to force field indices
    
    PARAMETERS:
    -----------
    pore : str
        Pore identifier (e.g., 'pore_1')
        Used in:
        - Output filename: pymol_complete_{pore}.pml
        - Directory paths: tmp_{pore}/
        - PyMOL object selection
    
    reference_fragment : str
        Fragment identifier for complete reference structure
        Must be a PyMOL object loaded in the current session
        Characteristics:
        - Contains all 129 atoms (complete)
        - Used as template for fitting
        - Typically a fragment with no missing sidechains
        Example: 'tmp_331'
    
    fragments_lst : list of str
        List of all fragment identifiers to fit
        Format: ['tmp_103', 'tmp_118', 'tmp_139', ...]
        Processing:
        - Reference fragment is saved first
        - All other fragments are fitted to reference
        - Fitted fragments are saved with updated coordinates
    
    RETURNS:
    --------
    None (creates file: pymol_complete_{pore}.pml)
    
    OUTPUT FILE:
    -----------
    pymol_complete_{pore}.pml - PyMOL macro script containing:
    
    1. Backbone selection definition:
       select backbone_NO, (elem N and not name NAB) + (elem O and not (name OAE + name OAD))
       
       Selection logic:
       - Include all nitrogen atoms EXCEPT NAB (crystallographic artifact)
       - Include all oxygen atoms EXCEPT OAE and OAD (artifact atoms)
       - Excludes side-chain atoms for cleaner alignment
       - Focuses on backbone geometry
    
    2. Output formatting:
       set retain_order, 1
       
       Ensures PDB atom order is preserved (important for force field)
    
    3. Reference structure saving:
       save tmp_{pore}/{reference}_sc.pdb, {reference}
       
       Saves the complete reference structure unchanged
    
    4. Fragment fitting and saving (loop):
       For each fragment (except reference):
       fit {reference} and backbone_NO, {fragment} and backbone_NO
       save tmp_{pore}/{fragment}_sc.pdb, {reference}
       
       This fits backbone atoms of fragment to reference backbone,
       then saves the REFERENCE structure (with fragment's fitted atoms)
    
    PYMOL COMMANDS EXPLAINED:
    -------------------------
    select backbone_NO, ...
        Create named selection of backbone atoms (N and O only)
        Excludes artifacts (NAB, OAE, OAD)
        Used for clean structural superposition
    
    fit structure1, structure2
        Align structure1 to structure2 using atom pairs
        Performs rotation and translation
        Requires both structures to have same atom count in selection
    
    save output.pdb, object
        Export PyMOL object to PDB format
        Retains atomic coordinates from current viewport
        Preserves coordinate information from fitting
    
    WORKFLOW EXECUTION:
    -------------------
    1. Load all fragment PDB files into PyMOL before running
       (typically: load tmp_{pore}/tmp_*.pdb)
    
    2. Execute in PyMOL command line or batch mode:
       pymol -c -r pymol_complete_{pore}.pml
    
    3. Script outputs: fitted PDB files with complete sidechains
    
    4. Next step: Atom mapping (map_and_remap_pdb_fragments)
    
    EXAMPLE:
    --------
    # Setup data
    pore = 'pore_1'
    complete_f = 'tmp_331'  # Complete reference
    fragments = ['tmp_103', 'tmp_118', 'tmp_139', ...]
    
    # Generate script
    generate_fitting_script(pore, complete_f, fragments)
    
    # Output created: pymol_complete_pore_1.pml
    
    # Run in PyMOL
    pymol -c -r pymol_complete_pore_1.pml
    
    # Result: tmp_pore_1/ directory filled with:
    #   - tmp_331_sc.pdb (reference)
    #   - tmp_103_sc.pdb (fitted)
    #   - tmp_118_sc.pdb (fitted)
    #   - ... (all fragments fitted)
    
    TECHNICAL NOTES:
    ----------------
    - fit command modifies PyMOL viewport coordinates
    - save command exports modified coordinates
    - backbone_NO selection prevents fitting on side-chain differences
    - Artifacts (NAB, OAE, OAD) are excluded for geometric accuracy
    - Output PDB files retain PyMOL session coordinates
    - Fitted structures are ready for atom mapping in next step
    
    IMPORTANT:
    ----------
    - Ensure all fragments are loaded in PyMOL before execution
    - Reference fragment must have all atoms present
    - Output directory tmp_{pore}/ must exist or PyMOL will fail
    - PyMOL retain_order setting ensures reproducible atom ordering
    """
    
    # ========================================================================
    # CREATE PYMOL MACRO FILE
    # ========================================================================
    # Open file for writing PyMOL commands
    # Naming: pymol_complete_{pore}.pml
    # Example: pymol_complete_pore_1.pml
    with open(f'pymol_complete_{pore}.pml', 'w') as fout:
        
        # ====================================================================
        # BACKBONE ATOM SELECTION
        # ====================================================================
        # Define selection for backbone atoms (N and O only)
        # Excludes artifacts and side-chain atoms
        # NAB: Crystallographic artifact nitrogen (excluded)
        # OAE, OAD: Crystallographic artifact oxygens (excluded)
        # These atoms cause geometry mismatches and should be ignored in fitting
        fout.write('select backbone_NO, (elem N and not name NAB) + (elem O and not (name OAE + name OAD))\n')
        
        # ====================================================================
        # OUTPUT FORMAT SETTING
        # ====================================================================
        # Retain atom order from PyMOL (important for force field compatibility)
        # Ensures output PDB files preserve atomic sequence
        fout.write('set retain_order, 1\n')
        
        # ====================================================================
        # SAVE REFERENCE FRAGMENT (UNCHANGED)
        # ====================================================================
        # Write command to save the complete reference structure
        # Reference is saved first without any fitting
        # This serves as baseline for comparison
        # Output: tmp_{pore}/{reference_fragment}_sc.pdb
        fout.write(f'save tmp_{pore}/{reference_fragment}_sc.pdb, {reference_fragment}\n')
        
        # ====================================================================
        # FIT AND SAVE ALL FRAGMENTS
        # ====================================================================
        # Loop through all fragments and fit each to reference
        for f in fragments_lst:
            # Skip the reference fragment itself (already saved)
            if f != reference_fragment:
                # ============================================================
                # FIT FRAGMENT TO REFERENCE
                # ============================================================
                # PyMOL fit command:
                #   fit target_object and selection, reference_object and selection
                #
                # Action: Align target_object to reference_object using backbone atoms
                # - Calculates optimal rotation and translation
                # - Minimizes RMSD of backbone atoms
                # - Modifies viewport coordinates of target_object
                #
                # fit {reference} and backbone_NO, {fragment} and backbone_NO
                #   Fit reference object (using backbone) to fragment object (using backbone)
                #
                # After fit command:
                # - reference object position = fitted state
                # - Fragment coordinates are used as the fit target
                # - Reference object backbone is superimposed on fragment backbone
                # - Reference sidechains are now positioned relative to fragment
                #
                # This effectively "completes" the fragment with reference sidechains
                
                # ============================================================
                # SAVE FITTED REFERENCE AS FRAGMENT OUTPUT
                # ============================================================
                # Save the FITTED REFERENCE as the output for this fragment
                # This captures:
                # - Fragment's backbone geometry (from fitting)
                # - Reference's sidechain geometry (from reference structure)
                # Result: Complete structure with proper sidechains
                #
                # Output: tmp_{pore}/{fragment}_sc.pdb
                
                fout.write(f'fit {reference_fragment} and backbone_NO, {f} and backbone_NO; save tmp_{pore}/{f}_sc.pdb, {reference_fragment}\n')



def generate_merging_script(pore, n_atoms_by_fragment, output_file):
    """
    Generates a Python script to merge PDB fragment files into a single output file.

    This script uses PyMOL to load PDB files, merge them, and insert 'TER' records
    at the end of each fragment based on the specified number of atoms per fragment.

    Args:
        pore (str): Name or identifier of the pore, used for naming temporary files and directories.
        n_atoms_by_fragment (int): Number of atoms per fragment. Used to determine where to insert 'TER' records.
        output_file (str): Path to the output file where the merged structure will be saved.

    Returns:
        None: The script is written to a file named `merge_{pore}_fragments.py`.
    """

    with open(f'merge_{pore}_fragments.py', 'w') as fout:

        # Write the script header
        fout.write("""#! python3
from pymol import cmd
import os

# Set PyMOL to retain the order of atoms as in the input file
cmd.do('set retain_order, 1')

# Load all PDB files matching the pattern in the specified directory
""")

        # Generate commands to load each PDB file
        for f in glob.glob(f'tmp_{pore}/*_mapped.pdb'):
            fout.write(f"cmd.do('load {f}')\n")

        # Continue the script: save, process, and write the output file
        fout.write(f"""# Save the loaded structure to a temporary file
cmd.do('save tmp_mapped.pdb')

# Open the temporary file for reading
with open('tmp_mapped.pdb', 'r') as fin:
    lines = fin.readlines()

# Open the output file for writing
with open('{output_file}', 'w') as fout:
    for line in lines:
        # Write all lines except 'TER' records
        if not line.startswith('TER'):
            fout.write(line)
        # Check if the current line is the last atom of a fragment
        if line.startswith('ATOM  ') or line.startswith('HETATM'):
            # Extract the atom serial number from the line (columns 7-11)
            atom_serial = int(line[6:11].strip())
            if atom_serial % {n_atoms_by_fragment} == 0:
                fout.write('TER\\n')

# Remove the temporary file
os.remove('tmp_mapped.pdb')

""")
        

def merge_gro_files(files_lst, fout):
    """
    Merges multiple GROMACS .gro files into a single output file.

    This function reads atom coordinates and box information from each input .gro file,
    renumbers the atoms sequentially, and writes the merged data to an output file.
    The box information from the last file in the list is used for the output.

    Args:
        files_lst (list): List of paths to the input .gro files to be merged.
        fout (str): Path to the output .gro file where the merged data will be written.

    Returns:
        None: The merged data is written to the specified output file.
    """
    atoms = []  # List to store all atom lines from input files
    total_atoms = 0  # Counter for total atoms across all files
    box = None  # Variable to store the box information from the last file

    for fname in files_lst:
        with open(fname) as f:
            lines = f.readlines()

        # Extract the number of atoms from the second line of the .gro file
        n_atoms = int(lines[1].strip())
        # Extract the atom lines (lines 3 to 2 + n_atoms)
        atom_lines = lines[2:2+n_atoms]

        # Get the box information from the last line of the file
        box = lines[-1].strip()

        for line in atom_lines:
            # Parse each atom line into its components
            resid = line[0:5]      # Residue ID
            resname = line[5:10]   # Residue name
            atomname = line[10:15] # Atom name
            atomnum = int(line[15:20])  # Original atom number (not used in output)
            coords = line[20:].rstrip()  # Atom coordinates

            # Increment the total atom count and create a new line with the updated atom number
            total_atoms += 1
            new_line = f"{resid}{resname}{atomname}{total_atoms:5d}{coords}\n"
            atoms.append(new_line)

    # Write the merged data to the output file
    with open(fout, "w") as out:
        out.write("Merged system\n")
        out.write(f"{len(atoms)}\n")
        out.writelines(atoms)
        out.write(f"{box}\n")

    print(f"Written {fout} with {len(atoms)} atoms.")
