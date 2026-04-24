# Creating Clean Pore Structures from Crystallographic Data

## Overview

This tutorial demonstrates how to extract and clean individual pore molecules from a dirty crystallographic structure for use in molecular dynamics simulations with the CHARMM36M oligourea force field.

**System Example:** H5 oligourea crystal structure (multiple pore molecules)

**Workflow:** Dirty PDB → Crystal expansion (VESTA) → Fragment splitting → Manual curation → Atom mapping → Clean individual pores

---

## Prerequisites


### Input Files
- `H5.pdb`: Original crystallographic structure
- `H5_atoms.map`: Atom index mapping file (relates original to mapped indices)
- Knowledge of fragment properties (e.g., 129 atoms per oligourea molecule)

---

## Step-by-Step Workflow

Open a `python` session

```python
# In Python:
import utils
import os
```

### Step 1: Create Crystal Sample (VESTA)

**Purpose:** Expand the crystallographic unit cell to create a larger sample containing multiple pore molecules.

**Software:** VESTA (Visualization, Electronic and Structure Analysis)

**Process:**
1. Open the original PDB file in VESTA: `H5.pdb`
2. Use the **Boundary** tool to define the supercell:
   - X-axis: 0 to 4 (expand 4 times in X)
   - Y-axis: 0 to 4 (expand 4 times in Y)
   - Z-axis: 0 to 1 (keep 1 unit in Z)
3. Save the expanded structure as: `H5_crystal_plus.pdb`

**Output:** `H5_crystal_plus.pdb` (larger sample with multiple pore molecules)

**Why this step?**
- Creates a system containing all distinct pores in the crystal
- Facilitates identification and isolation of individual pores
- Provides context for crystallographic symmetry

---

### Step 2: Split Structure into Molecular Fragments

**Purpose:** Automatically decompose the large structure into individual connected molecular fragments.

**Tool:** `split_fragments.py` (MDAnalysis-based)

**Command:**
```python
# In Python:
os.mkdir('tmp')
n_molecules = split_molecular_fragments('H5_crystal_plus.pdb' 'tmp/')
```

**Output:** 
- Directory `tmp/` containing ~500+ individual molecule PDB files
- Naming: `tmp_0.pdb`, `tmp_1.pdb`, ..., `tmp_N.pdb`

**Key Insight:**
- Each fragment represents a distinct molecule/chain
- Some fragment numbers may represent artifacts (hydrogens, duplicate objects)
- These will be manually curated in the next steps

---

### Step 3: Manual Curation in PyMOL (Interactive)

**Purpose:** Manually identify which fragments form complete pores and clean up crystallographic artifacts.

#### 3a. Identify Pore Fragments

**Process:**
1. Open PyMOL and load all fragment files:
   ```pymol
   load tmp/tmp_*.pdb
   ```
   
2. Visually inspect fragments to identify which ones form pores:
   - Look for 3D ring structures (characteristic of oligourea pores)
   - Pores are typically composed of 30-40 individual fragment objects
   - Non-pore fragments: isolated side chains, crystallographic artifacts

3. Select all fragments forming a pore in PyMOL by selecting chains

4. Get the object list:
   ```pymol
   # In PyMOL:
   print(cmd.get_object_list('pore_fragments'))
   ```

#### 3b. Create Pore Group

**Process:**
1. Convert the PyMOL selection list to a group command:
   ```python
   # In Python:
   pymol_list = ['tmp_103', 'tmp_118', 'tmp_139', ...]  # from print output
   group_command = 'group pore_1, ' + ' '.join(pymol_list)
   print(group_command)
   ```

2. Execute in PyMOL:
   ```pymol
   # In PyMOL:
   group pore_1, tmp_103 tmp_118 tmp_139 tmp_154 tmp_171 tmp_182 ...
   ```

#### 3c. Fix Crystallographic Artifacts

**Issue:** VESTA sometimes splits a single molecule into multiple objects (bonded fragments considered separate).

**Affected Fragments:** `tmp_7`, `tmp_22`, `tmp_43`, `tmp_58`, `tmp_75`, `tmp_86`

**Process:**
1. Identify the fragment pairs that should be merged:
   ```pymol
   # In PyMOL:
   select tmp_7
   # Visual inspection shows this contains 2 separate molecules bonded together
   ```
2. Delete bond with PyMOL pluging Builder

3. Remove unwanted hydrogens:
   ```pymol
   # In PyMOL:
   remove (name H*)  # Remove all hydrogens
   ```

4. Separate the two fragments into distinct objects.

5. Update the pore group:
   ```pymol
   # In PyMOL:
   group pore_1, pore_1 tmp_7a tmp_7b
   delete tmp_7  # Remove the original artifact
   ```

6. Repeat for all affected fragments (tmp_22, tmp_43, tmp_58, tmp_75, tmp_86)

**Repeat for other pores** if your structure contains multiple pores.

⚠️ **Don't close PyMOL**

---

### Step 4: Structure Alignment and Sidechain Completion

**Purpose:** Complete incomplete fragments by aligning them to a reference fragment with all atoms present.

**Why alignment is needed:**
- Some crystal fragments are missing sidechains (incomplete)
- Alignment to a complete reference fills in missing atoms
- Ensures uniform atom count (129 atoms/fragment)

```python
# In Python:
pore = 'pore_1'
complete_f = 'tmp_331' # a complete fragment with 129 atoms
# Get the list of objects in the group
# Execute in PyMOL commands panel: print(cmd.get_object_list('pore_1')
fragments = ['tmp_103', 'tmp_118', 'tmp_139', 'tmp_154', 'tmp_171', 'tmp_182', 
             'tmp_202', 'tmp_214', 'tmp_235', 'tmp_247', 'tmp_267', 'tmp_278', 
             'tmp_295', 'tmp_310', 'tmp_331', 'tmp_346', 'tmp_363', 'tmp_374', 
             'tmp_394', 'tmp_406', 'tmp_427', 'tmp_439', 'tmp_459', 'tmp_470', 
             'tmp_22a', 'tmp_22b', 'tmp_43a', 'tmp_43b', 'tmp_7a', 'tmp_7b', 
             'tmp_58a', 'tmp_58b', 'tmp_75a', 'tmp_75b', 'tmp_86a', 'tmp_86b']

utils.generate_fitting_script(pore, complete_f, fragments)

os.mkdir(f'tmp_{pore}')
```

**How it works:**
1. Select backbone atoms (N, O) excluding specific artifact atoms (NAB, OAE, OAD)
2. For each fragment, fit it to the reference using backbone alignment
3. Save the fitted fragment (which now has complete sidechain from reference)

Execute in PyMOL: `@pymol_complete_pore_1.pml`

**Output:** Directory `tmp_pore_1/` containing:
- `tmp_331_sc.pdb`: Complete reference
- `tmp_103_sc.pdb` through `tmp_86b_sc.pdb`: Aligned fragments with complete sidechains

---

### Step 5: Atom Index Mapping and PDB Remapping

**Purpose:** Map atom indices from crystal structure to oligourea force field numbering system.

**Input File:** `H5_atoms.map` (tab-separated)
```
xindex    index    name    resname    resid
1         1        C1      OUA        1
2         2        C2      OUA        1
...       ...      ...     ...        ...
```

```python
# In Python:
seg = 'POR1'  # Segment identifier for the pore
amap = utils.read_atom_map('H5_atoms.map')

# Call the function to map and remap all fragments
mapping_stats = utils.map_and_remap_pdb_fragments(
    fragments=fragments,
    pore=pore,
    amap=amap,
    seg=seg,
    reorder_fn=utils.reorder_pdb_by_atom_number
```

**Key operations:**
1. **Map atom indices**: Original crystal index → force field index
2. **Rename atoms**: Update atom names according to oligourea naming scheme
3. **Update chain IDs**: Assign unique chain letter (A-Z, then a-z) to each fragment
4. **Update segment ID**: All fragments get segment ID "POR1" (or appropriate pore ID)
5. **Reorder atoms**: Ensure atoms appear in the order specified by H5_atoms.map

**Output:** `tmp_pore_1/tmp_*_mapped.pdb` files with corrected indices and force field names

---

### Step 6: Merge Mapped Fragments into Single Pore File

**Purpose:** Combine all individual mapped fragments into one complete pore PDB file.

This will generate `merge_pore_1.py`:

```python
# In Python:
n_atoms_by_fragment = 129
utils.generate_merging_script(pore, n_atoms_by_fragment, output_file=f'{pore}_mapped.pdb')
```
Execute in a new terminal: `pymol -c -r merge_pore_1.py`

**Output:** `pore_1_mapped.pdb` - Complete pore structure with:
- All 30-40 oligourea molecules properly labeled
- Unique chain IDs (A-Z, a-z)
- Segment ID (POR1)
- Correct atom indices and names
- TER records between molecules

**Pore PDB properties:**
- ✓ No crystallographic artifacts
- ✓ Complete sidechains on all molecules
- ✓ Correct force field atom names/indices
- ✓ Proper chain and segment labeling
- ✓ Ready for GROMACS force field assignment (CHARMM36M_oligourea)

---

### Step 7: Generate the Pore Individual Gromacs Files

**Purpose:** Convert the pdb to gromacs format and generate the gromacs input files

Create temporary folders for each pore convertion where the folder contains the `charmm36m_oligourea.ff`

```bash
gmx_mpi pdb2gmx -f pore_1_mapped.pdb -o pore_1_mapped.gro
```

**Output:** 
- `pore_1_mapped.gro`: Gro structure
- `posres_Protein_chain_[].itp`: Positions restraints for all chains
- `topol_Protein_chain_[].itp`: Topology file for all chains
- `topol.top`: Global topology file

### Step 8: Create the Multipore Gromacs Files

**Purpose:** Create the multipore Gromacs input files 

#### 8a. Generate the multipore gro structure

```python
# In Python:
gro_lst = ['pore_1_mapped.gro', 'pore_2_mapped.gro', 'pore_3_mapped.gro',
           'pore_4_mapped.gro', 'pore_5_mapped.gro', 'pore_6_mapped.gro',
           'pore_7_mapped.gro']

utils.merge_gro_files(gro_lst, 'H5_crystal_plus.gro')
```
**Output:** `H5_crystal_plus.gro` - The multipore complete structure

#### 8b. Create the multipore `topol.top`

**Previous knowledge:** The number of fragments in the multipore

```
; Include forcefield parameters
#include "./charmm36m_oligourea.ff/forcefield.itp"

; Include chain topologies
#include "topol_chains.itp"

; Include water topology
#include "./charmm36m_oligourea.ff/tip3p.itp"

#ifdef POSRES_WATER
; Position restraint for each water oxygen
[ position_restraints ]
;  i funct       fcx        fcy        fcz
   1    1       1000       1000       1000
#endif

; Include topology for ions
#include "./charmm36m.ff/ions.itp"

[ system ]
; Name
Protein in water

[ molecules ]
; Compound        #mols
Protein             252
```

#### 8c. Create the multipore `posre_chain.itp` and `topol_chain.itp`

The multipore is an homo-oligomer. Each chain as the same topology. 

Copy/Paste a `posre_Protein_chain_[].itp` and a `topol_Protein_chain_[].itp` from one of the temporary folder of a pore and rename them to `posre_chain.itp` and `topol_chain.itp`

## Final Output

After completing all steps:

```
working_directory/
├── pore_1.pdb          # ✓ Clean, ready for simulation
├── pore_2.pdb          # ✓ (if structure has multiple pores)
└── ...
└── pore_1.gro
└── ...
└── H5_multipore.gro
└── topol_chains.itp
└── posre_chain.itp
└── topol.top
└── ...
```
## Next steps

What you have to do: solvation, inclusion in membrane...
