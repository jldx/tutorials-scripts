# Tutorial: Generating Universal Forcefields, Topologies, and Assembly Workflow for Oligourea Foldamers

**Goal:**
This tutorial guides you through the process of generating forcefields, topologies, and assembling workflows for simulating oligourea foldamers. The workflow is designed for oligoureas capped with Boc or NHMe (or uncapped), focusing on both terminal and internal residues.

## Step 1: Generate Initial Forcefields and Topologies

### Tools Required:
- [Charmm-Gui](https://www.charmm-gui.org/)

### Workflow:

1. **Ligand Reader & Modeler Module**
   - Draw your oligourea molecule and save the `smiles` file.
   - Generate temporary Cgenff files and download the `charmm-gui.tgz` archive.
   - Extract `ligandrm.pdb` and `drawing_3D.mol2` from the archive.

2. **Membrane Builder Module** (for lipid interactions)
   - Load `ligandrm.pdb` into the Membrane Builder.
   - Use `drawing_3D.mol2` to generate CHARMM top & par files:
     - **Option 1**: Use CHARMM General Force Field (Cgenff).
     - **Option 2**: Use Antechamber (gaff2).
   - Construct a bilayer membrane (e.g., 200 POPC / 200 POPS).
   - Add basic ion types (e.g., NaCl).
   - Save the Gromacs files and rename the output to `cgenff.tgz` or `gaff2.tgz` based on the forcefield used.

3. **Extract Relevant Files**
   - From the `gromacs/toppar` folder:
     - Rename `forcefield.itp` to `[GxG]_forcefield.itp` and move it to `GxG/ff_[]`.
     - Move `GxG.itp` to `GxG/mol_[]`.

---

## Step 2: Generate Overall Forcefields and Topologies for Residues

### Prerequisites:
- Manually construct the residue mapping (user-defined ID and atom naming). Refer to `GAG.map` for examples.

### Steps:

1. **Generate Residue Topologies**
   - Use the script `extract_residues_itp.py` to generate residue topologies from GxG.itp and GxG.map
   - Output: all `UUx.itp`.

2. **Generate Global Forcefields**
   - Use the script `merge_forcefields_itp.py` to generate the global forcefield. This will merge all individual `forcefield.itp` from all GxG to one.
   - Output: `forcefield.itp`.

---

## Step 3: Compare Generated Parameters with Existing Protein Parameters

### Objective:
Compare the parameters of residues generated as `GxG` with their equivalents in standard protein forcefields (Charmm36m, ff19SB).

### Data Sources:
- **Charmm36m**: Retrieve parameters from `src/charmm36_ljpme-jul2022.ff.tgz`:
  - `aminoacids.rtp`
  - `charm36m_ffnonbonded.itp`
  - `charm36m_ffbonded.itp`
- **Amber ff19SB**: Use Gromacs-formatted files from `src/1T45_ff19sb` (generated via Charmm-Gui).

### Compiled Data:
- Atom parameters (name, type, charge, sigma, epsilon) are compiled in `src/ff_atomtypes.csv`.

---

## Summary of Key Files and Scripts

| File/Script                | Description                                                                 |
|-----------------------------|-----------------------------------------------------------------------------|
| `residues_itp.py`          | Generates residue topologies.                                               |
| `forcefield_itp.py`        | Generates global forcefields.                                               |
| `ff_atomtypes.csv`         | Compiled atom parameters for comparison.                                    |
| `charmm36_ljpme-jul2022.ff.tgz` | Charmm36m parameters for amino acids.                                |
| `1T45_ff19sb`              | Amber ff19SB parameters in Gromacs format.                                 |

---

## Next Steps
- Use the generated forcefields and topologies in your molecular dynamics simulations.
- Validate the parameters by comparing simulation results with experimental data.