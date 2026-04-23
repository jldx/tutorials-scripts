# Tutorial: Creating a Membrane System from Scratch in GROMACS

This tutorial guides you through the process of creating a membrane system from scratch in GROMACS when you cannot use CHARMM-GUI.

---

## **Prerequisites**
- GROMACS installed and available in your system path.
- A PDB file of your molecule (`mol.pdb`).
- CHARMM36 force field files.
- Python and PyMOL installed for structure manipulation.
- Script `reorder_atoms_POPC.py` for reordering atoms in the membrane.

---

## **Step 1: Prepare the Molecule**

### **1.1 Center the Molecule**
Center your molecule in a cubic box with a 1.2 nm margin:

```bash
gmx_mpi editconf -f mol.pdb -o mol_center.gro -bt cubic -d 1.2 -center 0 0 0
```
- Note the `new box vectors` from the output. This will be used to determine the membrane size.

### **1.2 Generate Topology**
Generate the topology for your molecule:

```bash
gmx_mpi pdb2gmx -f mol_center.gro -water tip3 -ter -o mol.gro
```

### **1.3 Adjust Position Restraints**
Edit `posre.itp` and replace `1000` with `POSRES_FC_BB`.

---

## **Step 2: Create the Membrane**

### **2.1 Build the Lipid Bilayer**
Use CHARMM-GUI to build a POPC bilayer of size 225 x 225 Å (based on the box vectors from Step 1.1). Download and extract `charmm-gui.tgz` to get `step5_input.pdb`.

### **2.2 Clean the Membrane PDB**
Remove water, ions, and unnecessary lines from the membrane PDB:

```bash
sed -i '/TIP3/d' step5_input.pdb
sed -i '/CLA/d' step5_input.pdb
sed -i '/SOD/d' step5_input.pdb
sed -i '/TER/d' step5_input.pdb
mv step5_input.pdb membrane.pdb
```

### **2.3 Center and Remove Hydrogens**
Center the membrane and remove hydrogen atoms:

```bash
gmx_mpi editconf -f membrane.pdb -o membrane_center.pdb -center 0 0 0
awk '$3 !~ /^H/' membrane_center.pdb > tmp.pdb && mv tmp.pdb membrane_center_noH.pdb
```

### **2.4 Reorder Atoms**
Reorder the atoms in the membrane PDB:

```bash
python3 reorder_atoms_POPC.py membrane_center_noH.pdb membrane_center_noH_reordered.pdb
```

### **2.5 Generate Membrane Topology**
Generate the topology for the membrane:

```bash
gmx_mpi pdb2gmx -f membrane_center_noH_reordered.pdb -o membrane_center.pdb -water tip3p
```

---

## **Step 3: Merge Molecule and Membrane**

### **3.1 Merge Structures**
Use PyMOL to merge the molecule and membrane:

```bash
pymol merge_mol_membrane.pml
```

### **3.2 Center the Merged System**
Center the merged system:

```bash
gmx_mpi editconf -f mol_membrane.pdb -o mol_membrane_center.gro -center 0 0 0
```

---

## **Step 4: Create the System Topology**

### **4.1 Edit Molecule and Membrane Topologies**
In `mol.itp` and `POPC.itp`, update the `[ moleculetype ]` name and include position restraints:

```plaintext
; Include Position restraint file
#ifdef POSRES
#include "posre_mol.itp"
#endif
```

### **4.2 Create `topol.top`**
Create the system topology file:

```plaintext
; Include forcefield parameters
#include "./charmm36m.ff/forcefield.itp"

; Include chain topologies
#include "toppar/mol.itp"
#include "toppar/POPC.itp"

; Include water topology
#include "./charmm36m_oligo.ff/tip3p.itp"

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
Mol               360
POPC              793
```

### **4.3 Create Index File**
Generate an index file for the system:

```bash
echo "q" | gmx_mpi make_ndx -f mol_membrane_center.gro -o index.ndx
```

---

## **Step 5: Pre-Minimization**
Run a pre-minimization of the system to relax any steric clashes.

---

## **Step 6: Solvate and Add Ions**

### **6.1 Solvate the System**
Add water to the system:

```bash
gmx_mpi solvate -cp npt.gro -p topol.top -o mol_membrane_solv.gro
```

### **6.2 Prepare for Ion Addition**
Generate a portable binary run file for ion addition:

```bash
gmx_mpi grompp -f ions.mdp -c mol_membrane_solv.gro -p topol.top -o ions.tpr
```

### **6.3 Add Ions**
Add ions to neutralize the system and set the salt concentration:

```bash
echo "20" | gmx_mpi genion -s ions.tpr -o mol_membrane_solv_ions.gro -p topol.top -pname NA -nname CL -conc 0.15
```
- `-pname NA -nname CL`: Specifies sodium (NA) and chloride (CL) ions.
- `-conc 0.15`: Sets the salt concentration to 0.15 M.
- `echo "20"`: Selects the solvent group for ion replacement.