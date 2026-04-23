# Tutorial: Solvating a Molecule in Methanol with GROMACS

This tutorial guides you through the process of solvating a molecule in methanol using GROMACS, starting from a PDB structure.

---

## **Prerequisites**
- GROMACS installed and available in your system path.
- A PDB file of your molecule.
- The CHARMM36m force field files (including `charmm36m.ff/forcefield.itp`).
- A GROMACS topology file for methanol (`MEOH.itp`).

---

## **Step 1: Generate Topology and Coordinates**
Convert your PDB file to a GROMACS structure file (`.gro`) and generate a topology:

```bash
gmx pdb2gmx -f <pdb structure> -o <gro structure> -ter
```
- Replace `<pdb structure>` with your input PDB file.
- Replace `<gro structure>` with your desired output GRO file.

---

## **Step 2: Create an Octahedral Box**
Define an octahedral box around your molecule:

```bash
gmx editconf -f <gro structure> -o <gro structure>_box.gro -bt octahedron -d 1.2
```
- This creates a box with a 1.2 nm distance between the solute and the box edge.

---

## **Step 3: Insert Methanol Molecules**
Fill the box with methanol molecules:

```bash
gmx insert-molecules -f <gro structure>_box.gro -ci methanol.gro -nmol 3000 -radius 2.0 -o <gro structure>_methanol.gro
```
- `-ci methanol.gro`: Specifies the methanol coordinate file.
- `-nmol 3000`: Inserts 3000 methanol molecules.
- `-radius 2.0`: Sets the minimum distance between inserted molecules.

---

## **Step 4: Generate an Index File**
Create an index file for the solvated system:

```bash
echo "q" | gmx make_ndx -f <gro structure>_methanol.gro -o index.ndx
```
- This creates `index.ndx`, which is useful for defining groups in simulations.

---

## **Step 5: Update the Topology File**
Edit your `topol.top` file:

1. **Include Methanol Topology:**
   Add the following line below `#include "./charmm36m.ff/forcefield.itp"`:
   ```plaintext
   #include "MEOH.itp"
   ```

2. **Update the `[ molecules ]` Section:**
   Add the number of methanol molecules (e.g., 3000):
   ```plaintext
   [ molecules ]
   ; Compound        #mols
   Protein_A         1
   MEOH             3000
   ```

---

## **Step 6: Adjust Position Restraints (Optional)**
If you are using position restraints, edit `posre.itp` and replace `1000` with a variable like `POSRES_FC_BB` for flexibility.

