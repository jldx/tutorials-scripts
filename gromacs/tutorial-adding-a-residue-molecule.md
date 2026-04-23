# Tutorial: Adding a Residue/Molecule to GROMACS Files

---

## **1. Overview of Required Files**

To add a new residue or molecule, you may need to modify the following files:

- **`.rtp`**: Residue Topology Parameter file (defines atoms, bonds, angles, dihedrals, and impropers)
- **`.hdb`**: Hydrogen Database file (for automatic hydrogen addition)
- **`atomtypes.atp`**: Atom types file (if new atom types are introduced)
- **`ffnonbonded.itp`**: Non-bonded parameters file (if new atom types are introduced)
- **`ffbonded.itp`**: Bonded parameters file (if new bonded types are introduced)
- **`residuetypes.dat`**: Residue types file (classifies the residue, e.g., Protein, DNA, Ion)
- **`specbond.dat`**: Special bonds file (if the residue has special connectivity to other residues)
- **`forcefield.itp`**: Force field file (if new parameters are required)

---

## **2. Step-by-Step Process**

### **Step 1: Add the Residue to the `.rtp` File**

The `.rtp` file defines the topology of your residue. You can either copy an existing residue and modify it, or use an external tool to generate the topology and adapt it to the `.rtp` format.

**Example for Oligourea Alanine (UUA):**

```plaintext
[ UUA ]
; alanine
  [ atoms ]
        N      NH1 -0.4700   1
       HN        H  0.3100   1
       CA      CT1  0.0700   1
       HA      HB1  0.0900   1
       CB      CT3 -0.2700   2
      HB1      HA3  0.0900   2
      HB2      HA3  0.0900   2
      HB3      HA3  0.0900   2
       CK      CT2 -0.0200   3
      HK1      HA2  0.0900   3
      HK2      HA2  0.0900   3
       NU      NH1 -0.4700   3
       HU        H  0.3100   3
        C        C  0.5100   4
        O        O -0.5100   4
  [ bonds ]
        N    HN
        N    CA
       CA    HA
       CA    CK
       CK   HK1
       CK   HK2
       CK    NU
       NU    HU
       NU     C
        C     O
        C    +N
       CB    CA
       CB   HB1
       CB   HB2
       CB   HB3
  [ impropers ]
  ; ensure planarity of urea group
        C    NU    +N     O
        N    -C    CA    HN
```

- **`[atoms]`**: Lists atom names, types, charges, and groups.
- **`[bonds]`**: Defines bonds between atoms. `+N` refers to the N atom of the next residue.
- **`[impropers]`**: Defines improper dihedrals to maintain planarity or chirality.

---

### **Step 2: Add the Residue to the `.hdb` File**

The `.hdb` file allows GROMACS to automatically add hydrogens to your residue.

**Example for UUA:**

```plaintext
UUA        5
1       1       HN       N      CA      -C
1       1       HU      NU      C       CK
2       6       HK      CK      CA      NU
1       5       HA      CA      CK      CB     N
3       4       HB      CB      CA      C
```

- The first line (`UUA        5`) indicates the residue name and the number of hydrogen addition rules.
- Each subsequent line defines how a hydrogen is added: the number of hydrogens, the atom to which the hydrogen is attached, and the reference atoms for positioning.

---

### **Step 3: Add New Atom Types (if needed)**

If your residue introduces new atom types, add them to:

- **`atomtypes.atp`**: Defines atom types and their properties.
- **`ffnonbonded.itp`**: Defines non-bonded parameters (Lennard-Jones, etc.) for new atom types.

---

### **Step 4: Add New Bonded Parameters (if needed)**

If your residue requires new bonded parameters (bonds, angles, dihedrals), add them to:

- **`ffbonded.itp`**: Defines bonded parameters.

---

### **Step 5: Add the Residue to `residuetypes.dat`**

Classify your residue in `residuetypes.dat` (e.g., Protein, DNA, Ion).

**Example:**

```plaintext
UUA   Protein
```

---

### **Step 6: Update `specbond.dat` (if needed)**

If your residue has special connectivity to other residues (e.g., termini), update `specbond.dat`.

---

### **Step 7: Include Additional Parameters (if needed)**

If your residue requires additional parameters (e.g., for oligourea), include them in a separate `.itp` file and add an `#include` directive to `forcefield.itp`.

**Example:**

```plaintext
#include "oligourea_additional_parm.itp"
```

---

## **3. Summary Checklist**

1. Add the residue to the `.rtp` file.
2. Add hydrogen addition rules to the `.hdb` file.
3. Add new atom types to `atomtypes.atp` and `ffnonbonded.itp` (if needed).
4. Add new bonded parameters to `ffbonded.itp` (if needed).
5. Add the residue to `residuetypes.dat`.
6. Update `specbond.dat` for special connectivity (if needed).
7. Include additional parameters in `forcefield.itp` (if needed).

---

## **4. References**

- [GROMACS Manual: Topology Files](https://manual.gromacs.org/current/reference-manual/topologies/pdb2gmx-input-files.html)
- [GROMACS Manual: Force Fields](https://manual.gromacs.org/current/reference-manual/functions/force-fields.html)