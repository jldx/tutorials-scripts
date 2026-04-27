# Coarse-Graining Proteins with Martini: From All-Atom to Coarse-Grained Models

## Objective

This tutorial describes the process of converting high-resolution all-atom protein structures into coarse-grained Martini force field representations. The workflow enables efficient molecular dynamics simulations of larger systems by reducing computational complexity. Optionally, the coarse-grained protein can be embedded within a coarse-grained membrane model for lipid-protein interaction studies.

## Requirements

The following tools must be installed and available:

- **martinize2**: The main tool for converting all-atom protein structures to Martini coarse-grained representations ([Download](https://vermouth-martinize.readthedocs.io/en/latest/index.html))
- **insane**: Utility package for building membrane systems and embedding proteins (optional, only if membrane insertion is desired) ([Download](https://github.com/Tsjerk/Insane))


## Step 1: Create Martini Model of Protein Structure

Convert the all-atom protein structure to a coarse-grained Martini representation and generate GROMACS topology files using the following command:

```bash
martinize2 -f <topology> -o topol.top -x <martini_model_output> -dssp -p backbone -ff martini3001
```

**Parameter explanations:**
- `-f <topology>`: Input all-atom protein structure (PDB or GRO format)
- `-o topol.top`: Output GROMACS topology file
- `-x <martini_model_output>`: Output coarse-grained structure file
- `-dssp`: Include secondary structure prediction using DSSP (improves coil vs. secondary structure mapping)
- `-p backbone`: Preserve backbone connectivity information
- `-ff martini3001`: Specify Martini force field version (3.0.01)


## Step 2: Embed Coarse-Grained Protein into a Membrane System

To insert the coarse-grained protein into a pre-equilibrated coarse-grained lipid bilayer, use the insane package with the following command:

```bash
insane -f <martini_model_output> -o <martini_model_output_in_membrane> -box 30,30,60 -l POPE:75 -l POPG:25 -center -dm 4.5 -sol W -p topol.top
```

**Parameter explanations:**
- `-f <martini_model_output>`: Input coarse-grained protein structure (from Step 1)
- `-o <martini_model_output_in_membrane>`: Output system with embedded protein
- `-box 30,30,60`: Simulation box dimensions in nanometers (X, Y, Z)
- `-l POPE:75 -l POPG:25`: Lipid composition specification (75% POPE, 25% POPG in this example)
- `-center`: Center the protein in the simulation box
- `-dm 4.5`: Distance threshold for membrane inclusion (prevents lipids within 4.5 nm of protein from being placed)
- `-sol W`: Solvent type (W = water)
- `-p topol.top`: Path to the GROMACS topology file

**Example system:** This command creates a 75:25 POPE:POPG lipid bilayer with the coarse-grained protein embedded centrally.