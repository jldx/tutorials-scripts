# Homology Modeling with AlphaFold: Completing Crystallographic Structures

## Objective

This tutorial describes how to leverage AlphaFold predicted models to complete experimental crystallographic structures, particularly when regions are missing or poorly resolved due to conformational disorder or experimental limitations.

## Important Considerations and Limitations

1. **Disordered Region Modeling**: AlphaFold has inherent limitations when predicting disordered regions. These flexible segments may be better modeled separately or treated with additional caution. Consider whether independent structural modeling of these regions is more appropriate.

2. **Conformational Bias**: AlphaFold predictions reflect the statistical distribution of conformational states in the PDB training data. Be aware that the model may preferentially predict the most commonly observed conformations, potentially overlooking rare or functionally important states.

3. **Model Quality Assessment**: Always validate the quality and confidence of AlphaFold predictions before using them to complete experimental structures. Check predicted local distance difference test (pLDDT) scores and other confidence metrics.

4. **Experimental Structure Context**: The nature of your experimental data significantly affects the reliability of AF-based completion. For cryo-EM structures, carefully review the original publication to assess map quality and evaluate whether the experimental structure used for fitting was appropriate. Poor fitting can propagate errors into the completed model.

## Workflow Overview

The workflow follows a straightforward three-step process to integrate AlphaFold predictions with experimental structures:

### Step 1: Structural Alignment and Chain Nomenclature

Using PyMOL or equivalent structural alignment software, align the AlphaFold predicted models with the original crystallographic structures. Ensure consistent chain naming to facilitate downstream analysis and interpretation.

### Step 2: Vacuum Energy Minimization

Perform an energy minimization of the completed structure in vacuum using GROMACS, AMBER, or another molecular dynamics/molecular mechanics package. This step allows the integrated structure to relax into a more stable conformation following the superposition.

### Step 3: Quality Control and Validation

Assess the quality of the completed structure using Procheck, Molprobity, or similar validation tools to identify potential geometry issues, clashes, or other structural anomalies introduced during the completion process.

## Materials and Methods Documentation

When reporting homology modeling results, document the following information:

- **Sequence information**: UniProt accession code of the modeled sequence
- **Experimental structures**: PDB codes and experimental methods (X-ray crystallography, cryo-EM, NMR, etc.)
- **AlphaFold model identifier**: Database model code or version used for predictions
- **Alignment protocol**: Complete PyMOL script or equivalent alignment method used to superpose structures