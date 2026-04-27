#!/usr/bin/env python3

"""
================================================================================
COMPARATIVE PROTEIN STRUCTURE MODELING USING MODELLER
================================================================================

Purpose:
    This script performs comparative (homology) protein structure modeling
    using the MODELLER software package. It generates multiple protein models
    based on template structures and evaluates them using the DOPE (Discrete
    Optimized Protein Energy) scoring function.

Workflow:
    1. Define a custom modeling class with symmetry restraints
    2. Generate multiple comparative models from sequence-template alignment
    3. Evaluate models using DOPE scoring
    4. Rank and report the best models
    5. Organize output files into a models directory

Prerequisites:
    - MODELLER software installed and configured
    - Sequence-template alignment file (aln.ali)
    - Template structure (PDB format)
    - Target sequence identifier (UniProt code or similar)

Key Concepts:
    - Symmetry Restraints: Enforces identical backbone geometry between chains
    - DOPE Score: Empirical scoring function for model quality assessment
    - Model Selection: Automated ranking enables identification of best models

================================================================================
"""

# ===== IMPORTS =====
from modeller import *              # Core MODELLER classes and functions
from modeller.automodel import *    # Automated comparative modeling framework
import subprocess                   # For executing shell commands
import sys                          # For Python version compatibility

# ===== CONFIGURATION =====
log.verbose()    # Enable verbose logging for detailed modeling output


# ===== CUSTOM MODELING CLASS =====
# Override default MODELLER behavior to add domain-specific restraints
class MyModel(AutoModel):
    """
    Custom AutoModel subclass implementing domain-specific structural restraints.
    
    This class enforces symmetry constraints between chains and reports
    on restraint violations after model generation.
    """
    
    def special_restraints(self, aln):
        """
        Apply custom symmetry restraints between chain pairs.
        
        Constrains chains A and B to have identical backbone geometry
        (C-alpha atoms only) to reduce computational complexity and enforce
        biological symmetry. This is useful for proteins with homologous
        or symmetric subunits.
        
        Parameters:
            aln: MODELLER alignment object containing chain information
        """
        # Select C-alpha atoms from chain A
        s1 = Selection(self.chains['A']).only_atom_types('CA')
        
        # Select C-alpha atoms from chain B
        s2 = Selection(self.chains['B']).only_atom_types('CA')
        
        # Apply symmetry restraint with weight 1.0 (perfect symmetry)
        # Restricts distance calculations to CA atoms only, significantly
        # reducing the computational burden compared to all-atom restraints
        self.restraints.symmetry.append(Symmetry(s1, s2, 1.0))
    
    def user_after_single_model(self):
        """
        Report symmetry violations after model building.
        
        Called after each individual model is generated. Reports any
        symmetry restraint violations exceeding 1.0 Angstrom threshold,
        allowing identification of modeling issues.
        """
        self.restraints.symmetry.report(1.0)

# ===== ENVIRONMENT INITIALIZATION =====
# Create a new MODELLER environment for this modeling session
# The environment manages all configuration, file I/O, and computational settings
env = Environ()

# ===== ENVIRONMENT CONFIGURATION =====
# Configure how input PDB files are parsed and processed
env.io.hetatm = False       # Ignore HETATMs (non-standard residues)
env.io.hydrogen = True      # Include hydrogen atoms in models (for detailed energetics)


# ===== MODEL GENERATION SETUP =====
# Create AutoModel instance with custom MyModel class to perform
# comparative modeling from sequence-template alignment
a = MyModel(env,
            alnfile  = 'aln.ali',                      # Input: Sequence-template alignment
            knowns   = 'template',                     # Input: Template structure identifier
            sequence = 'P75831',                       # Input: Target sequence identifier (UniProt)
            assess_methods=(assess.DOPE))              # Use DOPE for model quality assessment

# ===== MODEL GENERATION PARAMETERS =====
a.starting_model = 1                # Generate models starting from index 1
a.ending_model   = 100              # Generate up to 100 models
                                    # (higher number = more sampling, longer computation)

# ===== PERFORM COMPARATIVE MODELING =====
# Execute the actual modeling pipeline:
# 1. Read template structures
# 2. Align target sequence to templates
# 3. Build 100 models with varying optimization starting conditions
# 4. Evaluate each model using DOPE scoring function
a.make()

# ===== OUTPUT CLEANUP =====
# MODELLER generates multiple intermediate and diagnostic files during modeling.
# These are removed to keep the working directory clean, retaining only final models.

# Remove MODELLER diagnostic files (.D* = distance constraints)
result = subprocess.run('rm P75831.D*', shell=True, check=True, text=True, capture_output=True)

# Remove MODELLER assessment files (.V* = violation reports)
result = subprocess.run('rm P75831.V*', shell=True, check=True, text=True, capture_output=True)

# Move all final model PDB files (.pdb) to a dedicated models directory
# for organized storage and easier downstream analysis
result = subprocess.run('mv P75831.*.pdb models', shell=True, check=True, text=True, capture_output=True)


# ===== MODEL EVALUATION AND RANKING =====
# Extract successfully generated models from the modeling run output.
# The output includes both successful models and failed attempts,
# so we filter to keep only complete, non-failed models.
ok_models = [x for x in a.outputs if x['failure'] is None]

# ===== DOPE SCORING RANKING =====
# DOPE (Discrete Optimized Protein Energy) score provides an empirical
# quality assessment for comparative models. Lower scores indicate
# better model quality and structural reliability.
key = 'DOPE score'

# Sort models by DOPE score (lowest/best scores first)
# Python 3 compatibility: modern sort with key function
if sys.version_info[:2] == (2,3):
    # Python 2.3 fallback (deprecated but maintained for compatibility)
    ok_models.sort(lambda a,b: cmp(a[key], b[key]))
else:
    # Python 3+ standard approach
    ok_models.sort(key=lambda a: a[key])


# ===== RESULTS REPORTING =====
# Write ranking results to file for downstream analysis and archival
with open('models_DOPE.txt', 'w+') as f:
    for m in ok_models:
        # Format: model_name [tab] DOPE_score (3 decimal places)
        # This enables easy parsing and sorting by analysis scripts
        f.write(f"{m['name']}\t{m[key]:.3f}\n")


# ===== BEST MODEL IDENTIFICATION =====
# Extract and report the top-ranked model
# This is typically used as the starting point for further refinement
# (e.g., MD simulation, structure validation, functional studies)
m = ok_models[0]
print("Top model: %s (DOPE score %.3f)" % (m['name'], m[key]))

