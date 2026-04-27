# Homology Modeling with MODELLER: Building Complete Protein Structures

## Objective

This tutorial describes how to use MODELLER, a leading computational tool for comparative protein structure modeling, to build complete three-dimensional structures based on sequence homology to experimentally determined template structures. This approach is particularly valuable when experimental structures are incomplete, contain missing loops, or lack certain structural regions.

## What is MODELLER?

MODELLER is a bioinformatics software package that automatically generates three-dimensional protein models based on experimental structures (templates) and amino acid sequence alignment. The tool combines template structure information with spatial restraints derived from homologous sequences to construct energetically favorable protein models.

## Installation

### Option 1: Conda (Recommended)

Install MODELLER from the Salilab conda channel:

```bash
conda install -c salilab modeller
```

Access the official package: https://anaconda.org/channels/salilab/packages/modeller/overview

### Option 2: Direct Installation

MODELLER can also be installed from the official Salilab website (www.salilab.org/modeller)

### License Key

MODELLER requires a valid license key for full functionality. A free academic license key is available:

```
License Key: MODELIRANJ
```

Configure this key in your MODELLER environment following the installation instructions.

## Modeling steps

### Extract the sequence of the template

Extract the amino acid sequence from your template PDB structure using the `sequence_from_pdb.py` script:

```bash
python3 sequence_from_pdb.py template.pdb template.fasta
```

**Parameters:**
- `template.pdb`: Input PDB structure file containing the template protein
- `template.fasta`: Output FASTA sequence file (will be created)

**Important Limitations:**
- This script extracts only standard 20 amino acids
- Post-translational modifications (PTMs) are NOT included in the output
- Non-standard residues and heteroatoms are ignored
- For information on handling modified residues, refer to:
  - `restyp.lib` (standard residue definitions)
  - `restyp_accelrys.lib` (extended residue library)

### Run multiple sequence alignment

Perform a multiple sequence alignment between your target sequence and the template sequence using a dedicated alignment tool:

**Recommended tools:**
- **Clustal Omega**: Web server at https://www.ebi.ac.uk/jdispatcher/msa/clustalo (recommended for most users)
- **BLAST**: For identifying related sequences in databases
- **MUSCLE**: For local alignment with custom parameters
- **T-Coffee**: For high-quality alignments with conservation information

**Alignment workflow:**
1. Extract target sequence (the protein you want to model)
2. Extract template sequence from PDB using `sequence_from_pdb.py`
3. Align both sequences using your chosen tool
4. Save alignment in PIR format

**Create the alignment input file (`aln.ali`):**

Convert your sequence alignment into PIR (Protein Information Resource) format as required by MODELLER. Detailed formatting instructions are available in the [official MODELLER PIR documentation](https://salilab.org/modeller/9v8/manual/node454.html).

Example PIR format:
```
>P1;template
structure:template.pdb:::::::0.0:0.0
MKTIIALSYIFCLVNADYKDDDDDKQN
*

>P1;target
sequence::::::::
MKTIIALSYIFCLVNADYKDDDDDKQN
*
```

### Run modeller

Configure the main modeling script with your specific parameters, then execute the comparative modeling:

**Configuration (in `modelling.py`):**

Edit the following parameters in `modelling.py` to match your project:

```python
alnfile  = 'aln.ali'          # Path to your alignment file (PIR format)
knowns   = 'template'         # Template structure identifier (as in aln.ali)
sequence = 'target'           # Target sequence identifier (as in aln.ali)
```

**Execute the modeling:**

```bash
python3 modelling.py
```

**What the script performs:**

1. **Alignment Reading**: Loads the sequence-template alignment from `aln.ali`
2. **Template Loading**: Retrieves template structures and identifies homologous regions
3. **Model Generation**: Creates 100 comparative models with varying optimization starting conditions
4. **Energy Scoring**: Evaluates each model using DOPE (Discrete Optimized Protein Energy) function
5. **Quality Ranking**: Sorts models by DOPE score (lower = better quality)
6. **Output Organization**: Writes results to `models_DOPE.txt` and organizes PDB files in `models/` directory

### Expected Output

The MODELLER script generates the following files and directories:

**Main outputs:**
- `models_DOPE.txt`: Tab-separated file with all models ranked by DOPE score (best to worst)
  - Format: `model_name [tab] dope_score`
  - Example: `target.B99990001.pdb  -0.56`
- `models/`: Directory containing all generated PDB structure files
  - Files named `target.B9999000X.pdb` (where X is the model number)

**Intermediate files (automatically cleaned up):**
- `.D*` files: Distance constraint files (removed after modeling)
- `.V*` files: Violation report files (removed after modeling)

## Quality Assessment

After model generation, validate the quality of predicted structures using independent assessment tools. Quality assessment is critical for identifying reliable models suitable for downstream applications.

### Understanding DOPE Scores

The DOPE (Discrete Optimized Protein Energy) score is built into MODELLER and provides a rapid quality estimate:

- **Lower scores indicate better quality**
- Scores are typically negative (range: -0.5 to -1.5 for proteins)
- Use the top 5-10% of models (lowest DOPE scores) for detailed analysis
- DOPE scores from `models_DOPE.txt` are already sorted for easy model selection

### Recommended Validation Tools

**Procheck**
- Analyzes backbone and side chain geometry
- Generates Ramachandran plots for phi-psi angle distribution
- Identifies steric clashes and unusual conformations
- Best for: Quick overall quality assessment

**Molprobity**
- Comprehensive validation with multiple criteria
- Detects clashes, geometry errors, and rotamer issues
- Provides detailed residue-by-residue reports
- Best for: Detailed validation and error identification

**QMEAN**
- Machine learning-based quality scoring
- Predicts local and global model quality
- Provides confidence estimates
- Best for: Independent quality confirmation

**DOPE Score (Built-in)**
- Already calculated during MODELLER run
- Quick initial model ranking
- Should be complemented with other validation tools

### Quality Criteria for Acceptance

**Phi-Psi angles (Ramachandran plot):**
- ≥90% of residues in favored regions
- <2% in disallowed regions
- Glycine and proline may have different distributions

**Steric clashes:**
- Minimal atom overlaps
- No severe van der Waals violations
- Side chains properly oriented

**Rotamer quality:**
- >90% residues in favored rotameric states
- Proper chi angle distributions
- Reasonable sidechain-backbone contacts

**Overall quality score:**
- DOPE score comparable to experimental structures of similar size
- Consistent across multiple validation tools
- No major geometry deviations from template