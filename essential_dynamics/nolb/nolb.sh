#!/bin/bash

################################################################################
# NORMAL MODE ANALYSIS USING NOLB (Non-Linear Normal Mode Analysis)
################################################################################
#
# Purpose:
#   This script performs Normal Mode Analysis (NMA) using NOLB.
# NOLB efficiently predicts protein dynamics and collective motions
#   by building and analyzing the protein's elastic network.
#
# Prerequisites:
#   - NOLB executable installed and accessible
#   - Input PDB structure file
#   - Sufficient disk space for output trajectory files
#
# Input Parameters:
#   - <input_structure>: PDB structure file
#   - <output_path>: Directory for output files
#   - <cutoff>: Distance cutoff for elastic network (typically 10-15 Å)
#   - <sampling>: Frequency of mode sampling (e.g., 100)
#
# Output Files Generated:
#   - Normal mode trajectories (PDB format)
#   - Eigenvalue spectrum (energetics of modes)
#   - Mode analysis reports (collectivity, amplitude)
#   - Logfile with analysis results
#
################################################################################

# ===== CONFIGURATION PARAMETERS =====
# Set these variables before running the script

# Path to NOLB executable (install NOLB or set path accordingly)
NOLB_PATH="<path to NOLB>"

# Input protein structure (PDB format)
INPUT_STRUCTURE="<input structure>"

# Output directory for results
OUTPUT_PATH="<output path>"

# Distance cutoff for elastic network model (in Angstroms)
# Typical values: 10-15 Å for all-atom, 13-15 Å for C-alpha only
CUTOFF="<cutoff>"

# Sampling frequency for mode trajectory generation
# Higher values increase output file size
SAMPLING="100"

# Output logfile for capturing results and diagnostic information
LOG_FILE="<log file>"


# ===== VALIDATION =====
# Check that required files and directories exist/can be created
cat "Validating NOLB setup...\n"

# Verify NOLB executable exists
if [ ! -f "$NOLB_PATH" ]; then
    echo "Error: NOLB executable not found at: $NOLB_PATH"
    echo "Please install NOLB or update NOLB_PATH variable"
    exit 1
fi

# Verify input structure exists
if [ ! -f "$INPUT_STRUCTURE" ]; then
    echo "Error: Input structure not found: $INPUT_STRUCTURE"
    exit 1
fi

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_PATH" || { echo "Error: Cannot create output directory"; exit 1; }

echo "Setup validation complete."
echo ""


# ===== RUN NOLB NORMAL MODE ANALYSIS =====
# Execute NOLB with specified parameters for elastic network model construction
# and normal mode analysis

cat "Starting NOLB Normal Mode Analysis...\n"
cat "  Input structure: $INPUT_STRUCTURE\n"
cat "  Output directory: $OUTPUT_PATH\n"
cat "  Distance cutoff: $CUTOFF Å\n"
cat "  Sampling frequency: $SAMPLING\n"
cat "  Log file: $LOG_FILE\n\n"

# Execute NOLB with parameters:
#   $NOLB_PATH: Path to NOLB executable
#   $INPUT_STRUCTURE: Input PDB file to analyze
#   -o: Output directory for results
#   -c: Distance cutoff for elastic network (controls network density)
#   -s: Sampling frequency for trajectory output
#   --analyze: Enable mode analysis (collectivity, etc.)
#   --collectivity: Calculate mode collectivity indices (measure of mode delocalization)
#   > $LOG_FILE: Redirect output to logfile
#
# Collectivity measures how many atoms participate in a given mode:
#   - High collectivity: Global motion (many atoms involved)
#   - Low collectivity: Local motion (few atoms involved)

"$NOLB_PATH" "$INPUT_STRUCTURE" \
    -o "$OUTPUT_PATH" \
    -c "$CUTOFF" \
    -s "$SAMPLING" \
    --analyze \
    --collectivity > "$LOG_FILE"

# Capture exit status for error handling
EXIT_STATUS=$?


# ===== RESULTS SUMMARY =====
# Report completion status and provide guidance for next steps

if [ $EXIT_STATUS -eq 0 ]; then
    echo ""
    echo "==============================================="
    echo "NOLB Analysis Complete - Success"
    echo "==============================================="
    echo ""
    echo "Output directory: $OUTPUT_PATH"
    echo "Contents should include:"
    echo "  - Mode trajectory files (PDB format)"
    echo "  - Eigenvalue spectrum"
    echo "  - Collectivity analysis"
    echo "  - Analysis report"
    echo ""
    echo "Next steps:"
    echo "  1. Review logfile: $LOG_FILE"
    echo "  2. Examine collectivity to identify biologically relevant modes"
    echo "  3. Visualize trajectories in PyMOL, VMD, or similar software"
    echo "  4. Compare modes across multiple structures if available"
    echo ""
else
    echo ""
    echo "==============================================="
    echo "Error: NOLB Analysis Failed"
    echo "==============================================="
    echo "Exit status: $EXIT_STATUS"
    echo "Check logfile for details: $LOG_FILE"
    exit $EXIT_STATUS
fi
