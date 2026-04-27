#!/bin/bash
# ============================================================
# Hydrogen bond analysis using cpptraj
# Computes HBonds between two residue ranges (both directions)
# and merges the output into single files.
# ============================================================

# --- Input files ---
topology='../../top.pdb'   # Topology file (.pdb)
trajectory='../../traj.xtc'       # Trajectory file (.xtc)

# --- Residue ranges (AMBER selection format) ---
partner_1='1-10'    # First  partner (donor & acceptor)
partner_2='10-100'  # Second partner (donor & acceptor)

# --- Output settings ---
output_prefix='prefix'           # Prefix for all final output files
script_file='hbonds.cpptraj'     # Temporary cpptraj input script

# ============================================================
# Build and run the cpptraj script
# HB1: partner_1 as acceptor, partner_2 as donor
# HB2: partner_2 as acceptor, partner_1 as donor
# Both directions are needed for a complete HBond picture.
# ============================================================
cat > "${script_file}" <<EOF
parm ${topology}
trajin ${trajectory}
hbond HB1 out tmp1.hbonds.numbers.dat \
    acceptormask :${partner_1}@N=,O= donormask :${partner_2}@N=,O= \
    dist 3.6 angle 150 \
    avgout tmp1.hbonds.avgout.dat \
    series uuseries tmp1.hbonds.series.npy
hbond HB2 out tmp2.hbonds.numbers.dat \
    acceptormask :${partner_2}@N=,O= donormask :${partner_1}@N=,O= \
    dist 3.6 angle 150 \
    avgout tmp2.hbonds.avgout.dat \
    series uuseries tmp2.hbonds.series.npy
run
quit
EOF

cpptraj -i "${script_file}"

# ============================================================
# Merge output files
# ============================================================

# avgout: plain row concatenation, skip header of tmp2
cat tmp1.hbonds.avgout.dat \
    <(tail -n +2 tmp2.hbonds.avgout.dat) \
    > "${output_prefix}.hbonds.avgout.dat"

# numbers: paste side by side, keep only last column of tmp2 and sum up the values
# (avoids duplicating the frame index column)
paste tmp1.hbonds.numbers.dat \
    <(awk '{print $NF}' tmp2.hbonds.numbers.dat) \
    > "tmp3.hbonds.numbers.dat"

awk 'NR==1 {print $1, "HB_total"} NR>1 {print $1, $2+$3}' tmp3.hbonds.numbers.dat > ${output_prefix}.hbonds.numbers.dat

# series: paste side by side, drop first column of tmp2
# (avoids duplicating the frame index column)
paste tmp1.hbonds.series.npy \
    <(awk '{$1=""; sub(/^ +/, ""); print}' tmp2.hbonds.series.npy) \
    > "${output_prefix}.hbonds.series.npy"

# ============================================================
# Cleanup temporary files
# ============================================================
rm "${script_file}" \
   tmp1.hbonds.avgout.dat tmp2.hbonds.avgout.dat \
   tmp1.hbonds.numbers.dat tmp2.hbonds.numbers.dat tmp3.hbonds.numbers.dat \
   tmp1.hbonds.series.npy tmp2.hbonds.series.npy

