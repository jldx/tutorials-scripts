# --- Input files ---
topology='../../top.pdb'   # Topology file (.pdb)
trajectory='../../traj.xtc'       # Trajectory file (.xtc)

# --- Output settings ---
output_prefix='prefix'           # Prefix for all final output files
script_file='secstruct.cpptraj'     # Temporary cpptraj input script

# ============================================================
# Build and run the cpptraj script
# ============================================================
cat > "${script_file}" <<EOF
parm ${topology}
trajin ${trajectory}

secstruct out ${output_prefix}.secstruct.dat sumout ${output_prefix}.secstruct.sumout ${output_prefix}.secstruct.assignout totalout ${output_prefix}.secstruct.totalout

run
quit
EOF

cpptraj -i "${script_file}"
