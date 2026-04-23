#!/usr/bin/env python3
"""
generate_equilibration_mdps.py (suggested explicit name)
---------------------------------------------------
Small helper script to generate a typical set of GROMACS MDP
files used during membrane system preparation and equilibration.

What this script generates:
 - mini.mdp             : steepest-descent energy minimization
 - npt_1.mdp..npt_N.mdp : staged NPT equilibration steps with decreasing
                                                    position restraints
 - npt_canonical.mdp    : final NPT (Parrinello-Rahman + Nose-Hoover)

Inputs / customization:
 - `tc_grps`, `comm_grps`: thermostat and COM removal group names
 - `ref_t`: target temperature (K)
 - `ref_p`: target pressure (bar)
 - `membrane`: True/False to enable semi-isotropic coupling and lipid restraints

Usage (suggested):
        python generate_equilibration_mdps.py

Notes:
 - This file currently uses simple top-level variables for configuration.
     For production use you may want to add an argparse CLI to set options
     at runtime (temperature, pressure, group names, output directory, ...).
"""

# ---------------------------------------------------------------------------
# User-configurable parameters (edit here or extend script to accept CLI args)
# ---------------------------------------------------------------------------
# Thermostat coupling groups (example for membrane systems)
tc_grps = "SOLU MEMB SOLV"
# Groups used for center-of-mass removal
comm_grps = "SOLU_MEMB SOLV"
# Does the system contain a membrane (affects pressure coupling and restraints)
membrane = True
# Target thermodynamic state
ref_t = 300
ref_p = 1

# Default coupling/time constants (can be tuned)
tau_t = 0.1
pcoupltype = "isotropic"
# Position-restraint schedules (kJ/mol/nm^2) for solute and membrane
restraints_mol = [2000.0, 1000.0, 500.0, 25.0, 50.0, 0.0]
restraints_membrane = [500.0, 250.0, 100.0, 50.0, 0.0]

if membrane:
        # For membranes we usually use semi-isotropic pressure coupling
        pcoupltype = "semiisotropic"


################################################################################
# Generate minimization mdp file
################################################################################

# Define restraints
restraints = f"-DPOSRES -DPOSRES_FC_BB={restraints_mol[0]}"

if membrane:
    restraints += f" -DPOSRES_FC_LIPID={restraints_membrane[0]} -DDIHRES -DDIHRES_FC={restraints_membrane[0]}"

# Write file
with open('mini.mdp', 'w') as fout:
    fout.write(f"""define                  = {restraints}
; Steepest Descent Energy Minimization
; Purpose: Remove steric clashes and geometry artifacts before equilibration
integrator              = steep                 ; steepest descent energy minimization
emtol                   = 100.0                 ; stopping criterion: convergence when F_max < emtol (kJ/mol/nm)
emstep                  = 0.01                  ; initial step size (nm) - conservative for stability
nsteps                  = 50000                 ; maximum number of minimization steps
; Note: Minimization typically converges before reaching nsteps
; Output control
nstxtcout               = 100                   ; save coordinates every 100 steps (for minimization trajectory)
nstvout                 = 0                     ; do NOT save velocities (no velocities in minimization)
nstfout                 = 0                     ; do NOT save forces (not relevant for minimization)
nstcalcenergy           = 10                    ; calculate energies every 10 steps
nstenergy               = 10                    ; save energies every 10 steps (monitor minimization progress)
nstlog                  = 10                    ; update log file every 10 steps (monitor convergence)
; Periodic boundary conditions
pbc                     = xyz                   ; 3D PDB
; Bond parameters
constraints             = all-bonds             ; constrain all bonds (including hydrogen)
constraint_algorithm    = LINCS                 ; Linear Constraint Solver (accurate and fast)
lincs_iter              = 1                     ; accuracy of LINCS (1 is sufficient for most systems)
lincs_order             = 4                     ; related to accuracy (4 is default and recommended)
; Neighbor searching (Verlet cutoff scheme)
cutoff-scheme           = Verlet                ; modern Verlet neighbor list scheme (recommended)
ns_type                 = grid                  ; search neighboring grid cells (efficient with Verlet)
nstlist                 = 20                    ; update neighbor list every 20 steps (40 fs) - good balance
rlist                   = 1.0                   ; cut-off distance for neighbor list (nm) - matches rcoulomb
; Electrostatics
rcoulomb                = 1.0                   ; short-range electrostatic cutoff (nm) - increased to 1.0 nm for better accuracy
coulombtype             = PME                   ; Particle Mesh Ewald for long-range electrostatics (gold standard)
pme_order               = 4                     ; cubic interpolation (default: good accuracy-speed balance)
fourierspacing          = 0.12                  ; grid spacing for FFT (nm) - reduced to 0.12 for better accuracy
; Van der Waals
rvdw                    = 1.0                   ; short-range van der waals cutoff (nm) - matches rcoulomb
vdwtype                 = Cut-off               ; plain cutoff with pair list radius rlist and rvdw
vdw-modifier            = Force-switch          ; smooth switch forces to zero at rvdw_switch
rvdw_switch             = 0.8                   ; where force switch begins (nm) - allows smooth transition
DispCorr                = EnerPres              ; account for cut-off vdW scheme (energy and pressure correction)
; Temperature coupling is NOT used for energy minimization
; Pressure coupling is NOT used for energy minimization
; Coordinates scaling and center of mass motion removal
refcoord_scaling        = com                   ; scale coordinates relative to center of mass
nstcomm                 = 1                     ; do NOT remove center of mass motion (minimization doesn't need it)
comm_mode               = linear                ; remove linear center of mass velocity
comm_grps               = {comm_grps}         ; groups to remove the center of mass           
""")

################################################################################
# Generate npt files for intermediate equilibration steps
################################################################################

tcoupl = "v-rescale"
pcoupl = "c-rescale"

temp_par = f"""tc_grps                 = {tc_grps}        ; three separate temperature coupling groups
tau_t                   = 0.1  0.1         ; time constant (ps) - reduced from 0.5 for better equilibration
ref_t                   = {ref_t}  {ref_t}         ; reference temperature (K) for each group
"""

if membrane:
        temp_par = f"""tc_grps                 = {tc_grps}        ; three separate temperature coupling groups
tau_t                   = 0.1  0.1  0.1         ; time constant (ps) - reduced from 0.5 for better equilibration
ref_t                   = {ref_t}  {ref_t}  {ref_t}         ; reference temperature (K) for each group
"""

for i in range(1, len(restraints_mol) - 1):

    # define restaints and number of time we have to put the temperature
    restraints = f"-DPOSRES -DPOSRES_FC_BB={restraints_mol[i]}"

    if membrane:
        restraints += f" -DPOSRES_FC_LIPID={restraints_membrane[i]} -DDIHRES -DDIHRES_FC={restraints_membrane[i]}"

    with open(f"npt_{i}.mdp", "w") as fout:
        # Write file
        fout.write(f"""define                  = {restraints}
; Run parameters
integrator              = md                    ; leap-frog integrator (efficient for equilibration)
dt                      = 0.002                 ; 2 fs timestep
nsteps                  = 500000                ; 1000 ps (1 ns) total simulation time
continuation            = yes                   ; continuing from prior simulation
; Output control
nstxtcout               = 5000                  ; save coordinates to XTC every 10 ps (compressed trajectory)
nstvout                 = 5000                  ; save velocities every 10 ps
nstfout                 = 5000                  ; save forces every 10 ps (disabled for production, kept for debugging)
nstcalcenergy           = 1000                  ; calculate energies every 2 ps (more frequent for monitoring)
nstenergy               = 1000                  ; save energies every 2 ps (for better energy analysis)
nstlog                  = 1000                  ; update log file every 2 ps (more frequent monitoring)
; Periodic boundary conditions
pbc                     = xyz                   ; 3D PDB
; Bond parameters
constraints             = all-bonds             ; constrain all bonds (including hydrogen)
constraint_algorithm    = LINCS                 ; Linear Constraint Solver (accurate and fast)
lincs_iter              = 1                     ; accuracy of LINCS (1 is sufficient for most systems)
lincs_order             = 4                     ; related to accuracy (4 is default and recommended)
; Neighbor searching (Verlet cutoff scheme)
cutoff-scheme           = Verlet                ; modern Verlet neighbor list scheme (recommended)
ns_type                 = grid                  ; search neighboring grid cells (efficient with Verlet)
nstlist                 = 20                    ; update neighbor list every 20 steps (40 fs) - good balance
rlist                   = 1.0                   ; cut-off distance for neighbor list (nm) - matches rcoulomb
; Electrostatics
rcoulomb                = 1.0                   ; short-range electrostatic cutoff (nm) - increased to 1.0 nm for better accuracy
coulombtype             = PME                   ; Particle Mesh Ewald for long-range electrostatics (gold standard)
pme_order               = 4                     ; cubic interpolation (default: good accuracy-speed balance)
fourierspacing          = 0.12                  ; grid spacing for FFT (nm) - reduced to 0.12 for better accuracy
; Van der Waals
rvdw                    = 1.0                   ; short-range van der waals cutoff (nm) - matches rcoulomb
vdwtype                 = Cut-off               ; plain cutoff with pair list radius rlist and rvdw
vdw-modifier            = Force-switch          ; smooth switch forces to zero at rvdw_switch
rvdw_switch             = 0.8                   ; where force switch begins (nm) - allows smooth transition
DispCorr                = EnerPres              ; account for cut-off vdW scheme (energy and pressure correction)
; Temperature coupling is on
; Separate thermostats for groups for better equilibration
tcoupl                  = {tcoupl}             ; Temperature coupling
nsttcouple              = 10                    ; apply thermostat every 10 steps (default is 1)
{temp_par}
; Pressure coupling is on
pcoupl                  = {pcoupl}              ; Pessure coupling
pcoupltype              = {pcoupltype}          ; semiisotropic: scale x-y (membrane) and z (solvent) separately
nstpcouple              = 5                     ; apply pressure coupling every 5 steps
tau_p                   = 1.0                   ; time constant (ps) - increased from 5.0 for stability during equilibration
ref_p                   = 1.0   1.0             ; reference pressure (bar) - 1 atm for both x-y and z
compressibility         = 4.5e-5    4.5e-5     ; isothermal compressibility (bar^-1) - from GROMACS manual for water
; Coordinates scaling and center of mass motion removal
refcoord_scaling        = com                   ; scale coordinates relative to center of mass
nstcomm                 = 100                   ; remove center of mass motion every 0.2 ps (every 100 steps)
comm_mode               = linear                ; remove linear center of mass velocity
comm_grps               = {comm_grps}        ; remove COM motion separately for (solute+membrane) and solvent
""")
        
################################################################################
# Generate npt file for canonical npt equilibration steps
################################################################################

tcoupl = "nose-hoover"
pcoupl = "parrinello-rahman"

with open("npt_canonical.mdp", "w") as fout:
    # Write file
    fout.write(f"""; Run parameters
integrator              = md                    ; leap-frog integrator (efficient for equilibration)
dt                      = 0.002                 ; 2 fs timestep
nsteps                  = 500000                ; 1000 ps (1 ns) total simulation time
continuation            = yes                   ; continuing from prior simulation
; Output control
nstxtcout               = 5000                  ; save coordinates to XTC every 10 ps (compressed trajectory)
nstvout                 = 5000                  ; save velocities every 10 ps
nstfout                 = 5000                  ; save forces every 10 ps (disabled for production, kept for debugging)
nstcalcenergy           = 1000                  ; calculate energies every 2 ps (more frequent for monitoring)
nstenergy               = 1000                  ; save energies every 2 ps (for better energy analysis)
nstlog                  = 1000                  ; update log file every 2 ps (more frequent monitoring)
; Periodic boundary conditions
pbc                     = xyz                   ; 3D PDB
; Bond parameters
constraints             = all-bonds             ; constrain all bonds (including hydrogen)
constraint_algorithm    = LINCS                 ; Linear Constraint Solver (accurate and fast)
lincs_iter              = 1                     ; accuracy of LINCS (1 is sufficient for most systems)
lincs_order             = 4                     ; related to accuracy (4 is default and recommended)
; Neighbor searching (Verlet cutoff scheme)
cutoff-scheme           = Verlet                ; modern Verlet neighbor list scheme (recommended)
ns_type                 = grid                  ; search neighboring grid cells (efficient with Verlet)
nstlist                 = 20                    ; update neighbor list every 20 steps (40 fs) - good balance
rlist                   = 1.0                   ; cut-off distance for neighbor list (nm) - matches rcoulomb
; Electrostatics
rcoulomb                = 1.0                   ; short-range electrostatic cutoff (nm) - increased to 1.0 nm for better accuracy
coulombtype             = PME                   ; Particle Mesh Ewald for long-range electrostatics (gold standard)
pme_order               = 4                     ; cubic interpolation (default: good accuracy-speed balance)
fourierspacing          = 0.12                  ; grid spacing for FFT (nm) - reduced to 0.12 for better accuracy
; Van der Waals
rvdw                    = 1.0                   ; short-range van der waals cutoff (nm) - matches rcoulomb
vdwtype                 = Cut-off               ; plain cutoff with pair list radius rlist and rvdw
vdw-modifier            = Force-switch          ; smooth switch forces to zero at rvdw_switch
rvdw_switch             = 0.8                   ; where force switch begins (nm) - allows smooth transition
DispCorr                = EnerPres              ; account for cut-off vdW scheme (energy and pressure correction)
; Temperature coupling is on
; Separate thermostats for groups for better equilibration
tcoupl                  = {tcoupl}             ; Temperature coupling
nsttcouple              = 10                    ; apply thermostat every 10 steps (default is 1)
{temp_par}
; Pressure coupling is on
pcoupl                  = {pcoupl}              ; Pessure coupling
pcoupltype              = {pcoupltype}          ; semiisotropic: scale x-y (membrane) and z (solvent) separately
nstpcouple              = 5                     ; apply pressure coupling every 5 steps
tau_p                   = 1.0                   ; time constant (ps) - increased from 5.0 for stability during equilibration
ref_p                   = 1.0   1.0             ; reference pressure (bar) - 1 atm for both x-y and z
compressibility         = 4.5e-5    4.5e-5     ; isothermal compressibility (bar^-1) - from GROMACS manual for water
; Coordinates scaling and center of mass motion removal
refcoord_scaling        = com                   ; scale coordinates relative to center of mass
nstcomm                 = 100                   ; remove center of mass motion every 0.2 ps (every 100 steps)
comm_mode               = linear                ; remove linear center of mass velocity
comm_grps               = {comm_grps}        ; remove COM motion separately for (solute+membrane) and solvent
""")
    