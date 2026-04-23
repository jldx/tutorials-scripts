# Umbrella Sampling with GROMACS: Free Energy Profile Tutorial

## Overview

This tutorial demonstrates how to calculate the **free energy profile (PMF)** for pulling a small molecule (CLM) through a lipid membrane bilayer (POPC) using **Umbrella Sampling** (US) combined with **Weighted Histogram Analysis Method** (WHAM).

**System:** CLM molecule + POPC membrane + water

**Collective Variable (CV):** Z-coordinate (perpendicular to membrane plane)

**Output:** Free Energy Profile showing energy barriers to membrane crossing

---

## Background

### What is Umbrella Sampling?

Umbrella Sampling is a free energy calculation technique that:
1. **Divides the reaction pathway** into multiple overlapping "windows"
2. **Applies harmonic potentials** (umbrellas) to sample each window thoroughly
3. **Combines results** using WHAM to reconstruct the unbiased free energy profile

### Why Umbrella Sampling?

Standard molecular dynamics rarely samples rare events (e.g., molecule crossing a membrane barrier). Umbrella Sampling biases the system to sample all relevant configurations:
- ✓ Explores high energy barriers
- ✓ Generates overlapping windows for statistical accuracy
- ✓ Computes free energy differences between states
- ✓ Captures entropic and enthalpic contributions

### Mathematical Framework

For each window $i$ with spring constant $k$ centered at position $z_0$:

$$V_i(z) = \frac{k}{2}(z - z_{0,i})^2$$

WHAM combines multiple biased simulations to recover the unbiased PMF:

$$A(z) = -k_B T \ln P(z) + \text{const}$$

where $P(z)$ is the unbiased probability distribution.

---

## Workflow Overview

```
┌─────────────────────────────────────────────────────────────────────┐
│ 1. STEERED MD (SMD)                                                 │
│    Generate initial pulling trajectory to identify CV range         │
│    Output: smd_pullx.xvg (CV vs. time)                              │
└──────────────────────────┬──────────────────────────────────────────┘
                           │
                           ▼
┌─────────────────────────────────────────────────────────────────────┐
│ 2. EXTRACT WINDOWS                                                  │
│    Select overlapping windows from SMD trajectory                   │
│    Extract conformations for each window                            │
│    Output: conf_*.gro (window starting structures)                  │
└──────────────────────────┬──────────────────────────────────────────┘
                           │
                           ▼
┌─────────────────────────────────────────────────────────────────────┐
│ 3. EQUILIBRATION (Optional but recommended)                         │
│    Briefly equilibrate each window structure                        │
│    Output: conf_*_eq.gro (equilibrated structures)                  │
└──────────────────────────┬──────────────────────────────────────────┘
                           │
                           ▼
┌─────────────────────────────────────────────────────────────────────┐
│ 4. UMBRELLA SAMPLING                                                │
│    Run US for each window with harmonic restraint                   │
│    Output: pullx.xvg, pullf.xvg for each window                     │
└──────────────────────────┬──────────────────────────────────────────┘
                           │
                           ▼
┌─────────────────────────────────────────────────────────────────────┐
│ 5. WHAM ANALYSIS                                                    │
│    Combine all windows to compute unbiased PMF                      │
│    Output: histo.dat, profile.xvg                                   │
└─────────────────────────────────────────────────────────────────────┘
```

---

## Detailed Steps

### Step 1: Steered MD (1_smd/)

**Purpose:** Generate a pulling trajectory to understand the CV range and identify barrier heights.

**MDP File:** `smd.mdp`
- Integrator: `md` (leap-frog)
- Duration: 4 ns
- Pull code: umbrella sampling with harmonic potential
- Pulling direction: Z-axis (perpendicular to membrane)
- Pulling velocity: 0.001 nm/ps (slow quasi-static pulling)
- Spring constant: 1000 kJ/(mol·nm²)

**Commands:**

```bash
cd 1_smd/

# Prepare binary trajectory file
gmx grompp -f smd.mdp -o smd.tpr \
    -c <topology> \
    -r <reference structure> \
    -p topol.top \
    -n index.ndx \
    -maxwarn 1

# Run the steered MD simulation
gmx mdrun -v -deffnm smd
```

**Outputs:** (additional to Gromacs outputs for simulation)
- `smd_pullx.xvg`: Reaction coordinate (Z-position) vs. time
- `smd_pullf.xvg`: Pulling force vs. time

**Interpreting the Output:**

The file `smd_pullx.xvg` shows how Z-coordinate changes during pulling:
- **X-axis:** Time (ps)
- **Y-axis:** Z-distance (nm) between solute and reference groups
- **Key observations:**
  - Where does Z increase fastest? (Low friction)
  - Where does Z plateau? (Energy barrier)
  - Plateaus indicate membrane crossing regions

The file `smd_pullf.xvg` shows the applied force:
- **Peaks** in force indicate energy barriers
- **Barrier height** ≈ max force × pulling velocity × (# barriers)

---

### Step 2: Extract Umbrella Windows (2_select_umbrella_windows/)

**Purpose:** Divide the pulling pathway into overlapping windows for umbrella sampling.

**Files Provided:**
- `select_umbrella_windows.ipynb`: Jupyter notebook for window selection

**Key Point:** Windows should overlap by ~0.1–0.2 nm to ensure proper statistical coverage.

---

### Step 3: Umbrella Sampling Runs (3_umbrella_sampling/)

#### 3a. Equilibration (Optional)

**Purpose:** Gently equilibrate each window structure before production US simulation.

**Directory:** `3_umbrella_sampling/a_equilibration/`

**Setup Script:** `set_up_run_umbrella_windows_equilibration.sh`

**Files Provided:**
- `equilibration.mdp` and `run_equilibration.pbs`: initial raw files
- `checking_equilibration.ipynb`: Jupyter notebook to check the windows equilibration before proceeding to umbrella sampling

#### 3b. Umbrella Sampling (Production)

**Purpose:** Run restricted MD in each window to sample the conformational space with an umbrella potential.

**Directory:** `3_run_umbrella/b_umbrella_sampling/`

**Setup Script:** `set_up_run_umbrella_sampling.sh`

**Files Provided:**
- `umbrella.mdp` and `run_umbrella.pbs`: initial raw files

---

### Step 4: WHAM Analysis (4_fep/)

#### 4a. Run the analysis

**Purpose:** Combine all windowed umbrella sampling simulations to compute the unbiased free energy profile.

**Files Needed:**
- `*.tpr`: Binary run files from Step 3b (contain spring constants)
- `*_pullx.xvg`: Pulled coordinate from each window
- `*_pullf.xvg`: Pulling force from each window (optional)

**Script:** `run_wham.sh`

**WHAM Output Files:**

1. **histo.dat**: Histogram of CV sampling
   - Shows which windows sampled which regions
   - Check for gaps (indicate insufficient overlap)

2. **profile.xvg**: Final free energy profile (PMF)
   - Column 1: CV (Z-coordinate, nm)
   - Column 2: Free energy (kJ/mol)
   - Column 3: Error estimate (kJ/mol)

3. **histogram.xvg**: Energy landscape visualization

#### 4b. Visualize Results

1. **Check window overlap:**
   ```python
   import numpy as np
   import matplotlib.pyplot as plt
   
   # Read histogram
   histo = np.loadtxt('histo.dat')
   
   plt.figure(figsize=(12, 6))
   plt.imshow(histo.T, aspect='auto', origin='lower')
   plt.xlabel('CV (nm)')
   plt.ylabel('Window index')
   plt.title('Window sampling coverage (should show overlap)')
   plt.colorbar()
   plt.show()
   ```

2. **Plot the PMF:**
   ```python
   # Read PMF
   pmf = np.loadtxt('profile.xvg')
   cv = pmf[:, 0]
   free_energy = pmf[:, 1]
   error = pmf[:, 2]
   
   plt.figure(figsize=(12, 6))
   plt.plot(cv, free_energy, 'b-', linewidth=2, label='PMF')
   plt.fill_between(cv, free_energy - error, free_energy + error, alpha=0.3, label='Bootstrap error')
   plt.xlabel('Z-coordinate (nm)', fontsize=12)
   plt.ylabel('Free Energy (kJ/mol)', fontsize=12)
   plt.title('Free Energy Profile: CLM crossing POPC Membrane', fontsize=14)
   plt.legend()
   plt.grid(alpha=0.3)
   plt.show()
   ```

3. **Identify barriers:**
   ```python
   # Find energy maxima (barriers)
   barrier_indices = np.argsort(free_energy)[-3:]  # Top 3 barriers
   for idx in barrier_indices:
       print(f"Barrier at Z = {cv[idx]:.2f} nm, ΔG = {free_energy[idx]:.2f} kJ/mol")
   ```

---

## Troubleshooting Guide

### Common Issues and Solutions

#### 1. **Large force fluctuations in pullf.xvg**
**Symptom:** Force oscillates wildly around 0
**Cause:** Spring constant too weak or sampling too fast
**Solution:**
- Increase `pull_coord1_k` (e.g., 1000 → 2000 kJ/(mol·nm²))
- Decrease production run temperature slightly
- Increase simulation time per window

#### 2. **Gaps in WHAM histogram**
**Symptom:** Horizontal white lines in histo.dat (no sampling in some regions)
**Cause:** Windows too far apart or insufficient overlap
**Solution:**
- Add intermediate windows between large gaps
- Reduce window spacing (e.g., 0.2 nm → 0.1 nm)
- Increase simulation time per window

#### 3. **Unrealistic free energy jumps**
**Symptom:** Large discontinuities in PMF between adjacent windows
**Cause:** Poor overlap or systematic bias in windows
**Solution:**
- Verify histogram for continuous coverage
- Re-run umbrella sampling with longer duration
- Check that initial window positions match the target umbrellas

#### 4. **WHAM convergence issues**
**Symptom:** "Warning: WHAM did not converge"
**Cause:** Insufficient overlap or too many windows
**Solution:**
- Reduce number of bootstrap samples (`-nBootstrap 50`)
- Increase window overlap
- Re-extract windows with finer spacing

#### 5. **Memory issues with WHAM**
**Symptom:** WHAM crashes with segmentation fault
**Cause:** Too many windows or too long trajectories
**Solution:**
- Skip every Nth frame: `gmx wham ... -skipPoints 10`
- Reduce bootstrap samples
- Run on high-memory node

---

## Best Practices

### Simulation Setup
✓ **Always equilibrate the system** before SMD/US (NVT → NPT)
✓ **Use semiisotropic pressure coupling** for membrane systems (x-y vs. z)
✓ **Define clear groups** in index.ndx (SOLU, CLM, MEMB, SOLV)
✓ **Save forces frequently** in SMD (needed for WHAM analysis)

### Window Selection
✓ **Aim for 10–20% overlap** between adjacent windows
✓ **Use uniform spacing** across the CV range
✓ **Ensure first and last windows** sample stable regions (force ≈ 0)
✓ **Test with coarse spacing first**, then refine

### Umbrella Sampling
✓ **Use weak temperature coupling** (tau_t = 0.5–1.0 ps)
✓ **Run each window ≥ 1 ns** for good statistics (≥ 2 ns recommended)
✓ **Monitor pullx.xvg** during production (should show fluctuations around umbrella center)
✓ **Keep spring constant consistent** across windows

### WHAM Analysis
✓ **Check histogram coverage** before trusting PMF
✓ **Use bootstrap errors** to assess uncertainty
✓ **Compare with alternative methods** (e.g., Jarzynski equality from SMD)
✓ **Document convergence criteria** used

---

## Advanced Topics

### Calculating Work from SMD

From `steer_md_pullf.xvg`, estimate free energy change:

$$W_{\text{SMD}} \approx \int_0^T F(t) \, dz = \int_0^T F(t) \cdot v \, dt$$

where $v = 0.001$ nm/ps is the pulling velocity.

```bash
# Integrate pulling force to get work done
gmx analyze -f steer_md_pullf.xvg -integrate
```

Compare SMD-derived work with umbrella sampling results (Jarzynski equality).

### Pulling Rate Effects

The relationship between pulling velocity and free energy:

$$\Delta G_{\infty} = G_0 - \frac{k_B T}{v_f} \ln \left\langle e^{-\frac{W(v_f)}{k_B T}} \right\rangle$$

Running SMD at **multiple velocities** and extrapolating to $v \to 0$ gives true equilibrium PMF.

### Multi-dimensional CVs

For more complex processes (e.g., protein unfolding + translocation):
- Define multiple pull coordinates
- Use 2D or 3D WHAM
- Increase computational cost accordingly

---

## Example Data Interpretation

Assuming the PMF calculation yields:

| Z-coordinate (nm) | ΔG (kJ/mol) | Structure |
|---|---|---|
| 0.0 | 0 | Fully inside membrane (hydrophobic core) |
| 0.5 | +45 | Partial insertion through headgroup |
| 1.0 | +65 | Highest barrier (membrane crossing) |
| 1.5 | +40 | Exiting on other side |
| 2.0 | +5 | At membrane interface |
| 3.0 | 0 | Fully in bulk water |

**Interpretation:**
- Energy barrier ≈ 65 kJ/mol (required to pull molecule through membrane)
- Favorable energy in bulk water (ΔG = 0)
- Unfavorable energy inside membrane (ΔG > 0)
- Roughly **asymmetric profile** (different barriers entering vs. exiting)
