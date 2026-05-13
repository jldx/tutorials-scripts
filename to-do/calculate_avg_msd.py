#! python3
import joblib
import pandas as pd
import MDAnalysis as mda
import numpy as np
import time


def calculate_average_msd_mol_nd(coords):
    
    # Calculate displacement for each tau
    displ = np.array([0.0])
    taus = np.array(0)
    n_frames = coords.shape[0]

    for i in range(0, n_frames, 10):
        tau = 0
        for j in range(i, n_frames - 1):         
            dist = np.linalg.norm(coords[i, :] - coords[j, :])
            taus = np.insert(taus, -1, tau)
            displ = np.insert(displ, -1, dist)
            tau += 1

    tau, displ = taus[:-1], displ[:-1]
    
    dt_xyz = pd.DataFrame({'tau': tau, 'displ': displ, 'squared_displ': displ**2})
    dt_xyz_avg = dt_xyz.groupby('tau').mean().reset_index()

    return dt_xyz, dt_xyz_avg

# ------------------------------------------------------------------------------

top_file = '/data/ledoux/awc-water-chloride-box/water_ions.gro'
traj_file = '/data/ledoux/awc-water-chloride-box/water_ions_25ns_nojump.xtc'
u = mda.Universe(top_file, traj_file)
water = u.select_atoms('name OW')
n_mols = water.atoms.n_atoms
n_frames = u.trajectory.n_frames

# Get the coordinates of all molecules
# Will be replaced by the reading of the joblib file from the permeation analysis
x_mols = np.zeros((n_mols, n_frames))
y_mols = np.zeros((n_mols, n_frames))
z_mols = np.zeros((n_mols, n_frames))

for i, frame in enumerate(u.trajectory):
    x_mols[:, i] = water.positions[:, 0]
    y_mols[:, i] = water.positions[:, 1]
    z_mols[:, i] = water.positions[:, 2]


displ_xyz_mols_lst = []
displ_avg_xyz_mols_lst = []


start_time = time.time() # Processing time monitoring

for mol_num in range(n_mols):
    
    if mol_num % 10 == 0:
        print(f'{mol_num} mols / {n_mols} - Processed in {(time.time() - start_time) // 60} min')
        start_time = time.time() # Processing time monitoring
        
    # Convert mol coordinates to array n x 3
    coords = list(zip(x_mols[mol_num], y_mols[mol_num], z_mols[mol_num]))
    coords_mol = np.array(coords)

    # Calculate displacement and average displacement
    displ_xyz, displ_avg_xyz = calculate_average_msd_mol_nd(coords_mol)
    displ_xyz['mol_num'] = mol_num
    displ_avg_xyz['mol_num'] = mol_num
    
    displ_xyz_mols_lst.append(displ_xyz)
    displ_avg_xyz_mols_lst.append(displ_avg_xyz)
    
displ_xyz_mols = pd.concat(displ_xyz_mols_lst, axis=0).reset_index(drop=True)
displ_avg_xyz_mols = pd.concat(displ_avg_xyz_mols_lst, axis=0).reset_index(drop=True)

joblib.dump([displ_xyz_mols, displ_avg_xyz_mols], 'displacement_xyz.joblib')