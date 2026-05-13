#!python3

import pandas as pd
import numpy as np
import joblib
import sys

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

    return dt_xyz


# ------------------------------------------------------------------------------

perm_file = sys.argv[1]
diff_file = sys.argv[2]

perm = joblib.load(perm_file)
n_mols = perm.shape[0]
x_mols, y_mols, z_mols  = perm['permx'], perm['permy'], perm['permz']

# ------------------------------------------------------------------------------

displ_xyz_mols_lst = []
# displ_avg_xyz_mols_lst = []

for mol_num in range(n_mols):
        
    # Convert mol coordinates to array n x 3
    coords = list(zip(x_mols[mol_num], y_mols[mol_num], z_mols[mol_num]))
    coords_mol = np.array(coords)

    # Calculate displacement
    displ_xyz = calculate_average_msd_mol_nd(coords_mol)
    displ_xyz['mol_num'] = mol_num
    displ_xyz_mols_lst.append(displ_xyz)
    
displ_xyz_mols = pd.concat(displ_xyz_mols_lst, axis=0).reset_index(drop=True)

joblib.dump(displ_xyz_mols, diff_file)


