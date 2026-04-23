#!python3


import numpy as np
import pandas as pd
import multiprocessing
import MDAnalysis as mda
from MDAnalysis.analysis.hydrogenbonds import HydrogenBondAnalysis

import warnings
warnings.filterwarnings('ignore')

def hbonds_analysis(universe, start, stop, step, donors_sel_str, hydrogens_sel_str,
                    acceptors_sel_str, d_a_cutoff, d_h_a_angle_cutoff):
    """
    Performs hydrogen bond analysis on a given range of frames.

    Parameters:
    universe (MDAnalysis Universe): The MDAnalysis Universe object containing the molecular structure.
    start (int): The starting frame for the analysis.
    stop (int): The stopping frame for the analysis.
    step (int): The step size for frame sampling.
    donors_sel_str (str): Selection string for donor atoms.
    hydrogens_sel_str (str): Selection string for hydrogen atoms.
    acceptors_sel_str (str): Selection string for acceptor atoms.
    d_a_cutoff (float): Distance cutoff between donor and acceptor atoms for hydrogen bond formation.
    d_h_a_angle_cutoff (float): Angle cutoff for donor-hydrogen-acceptor angle for hydrogen bond formation.

    Returns:
    numpy.ndarray: Array of hydrogen bond results with columns ['frame', 'donor_index', 'hydrogen_index',
                                                                  'acceptor_index', 'DA_distance', 'DHA_angle'].
    """
    hbonds = HydrogenBondAnalysis(universe=universe,
                                  donors_sel=donors_sel_str,
                                  hydrogens_sel=hydrogens_sel_str,
                                  acceptors_sel=acceptors_sel_str,
                                  d_a_cutoff=d_a_cutoff,
                                  d_h_a_angle_cutoff=d_h_a_angle_cutoff,
                                  update_selections=False)
    hbonds.run(start=start, stop=stop, step=step, verbose=False)
    return hbonds


def hbonds_analysis_worker(task):
    """
    Worker function to unpack task tuple and call hbonds_analysis.

    Parameters:
    task (tuple): A tuple containing the parameters for hbonds_analysis.

    Returns:
    numpy.ndarray: Array of hydrogen bond results from the hbonds_analysis function.
    """
    universe, start, stop, step, donors_sel_str, hydrogens_sel_str, \
        acceptors_sel_str, d_a_cutoff, d_h_a_angle_cutoff = task
    return hbonds_analysis(universe, start, stop, step, donors_sel_str, hydrogens_sel_str,
                           acceptors_sel_str, d_a_cutoff, d_h_a_angle_cutoff)


def calculate_hbonds(universe, step, donors_sel_str, hydrogens_sel_str,
                     acceptors_sel_str, d_a_cutoff, d_h_a_angle_cutoff, num_processes):
    """
    Calculates hydrogen bonds in parallel using multiple processes.

    Parameters:
    universe (MDAnalysis Universe): The MDAnalysis Universe object containing the molecular structure.
    step (int): The step size for frame sampling.
    donors_sel_str (str): Selection string for donor atoms.
    hydrogens_sel_str (str): Selection string for hydrogen atoms.
    acceptors_sel_str (str): Selection string for acceptor atoms.
    d_a_cutoff (float): Distance cutoff between donor and acceptor atoms for hydrogen bond formation.
    d_h_a_angle_cutoff (float): Angle cutoff for donor-hydrogen-acceptor angle for hydrogen bond formation.
    num_processes (int): Number of parallel processes to use for the analysis.

    Returns:
    numpy.ndarray: Combined array of hydrogen bond results from all processes.
    """
    num_frames = len(universe.trajectory)
    chunk_size = num_frames // num_processes
    tasks = [(universe, i * chunk_size, min((i + 1) * chunk_size, num_frames), step, donors_sel_str, \
              hydrogens_sel_str, acceptors_sel_str, d_a_cutoff, d_h_a_angle_cutoff) for i in range(num_processes)]

    with multiprocessing.Pool(processes=num_processes) as pool:
        results = pool.map(hbonds_analysis_worker, tasks)

    return results
    
    
def add_atoms_residues_details(universe, df):
    """
    Adds detailed information about atoms and residues to a DataFrame of hydrogen bonds.

    Parameters:
    universe (MDAnalysis.Universe): The MDAnalysis Universe object containing the molecular structure.
    df (pandas.DataFrame): DataFrame containing hydrogen bond information with columns:
        - 'donor_index': Index of the donor atom.
        - 'hydrogen_index': Index of the hydrogen atom.
        - 'acceptor_index': Index of the acceptor atom.

    Returns:
    pandas.DataFrame: The input DataFrame with additional columns:
        - 'donor_atom': Name of the donor atom.
        - 'hydrogen_atom': Name of the hydrogen atom.
        - 'acceptor_atom': Name of the acceptor atom.
        - 'donor_elem': Element symbol of the donor atom (first letter of the atom name).
        - 'acceptor_elem': Element symbol of the acceptor atom (first letter of the atom name).
        - 'hbond_detail': Detailed representation of the hydrogen bond as D-H···A.
        - 'hbond_detail_elem': Element-level representation of the hydrogen bond as D-H···A.

    Example:
    >>> import MDAnalysis as mda
    >>> import pandas as pd
    >>> u = mda.Universe('topology.psf', 'trajectory.dcd')
    >>> data = {
    ...     'donor_index': [1, 2],
    ...     'hydrogen_index': [3, 4],
    ...     'acceptor_index': [5, 6]
    ... }
    >>> hbonds_df = pd.DataFrame(data)
    >>> add_atoms_residues_details(u, hbonds_df)
    """
    df['donor_atom'] = universe.atoms.names[df['donor_index']]
    df['hydrogen_atom'] = universe.atoms.names[df['hydrogen_index']]
    df['acceptor_atom'] = universe.atoms.names[df['acceptor_index']]
    df['donor_elem'] = [name[0] for name in df['donor_atom']]
    df['acceptor_elem'] = [name[0] for name in df['acceptor_atom']]
    df['hbond_detail'] = [f'{d}-H···{a}' for d, a in zip(df['donor_atom'], df['acceptor_atom'])]
    df['hbond_detail_elem'] = [f'{d}-H···{a}' for d, a in zip(df['donor_elem'], df['acceptor_elem'])]

    return df



def create_hbonds_df(universe, hbonds):
    """
    Creates a pandas DataFrame with detailed hydrogen bond information from the analysis results.

    Parameters:
    universe (MDAnalysis Universe): The MDAnalysis Universe object containing the molecular structure.
    hbonds (numpy.ndarray): Array of hydrogen bond results.

    Returns:
    pandas.DataFrame: A DataFrame containing hydrogen bond results with enriched information.
    """
    df = pd.DataFrame(hbonds,
                      columns=['frame', 'donor_index', 'hydrogen_index',
                               'acceptor_index', 'da_distance', 'dha_angle'])

    df[['da_distance', 'dha_angle']] = df[[
        'da_distance', 'dha_angle']].astype('float16')
    df[['frame', 'donor_index', 'hydrogen_index', 'acceptor_index']] = df[[
        'frame', 'donor_index', 'hydrogen_index', 'acceptor_index']].astype(int)

    df['donor_residue_index'] = universe.atoms.resnums[df['donor_index']]
    df['acceptor_residue_index'] = universe.atoms.resnums[df['acceptor_index']]
    df['donor_atom'] = universe.atoms.names[df['donor_index']]
    df['hydrogen_atom'] = universe.atoms.names[df['hydrogen_index']]
    df['acceptor_atom'] = universe.atoms.names[df['acceptor_index']]
    df['donor_elem'] = [name[0] for name in df['donor_atom']]
    df['acceptor_elem'] = [name[0] for name in df['acceptor_atom']]

    df['hbond_detail'] = [f'{d}-H···{a}' for d,
                          a in zip(df['donor_atom'], df['acceptor_atom'])]
    df['hbond_detail_elem'] = [f'{d}-H···{a}' for d,
                               a in zip(df['donor_elem'], df['acceptor_elem'])]

    return df
    
   
def calculate_hbonds_avg(universe, hbonds_df):
    """
    Computes the average properties of hydrogen bonds over the trajectory.

    Parameters:
    universe (MDAnalysis Universe): The MDAnalysis Universe object containing the molecular structure.
    hbonds_df (pandas.DataFrame): A DataFrame containing hydrogen bond results with enriched information.

    Returns:
    pandas.DataFrame: A DataFrame with columns ['donor_index', 'hydrogen_index', 'acceptor_index', 'n_frames', 'frac', 
                                                 'da_distance', 'dha_angle'], representing the average properties of 
                                                 hydrogen bonds over the trajectory.
    """
    grouping = hbonds_df.groupby(['donor_index', 'hydrogen_index', 'acceptor_index'])

    # Count occurrence and fraction of each Hbond
    count = grouping.size().reset_index(name='n_frames')
    frac = (count['n_frames'] / universe.trajectory.n_frames).astype('float16')
    count['frac'] = frac

    # Calculate mean of DA distance and DHA angle
    mean_dist_ang = grouping[['da_distance', 'dha_angle']].mean().reset_index()

    # Merge the size and mean calculations
    df = pd.merge(count, mean_dist_ang, on=['donor_index', 'hydrogen_index', 'acceptor_index'])
    
    df = df.sort_values(by='frac', ascending=False).reset_index(drop=True)

    return df
   
   
def create_time_series_dict(universe, hbonds_df):
    """
    Creates a dictionary with time series data for each hydrogen bond interaction over the trajectory.

    Parameters:
    universe (MDAnalysis Universe): The MDAnalysis Universe object containing the molecular structure.
    hbonds_df (pandas.DataFrame): A DataFrame containing hydrogen bond results with enriched information.

    Returns:
    dict: A dictionary where keys are tuples representing hydrogen bond interactions (donor_index, acceptor_index),
          and values are numpy arrays representing the time series of the hydrogen bond occurrence over the trajectory.
    """
    time_series_dict = {}

    for _, hb in hbonds_df.iterrows():
        f = hb['frame']
        atm1 = hb['donor_index']
        atm2 = hb['acceptor_index']

        # Possible duplicated interactions (donor-acceptor pairs)
        hb_key1 = (atm1, atm2)
        hb_key2 = (atm2, atm1)

        # Determine the key to use for the time series
        if hb_key1 in time_series_dict:
            hb_key = hb_key1
        elif hb_key2 in time_series_dict:
            hb_key = hb_key2
        else:
            hb_key = hb_key1
            time_series_dict[hb_key] = np.zeros(universe.trajectory.n_frames)

        # Increment the count for the frame where the hydrogen bond is observed
        time_series_dict[hb_key][f] += 1
        
    return time_series_dict



def series_dict_to_df(ser_dict):
    """
    Converts a dictionary of time series data into a pandas DataFrame.

    The dictionary is expected to have keys that are tuples representing pairs of indices, and values that are arrays or lists of time series data corresponding to those pairs.

    Parameters:
    ser_dict (dict): A dictionary where:
        - Keys are tuples, each representing a pair of indices (e.g., (donor_index, acceptor_index)).
        - Values are lists or arrays of time series data for the corresponding pair of indices.

    Returns:
    pandas.DataFrame: A DataFrame with columns:
        - 'donor_index': Index of the donor atom.
        - 'acceptor_index': Index of the acceptor atom.
        - Columns for each time frame's data.

    Notes:
    - The DataFrame will have one row for each key in the dictionary, with columns representing the time series data across different frames.
    - The 'index' column is dropped after extracting 'donor_index' and 'acceptor_index'.

    Example:
    >>> ser_dict = {
    ...     (1, 5): [0, 1, 0, 2],
    ...     (2, 6): [1, 0, 1, 0]
    ... }
    >>> df = series_dict_to_df(ser_dict)
    >>> print(df)
       donor_index  acceptor_index  0  1  2  3
    0            1               5  0  1  0  2
    1            2               6  1  0  1  0
    """
    df = pd.DataFrame.from_dict(ser_dict, orient='index').reset_index()
    df[['donor_index', 'acceptor_index']] = pd.DataFrame(df['index'].tolist())
    df = df.drop('index', axis=1)
    df = df.astype(int)
    
    return df


def find_water_bridges(series_df):
    """
    Identifies water bridges from a DataFrame of hydrogen bond interactions.

    Parameters:
    series_df (pandas.DataFrame): A DataFrame containing hydrogen bond interactions with the following columns:
        - 'donor_index': Index of the donor atom.
        - 'acceptor_index': Index of the acceptor atom.
        - Other columns representing hydrogen bond occurrence over time.

    Returns:
    pandas.DataFrame: A DataFrame containing water bridges with columns:
        - 'atom_index_1': Index of the first atom in the bridge.
        - 'atom_index_2': Index of the second atom in the bridge.
        - 'atom_index_3': Index of the third atom in the bridge.
        - Additional columns indicating the occurrence of the bridge over time.

    Example:
    >>> series_df = pd.DataFrame({
    ...     'donor_index': [1, 2, 3],
    ...     'acceptor_index': [2, 3, 1],
    ...     0: [1, 0, 1],
    ...     1: [0, 1, 0]
    ... })
    >>> bridges = find_water_bridges(series_df)
    >>> print(bridges)
       atom_index_1  atom_index_2  atom_index_3  0  1
    0             1             2             3  1  0
    1             3             1             2  1  0
    """
    bridges = []

    for k, row in series_df.iterrows():
        has_bridges = series_df[series_df['donor_index'] == row['acceptor_index']]

        if len(has_bridges):
            # Remove donor and acceptor atom index columns that are at the end of the df
            row_series = row.iloc[:-2]

            for _, b in has_bridges.iterrows():
                b_series = b.iloc[:-2] & row_series
                b_series = b_series.astype(int)

                # If events occur at the same time as the considered row (indicating at least 1 bridge event)
                if np.sum(b_series) != 0:
                    bridges.append(
                        [row['donor_index'], row['acceptor_index'], b['acceptor_index']] + list(b_series))

    bridges = pd.DataFrame(bridges)
    bridges = bridges.rename(columns={0: 'atom_index_1', 1: 'atom_index_2', 2: 'atom_index_3'})

    return bridges


def add_bridge_details(universe, bridge_df):
    """
    Adds detailed information about atoms and residues to a DataFrame containing bridge interactions.

    Parameters:
    universe (MDAnalysis Universe): The MDAnalysis Universe object containing the molecular structure.
    bridge_df (pandas.DataFrame): A DataFrame with columns:
        - 'atom_index_1': Index of the first atom in the bridge.
        - 'atom_index_2': Index of the second atom in the bridge.
        - 'atom_index_3': Index of the third atom in the bridge.

    Returns:
    pandas.DataFrame: The input DataFrame with additional columns:
        - 'atom_1', 'atom_2', 'atom_3': Names of the atoms corresponding to the indices.
        - 'residue_1', 'residue_2', 'residue_3': Residue numbers of the atoms.
        - 'elem_1', 'elem_2', 'elem_3': Elements of the atoms (first letter of the atom name).
        - 'bridge_detail': String representation of the bridge in the format 'atom1-atom2-atom3'.
        - 'bridge_detail_elem': String representation of the bridge elements in the format 'elem1-elem2-elem3'.

    Example:
    >>> universe = MDAnalysis.Universe('structure.pdb', 'trajectory.dcd')
    >>> bridge_df = pd.DataFrame({
    ...     'atom_index_1': [1, 4],
    ...     'atom_index_2': [2, 5],
    ...     'atom_index_3': [3, 6]
    ... })
    >>> detailed_bridge_df = add_bridge_details(universe, bridge_df)
    >>> print(detailed_bridge_df)
       atom_index_1  atom_index_2  atom_index_3 atom_1 atom_2 atom_3  residue_1  residue_2  residue_3 elem_1 elem_2 elem_3   bridge_detail bridge_detail_elem
    0             1             2             3      N      CA      C        100        100        100      N      C      C      N-CA-C             N-C-C
    1             4             5             6      N      CA      C        101        101        101      N      C      C      N-CA-C             N-C-C
    """
    df = bridge_df.copy()
    
    df['atom_1'] = universe.atoms.names[df['atom_index_1']]
    df['atom_2'] = universe.atoms.names[df['atom_index_2']]
    df['atom_3'] = universe.atoms.names[df['atom_index_3']]

    df['residue_1'] = universe.atoms.resnums[df['atom_index_1']]
    df['residue_2'] = universe.atoms.resnums[df['atom_index_2']]
    df['residue_3'] = universe.atoms.resnums[df['atom_index_3']]

    df['elem_1'] = [name[0] for name in df['atom_1']]
    df['elem_2'] = [name[0] for name in df['atom_2']]
    df['elem_3'] = [name[0] for name in df['atom_3']]

    df['bridge_detail'] = [f'{a1}-{a2}-{a3}' for a1, a2, a3 in zip(df['atom_1'], df['atom_2'], df['atom_3'])]
    df['bridge_detail_elem'] = [f'{a1}-{a2}-{a3}' for a1, a2, a3 in zip(df['elem_1'], df['elem_2'], df['elem_3'])]
    
    return df

def count_bridges_type(bridges):
    """
    Counts the occurrences of each type of water bridge based on the bridge detail and returns a DataFrame.

    Parameters:
    bridges (pandas.DataFrame): A DataFrame containing water bridges with a column 'bridge_detail' that describes the bridge.

    Returns:
    pandas.DataFrame: A DataFrame with two columns - 'bridge' and 'count', where 'bridge' is the bridge detail and 'count' is the number of occurrences of that bridge type.

    Example:
    >>> bridges = pd.DataFrame({
    ...     'bridge_detail': ['A-B-C', 'A-B-C', 'D-E-F', 'A-B-C']
    ... })
    >>> counts_df = count_bridges_type(bridges)
    >>> print(counts_df)
       bridge  count
    0  A-B-C      3
    1  D-E-F      1
    """
    count = bridges.groupby('bridge_detail').size()
    
    df = pd.DataFrame({'bridge': count.index, 'count': count.values}) 
    
    return df


