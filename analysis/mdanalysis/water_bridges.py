#!python3

import numpy as np
import pandas as pd
import MDAnalysis as mda

import warnings
warnings.filterwarnings('ignore')

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


