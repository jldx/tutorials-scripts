import pandas as pd
import numpy as np

def sort_df(df, column_name):
    """
    Sorts a DataFrame by a specified column and resets the index.

    Parameters:
    df (pandas.DataFrame): The DataFrame to be sorted.
    column_name (str): The column name by which to sort the DataFrame.

    Returns:
    pandas.DataFrame: The sorted DataFrame with a reset index.

    Example:
    >>> data = {'A': [3, 1, 2], 'B': [9, 8, 7]}
    >>> df = pd.DataFrame(data)
    >>> sorted_df = sort_df(df, 'A')
    >>> sorted_df
       A  B
    0  1  8
    1  2  7
    2  3  9
    """
    df = df.sort_values(by=column_name).reset_index(drop=True)
    return df
 
def outliers_idx(series, std_thresh):
   """
   Identifies the indices of outliers in a series based on a specified standard deviation threshold.

   Parameters:
   series (array-like): An array-like object (such as a list, numpy array, or pandas Series) of numerical data.
   std_thresh (float): The standard deviation threshold for identifying outliers.

   Returns:
   numpy.ndarray: An array of indices corresponding to the outliers in the input series.

   Example:
   >>> series = [10, 12, 14, 15, 17, 18, 19, 100]
   >>> outliers_idx(series, 2)
   array([7])
   """
   mean = np.mean(series)
   std_dev = np.std(series)
   outliers = np.where((series < mean - std_thresh * std_dev) | 
                     (series > mean + std_thresh * std_dev))[0]
   return outliers
   
def extract_coord(u, select, coord_axis):
    """
    Extracts the coordinates of selected atoms along a specified axis from a trajectory.

    Parameters:
    u (MDAnalysis.Universe): The MDAnalysis universe containing the trajectory.
    select (MDAnalysis.AtomGroup): The atom group selection for which coordinates are to be extracted.
    coord_axis (int): The axis along which to extract coordinates (0 for x, 1 for y, 2 for z).

    Returns:
    np.ndarray: An array of dimensions (n_frames, 1, n_coord) containing the extracted coordinates.

    Example:
    >>> u = mda.Universe(topology, trajectory)
    >>> select = u.select_atoms("name OW")
    >>> coord_axis = 0
    >>> coords = extract_coord(u, select, coord_axis)
    """
    pos = []
    for ts in u.trajectory[1:]:
        pos.append([select.positions[:, coord_axis]])

    return np.array(pos)

