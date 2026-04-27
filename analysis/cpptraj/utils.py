import numpy as np
import pandas as pd

def read_hbonds_details(file):
    """
    Read hydrogen bond details from the avgout action of Cpptraj command hbonds and process them.

    Args:
        file (str): Path to the file containing hydrogen bond details.

    Returns:
        pandas.DataFrame: Processed data containing hydrogen bond details.
            - hbond: Cpptraj naming convention of hbonds in the hbond series results
            - donor, donorh, acceptor: Cpptraj string representing the hbond 
                                       (:u[donor/donorh/acceptor]_u[1/2]@atm[1/2])
            - n_frames: number of frames the hbond exists
            - frac: fraction of frames the hbond exists
            - dist: average distance
            - ang: average angle
            - u1 and u2: units index in hbond
            - atm1 and atm2: respective atoms of u1 and u2 in hbond
            - u_donor and u_acceptor: units in hbond
    """
    # Read and clean data
    dat = pd.read_table(file,
                        delim_whitespace=True,
                        comment='#',
                        header=None,
                        names=['acceptor', 'donorh', 'donor', 'n_frames',
                               'frac', 'dist', 'ang']
                        )

    # Split 'donor' column into unit, donor index, and donor atom parts
    dat[['u_donor', 'u1', 'atm1']] = \
        dat['donor'].str.split(pat=':|@|_', expand=True)
    
    # Split 'acceptor' column into unit, acceptor index, and acceptor atom parts
    dat[['u_acceptor', 'u2', 'atm2']] = \
        dat['acceptor'].str.split(pat=':|@|_', expand=True)

    # Create hbonds naming convention for coherence with the naming of the hbond series
    donorh = [d.split('@')[-1] for d in dat['donorh']] # Extract donor hydrogen atom name
    dat['hbond'] = ['-'.join([dat['acceptor'].iloc[hb], dat['donor'].iloc[hb], donorh[hb]])
                    for hb in range(dat.shape[0])]
    
    # Assign data types to each column
    dat[['donor', 'donorh', 'acceptor', 'u_donor', 'atm1', 'u_acceptor', 'atm2']] = dat[['donor', 'donorh', 'acceptor', 'u_donor', 'atm1', 'u_acceptor', 'atm2']].astype(str)
    dat[['u1', 'u2', 'n_frames']] = dat[['u1', 'u2', 'n_frames']].astype(int)
    dat[['dist', 'ang', 'frac']] = dat[['dist', 'ang', 'frac']].astype(float)
    
    # Separate acceptor and donor elements
    dat['atm1_elem'] = [i[0] for i in list(dat['atm1'])]
    dat['atm2_elem'] = [i[0] for i in list(dat['atm2'])]
    
    # Define H-bond type
    dat['hbond_type'] = [i['atm1_elem'] + '-H···' + i['atm2_elem'] for k, i in dat.iterrows()]
    dat['hbond_type_detail'] = [i['atm1'] + '-H···' + i['atm2'] for k, i in dat.iterrows()]

    return dat

def read_hbonds_series(file_path, details):
    """
    Read hydrogen bond series from a file and process them.

    Args:
        file_path (str): Path to the file containing hydrogen bond series.
        details (pandas.DataFrame): DataFrame containing hbond details.
                                    From read_hbonds_details().

    Returns:
        pandas.DataFrame: Processed data containing hydrogen bond series.
    """
    # Read data from the file and transpose
    mat = pd.read_table(file_path, delim_whitespace=True).T.reset_index()

    # Remove the "#Frame" index
    mat = mat.loc[1:, :]

    # Select the filtered hbonds from read_hbonds_details()
    mat = mat[mat['index'].isin(details['hbond'])]


    # # Split the index column into separate columns for units and atoms
    mat[['u_acceptor', 'u1', 'atm1', 'u_donor', 'u2', 'atm2', 'u_donorh']] = \
        mat['index'].str.split(pat='-|@|_', expand=True)

    # # Reorder data such that every column not corresponding to the series is at the beginning
    tmp1 = mat[['index', 'u_acceptor', 'u1', 'atm1', 'u_donor', 'u2', 'atm2', 'u_donorh']]
    tmp2 = mat.drop(['index', 'u_acceptor', 'u1', 'atm1', 'u_donor', 'u2', 'atm2', 'u_donorh'], axis=1)

    mat = pd.concat([tmp1, tmp2], axis=1)

    # Convert unit columns to integer type
    mat[['u1','u2']] = mat[['u1','u2']].astype(int)

    mat = mat.rename({'index': 'hbond'}, axis=1)

    return mat

def contact_map(details, n_units): # Common function between contacts and hbonds
    """
    Generate a contact map from a DataFrame of contacts.

    Args:
        details (pandas.DataFrame): DataFrame containing contact details.
                                    From read_[contacts/hbonds]_details().
        n_units (int): Number of units.

    Returns:
        numpy.ndarray: Contact map representing the number of frames a contact exists
                       between units.
    """
    mat = np.zeros((n_units, n_units))

    # Loop through each row of the contacts DataFrame
    for _, v in details.iterrows():
        # Increment the contact map for the corresponding units
        mat[v['u1'] - 1, v['u2'] - 1] += v['n_frames']
        mat[v['u2'] - 1, v['u1'] - 1] += v['n_frames']

    return mat