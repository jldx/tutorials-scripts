#! python3

# Processes all ITP files in the input folder and generates the overall forcefield
# Usage: python3 merge_forcefields.itp <residues forcefield itp files folder> <output file itp>

import argparse
import os
import sys

import pandas as pd
from utils import *

def update_content(forcefield, header):
    df = (pd
          .concat(forcefield[header], axis=0)
          .fillna(value='')
          .drop_duplicates())
    
    if header == 'atomtypes':
        df = (df.iloc[:, :7]
              .rename(columns={0: 'name', 1:'number', 2: 'mass', 3: 'charge', 
                               4: 'type', 5: 'sigma', 6: 'epsilon'})
              .assign(number=lambda df: df['number'].astype(float))
              .assign(mass=lambda df: df['mass'].astype(float))
              .assign(charge=lambda df: df['charge'].astype(float))
              .assign(sigma=lambda df: df['sigma'].astype(float))
              .assign(epsilon=lambda df: df['epsilon'].astype(float))
              .groupby(by=['name', 'number', 'mass', 'type', 'sigma', 'epsilon']).mean()
              .reset_index())[['name', 'number', 'mass', 'charge', 'type', 'sigma', 'epsilon']]

    return df



def updated_content_line(header, content):
    if header == 'atomtypes':
        line = (f'{content[0]:>7} {int(content[1]):>6} {float(content[2]):>10.4f} '
                f'{float(content[3]):>10.3f} {content[4]:>5} {float(content[5]):>20.11e} '
                f'{float(content[6]):>15.6e} \n')
        
    elif header == 'nonbond_params':
        line = (f'{content[0]:>7} {content[1]:>7} {int(content[2]):>5}  {float(content[3]):>4.11e}  '
                f'{float(content[4]):>4.11e} \n')
            
    elif header in ['bondtypes', 'pairtypes']:
        line = (f'{content[0]:>7} {content[1]:>7} {int(content[2]):>5}  {float(content[3]):>4.6e}  '
                f'{float(content[4]):>4.6e} \n')
            
    elif header == 'angletypes':
        line = (f'{content[0]:>7} {content[1]:>7} {content[2]:>7} {int(content[3]):>5}  '
                f'{float(content[4]):>4.7e}  {float(content[5]):>4.7e}  {float(content[6]):>4.7e} '
                f'{float(content[7]):>4.7e} \n')
        
    else:
        try:
            line = (f'{content[0]:>7} {content[1]:>7} {content[2]:>7} {content[3]:>7} '
                    f'{int(content[4]):>5} {float(content[5]): 4.6e} {float(content[6]): 4.6e} '
                    f'{content[7]:>6} \n') 
        except:
            # Case where the multiplicity is absent for impropers
            line = (f'{content[0]:>7} {content[1]:>7} {content[2]:>7} {content[3]:>7} '
                    f'{int(content[4]):>5} {float(content[5]): 4.6e} {float(content[6]): 4.6e} \n') 
        
    return line
  

def parse_arguments():
    """
    Parse command line arguments for the forcefield name and, input and output files.
    
    Returns
    -------
    argparse.Namespace
        Parsed arguments containing validated ones
    
    Raises
    ------
    SystemExit
        If any specified folder does not exist.
    """
    parser = argparse.ArgumentParser(description="Process folders and residue prefix.")
    parser.add_argument('-i', '--input', required=True, help='Path to the input folder')
    parser.add_argument('-o', '--output', required=True, help='Path to the forcefield file')
    args = parser.parse_args()
    
    return args 

def main():
    """
    Main function to process molecular topology files.
    
    Processes all ITP files in the input folder and generates the overall forcefield
    The function is currently hardcoded for specific residue naming conventions.
    """
    # Parse command line arguments
    args = parse_arguments()
    input_folder = args.input
    output_file = args.output

    # List fo the molecule names to parse of format 'GxG' with an exception of GHE
    molnames = [f'G{aa}G' for aa in ['A', 'C', 'D', 'E', 'F', 'G', 'I', 'K', 'L', 
                                     'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']] + ['GHE']
    
    # Add caps
    molnames += ['BAN', 'BGN', 'BPN']

    # Dictionary containing ITP file section headers with their format descriptions
    # Each key represents a section name, and the value contains the section header
    itp_headers = {
        # Main header of a forcefield ITP file
        'defaults': '[ defaults ]\n; nbfunc	comb-rule	gen-pairs	fudgeLJ	fudgeQQ\n1	2	yes	1.000000	1.000000\n',
        
        # Atoms section: defines individual atoms in the forcefield
        'atomtypes': '\n[ atomtypes ]\n; name	at.num	mass	charge	ptype	sigma	epsilon	\n',
        
        # Non-bonded parameters
        'nonbond_params': '\n[ nonbond_params ]\n; i	j	func	sigma	epsilon\n',
        
        # Bonds section: defines covalent bonds between atoms
        'bondtypes': '\n[ bondtypes ]\n; i	j	func	b0	Kb\n',
        
        # Pairs section: defines 1-4 interactions (typically for non-bonded interactions)
        'pairtypes': '\n[ pairtypes ]\n; i	j	func	sigma1-4	epsilon1-4\n',
        
        # Dihedrals section: defines proper dihedral angles between four atoms
        'angletypes': '\n[ angletypes ]\n; i	j	k	func	th0	Kth	s0	Kub\n',
        
        # Dihedrals section: defines proper dihedral angles between four atoms
        'dihedraltypes': '\n[ dihedraltypes ]\n; i	j	k	l	func	phi0	Kphi	mult\n',
        
        # Impropers section: defines improper dihedral angles (out-of-plane bending)
        # Note: Uses same [ dihedrals ] header as proper dihedrals in ITP files
        'impropers': '\n[ dihedraltypes ]\n; i	j	k	l	func	q0	Kq\n'}
        
    headers_names = list(itp_headers.keys())[1:]

    forcefield = {header: [] for header in headers_names}

    for mol in molnames:
        f = f'{input_folder}/{mol}_forcefield.itp'
        ff = itp_parser(f)
        del ff['defaults']
        
        for header in headers_names:
            forcefield[header].append(pd.DataFrame(ff[header]))
            
    with open(output_file, 'w') as fout:
        
        for header in headers_names:
            updated_section = update_content(forcefield, header)

            if header == 'atomtypes':
                fout.write(itp_headers['defaults'])
            
            fout.write(itp_headers[header])
                
            for k, content in updated_section.iterrows():
                line = updated_content_line(header, content)
                fout.write(line)
                
if __name__ == '__main__':
    main()
