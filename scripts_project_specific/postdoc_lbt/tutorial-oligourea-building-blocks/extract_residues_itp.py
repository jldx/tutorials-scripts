#!python3

# Description: create a Gromacs ITP file of a residue of interest included in a 
# 'peptide' from the ITP file of the 'peptide'.

# Usage example: python3 residues_itp.py -i input_folder -o output_folder -m residues_mapping_folder

# - The ITP in input_folder should have the same name as the mapping file in residue_mapping_folder 
# - The output file is the name of the residue included in GxG
# - E.g: For GAG: from GAG.itp and GAG.map -> UUA.itp
#        The UU prefix is hardcoded (see function main())

import argparse
import os
import sys

import pandas as pd
import glob
import re
from utils import *


def map_parser(molmap_file_path):
    df = (pd.read_table(molmap_file_path, sep='\\s+', comment='#', header=None)
          .rename(columns={0: 'index_ori', 1: 'name_ori', 2: 'name', 3: 'name_ori_index'})
          .astype({'index_ori': str}))
    return df


def molmap_check_in_mol(content_list):
    for i in content_list:
        
        if i == '/':
            return False

        try:
            int(i)
            return True

        except ValueError:
            pass

    return False

def update_content(header, section_content, map_ids, map_names, n_ids, resname):
    
    if header == 'atoms':

        section_content = (pd.DataFrame(section_content)
                   .iloc[:, [0, 1, 2, 3, 4, 6, 7]]
                   .rename(columns={0: 'index_ori', 1: 'atomtype', 2: 'resid', 
                                    3: 'resname', 4: 'name_ori', 6: 'charge', 7: 'mass'})
                   .assign(index=lambda df: df['index_ori'].map(map_ids))
                   .assign(name_res=lambda df: df['index_ori'].map(map_names))
                   .assign(index_ori=lambda df: df['index_ori'].astype(int)))
        
        section_content_new = (section_content[section_content['index'].str.isnumeric()]
                               .assign(index=lambda df: df['index'].astype(int))
                               .assign(charge_group=lambda df: df['index'])
                               .assign(resname=resname)
                               .sort_values('index')
                               .reset_index(drop=True)
                               [['index', 'atomtype', 'resid', 'resname', 'name_res', 
                                 'charge_group', 'charge', 'mass']]
                               .values.tolist())
        
    else:

        section_content_new = []

        for content in section_content:

            content_ids = content[:n_ids]
            content_param = content[n_ids:]

            content_ids_new = [map_ids[i] for i in content_ids]
            content_ids_new_param = content_ids_new + content_param
            
            in_mol = molmap_check_in_mol(content_ids_new)
            is_duplicate = content_ids_new_param in section_content_new
                
            if in_mol & (not is_duplicate):
                section_content_new.append(content_ids_new_param)

    return section_content_new


def update_sections(resname, top, molmap):
    
    header_n_ids = {'atoms': 0, 'bonds': 2, 'pairs': 2,
                      'angles': 3, 'dihedrals': 4, 'impropers': 4}
    header_names = list(header_n_ids.keys())
    
    map_names = df_to_dict(molmap, 'index_ori', 'name')
    map_ids = df_to_dict(molmap, 'index_ori', 'name_ori_index')
    
    sections_updated = {}

    for header in header_names:

        section_content = top[header]
        n_ids = header_n_ids[header]

        section_updated = update_content(
            header, section_content, map_ids, map_names, n_ids, resname)

        sections_updated[header] = section_updated
        
    return sections_updated


def write_mol_itp(itp_file_path, itp_headers, sections):
    with open(itp_file_path, 'w') as fout:

        header_names = list(itp_headers.keys())[:-1]

        fout.write(itp_headers['moleculetype'])

        for header in header_names:

            fout.write(itp_headers[header])

            for c in sections[header]:

                if header == 'atoms':
                    fout.write(f'{c[0]:>6} {c[1]:>10} {c[2]:>6} {c[3]:>8} {c[4]:>6}')
                    fout.write(f'{c[5]:>7} {float(c[6]):>12.6f} {float(c[7]):>10.4f}\n')

                elif header in ['bonds', 'pairs']:
                    fout.write(f'{c[0]:>5} {c[1]:>5} {c[2]:>5}\n')

                elif header == 'angles':
                    fout.write(f'{c[0]:>5} {c[1]:>5} {c[2]:>5} {c[3]:>5}\n')

                else:
                    fout.write(f'{c[0]:>5} {c[1]:>5} {c[2]:>5} {c[3]:>5} {c[4]:>5}\n')     


def parse_arguments():
    """
    Parse command line arguments for input, output, and mapping folders.
    
    Returns
    -------
    argparse.Namespace
        Parsed arguments containing validated folder paths.
    
    Raises
    ------
    SystemExit
        If any specified folder does not exist.
    """
    parser = argparse.ArgumentParser(description="Process folders and residue prefix.")
    parser.add_argument('-i', '--input', required=True, help='Path to the input folder')
    parser.add_argument('-o', '--output', required=True, help='Path to the output folder')
    parser.add_argument('-m', '--mapping', required=True, help='Path to the mapping folder')
    args = parser.parse_args()
    
    return args

def main():
    """
    Main function to process molecular topology files.
    
    Processes all ITP files in the input folder, applies atom mappings,
    and generates updated topology files with new residue names.
    The function is currently hardcoded for specific residue naming conventions.
    """
    # Parse command line arguments
    args = parse_arguments()
    input_folder = args.input
    output_folder = args.output
    mapping_folder = args.mapping

    # Dictionary containing ITP file section headers with their format descriptions
    # Each key represents a section name, and the value contains the section header
    itp_headers = {
        # Atoms section: defines individual atoms in the molecule
        'atoms': '\n[ atoms ]\n; nr type resnr residu atom cgnr charge mass\n',
        
        # Bonds section: defines covalent bonds between atoms
        'bonds': '\n[ bonds ]\n; ai aj funct b0 Kb\n',
        
        # Pairs section: defines 1-4 interactions (typically for non-bonded interactions)
        'pairs': '\n[ pairs ]\n; ai aj funct c6 c12 or\n; ai aj funct fudgeQQ q1 q2 c6 c12\n',
        
        # Angles section: defines bond angles between three atoms
        'angles': '\n[ angles ]\n; ai aj ak funct th0 cth S0 Kub\n',
        
        # Dihedrals section: defines proper dihedral angles between four atoms
        'dihedrals': '\n[ dihedrals ]\n; ai aj ak al funct phi0 cp mult\n',
        
        # Impropers section: defines improper dihedral angles (out-of-plane bending)
        # Note: Uses same [ dihedrals ] header as proper dihedrals in ITP files
        'impropers': '\n[ dihedrals ]\n; ai aj ak al funct q0 cq\n'
    }
    
    # Process all files in input folder
    files_list = glob.glob(f'{input_folder}/*')
    for f in files_list:
        # Extract molecule name from file path
        f_splitted = re.split('/|\.', f)
        molname = f_splitted[-2]
        
        # Generate residue name (hardcoded for specific naming convention)
        # resname can be changed to accommodate other residues
        
        resname = f'UU{molname[1]}' # For a residue IN a molecule
        #resname = molname # For the caps Boc and NHMe
        
        # Update moleculetype section with new residue name
        itp_headers['moleculetype'] = f'[ moleculetype ]\n; name nrexcl\n{resname} 3\n'
        
        # Process topology file
        top = itp_parser(f)
        molmap = map_parser(f'{mapping_folder}/{molname}.map')
        sections_updated = update_sections(resname, top, molmap)
        
        # Write updated topology file
        write_mol_itp(
            f'{output_folder}/{resname}.itp', itp_headers, sections_updated)


if __name__ == '__main__':
    main()
