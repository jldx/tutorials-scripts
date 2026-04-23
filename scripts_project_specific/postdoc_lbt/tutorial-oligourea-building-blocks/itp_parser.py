#! python3

import pandas as pd

def itp_parser(itp_file_path):
    """
    Parse a GROMACS topology (.itp) file and extract its contents.

    Args:
        itp_file_path (str): Path to the GROMACS topology file (.itp format).

    Returns:
        dict: A dictionary containing the parsed topology data where:
            - Keys (str): Section headers from the .itp file. Note that 'dihedrals' 
              sections are automatically separated into 'dihedrals' and 'impropers' 
              based on their function type.
            - Values (list of list): Each value is a list containing the data rows 
              for that section, with each row represented as a list of strings.

    Example:
        >>> topology_data = parse_itp_file("protein.itp")
        >>> print(topology_data.keys())
        dict_keys(['atoms', 'bonds', 'angles', 'dihedrals', 'impropers'])
        >>> print(topology_data['atoms'][0])  # First atom entry
        ['1', 'CA', '1', 'ALA', 'CA', '1', '0.0000', '12.01']
    """
    # Initialize storage for sections and tracking variables
    sections = {}  # Dictionary to store section name -> content mapping
    previous_header = None  # Track the previous section title for duplicate detection
    current_header = None   # Current section being processed
    current_content = []   # Content lines for the current section
    
    with open(itp_file_path, 'r') as fin:
        for line in fin:
            # Check if line is a section header
            is_header = line.strip().startswith('[') and line.strip().endswith(']')
            # Check if line contains actual data
            is_content = line[0] not in ['#', ';', '\n']
            
            if is_header:
                # Save the previous section's content before starting a new one
                if current_header is not None:
                    sections[current_header] = current_content
                    current_content = []  # Reset content list for new section
                
                # Extract header name and handle duplicate header name:
                # case of dihedrals and impropers sharing the same header name
                previous_header = current_header
                current_header = line.strip()[1:-1].strip()  # Remove brackets and whitespace

                if previous_header == current_header:
                    current_header = 'impropers'
            
            # Add content lines to current section
            elif current_header is not None and is_content:
                current_content.append(line.split())
    
    # Save the last section's content
    if current_header is not None:
        sections[current_header] = current_content
    
    return sections
    
    
def df_to_dict(df, key_name, value_name):
    return {a[key_name]: a[value_name] for k, a in df.iterrows()}
    
    


