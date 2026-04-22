#! python3

import sys
import os
import psutil
import argparse
import time
from multiprocessing import Pool, cpu_count
import numpy as np
import pandas as pd
import MDAnalysis as mda
import joblib

def welcome() -> None:
    welcome = r"""
--------------------------------------------------------------------------------

        =()=
    ,/'\_||_                       PERMEATIONS COUNTER
    ( (___  `.            Python programm computing permeation events 
    `\./  `=='        of a molecule through a membrane from MD simulations  
          |||                
          |||             Version 2, last updated: 15 oct. 2025
          |||            

         Version 1 by Arthur Hardiagon        Version 2 by Julie Ledoux            
                    Laboratory of Theoretical Biochemistry
                        CNRS UMR 8266, Paris, France

Original code available at https://github.com/ahardiag/perm-md-count

--------------------------------------------------------------------------------
"""
    print(welcome)

def print_help() -> None:
    help_text = """
Usage: python3 permeations_counter.py --topology <topology_file> --trajectory <trajectory_file> --output <output_prefix>

Arguments:
  --topology     Path to topology file (MDAnalysis supported format)
  --trajectory   Path to trajectory file (MDAnalysis supported format)
  --output       Prefix for the output files

Example:
  python permeations_counter.py --topology top.pdb --trajectory traj.xtc --output result

Note: Make sure to execute the script in the provided conda environment

Additional ressources:
    For MDAnalysis supported files format, see: https://userguide.mdanalysis.org/stable/formats/index.html
"""
    print(help_text)

def parse_and_validate_args() -> 'argparse.Namespace':
    """
    Parses command-line arguments for topology, trajectory, and output file paths.

    Validates that the topology and trajectory files exist.
    If the --help flag is provided, displays usage information and exits.

    Returns:
        argparse.Namespace: Parsed arguments with attributes 'topology', 'trajectory', and 'output'.

    Exits:
        If required files do not exist or if help is requested.
    """
    parser = argparse.ArgumentParser(description="Python programm computing permeation events through a membrane", add_help=False)
    parser.add_argument('--topology', type=str, required=False, help='Path to topology file')
    parser.add_argument('--trajectory', type=str, required=False, help='Path to trajectory file')
    parser.add_argument('--output', type=str, required=False, help='Prefix for the output files')
    parser.add_argument('--help', action='store_true', help='Show this help message and exit')

    args, unknown = parser.parse_known_args()

    if args.help:
        print_help()
        sys.exit(0)

    # Validate that topology and trajectory files exist
    if not os.path.isfile(args.topology):
        print(f"Error: Topology file '{args.topology}' does not exist.", file=sys.stderr)
        sys.exit(1)
    if not os.path.isfile(args.trajectory):
        print(f"Error: Trajectory file '{args.trajectory}' does not exist.", file=sys.stderr)
        sys.exit(1)

    return args

def convert_seconds_to_hours_min(n_seconds: int) -> str:
    """
    Converts a number of seconds to a string formatted as hours, minutes, and seconds.

    Args:
        n_seconds (int): The total number of seconds to convert.

    Returns:
        str: A string in the format "{hours}h {minutes}min {seconds}s".
    """
    hours, remainder = divmod(n_seconds, 3600)
    minutes, seconds = divmod(remainder, 60)
    return f"{int(hours)}h {int(minutes)}m {int(seconds)}s"

def have_enough_memory(n_processes: int, file_size: int) -> bool:
    """
    Check if the required memory for running n_processes, each handling a file of size file_size bytes,
    can fit in the available system RAM.

    Args:
        n_processes (int): Number of processes to run.
        file_size (int): Size of each file in bytes.

    Returns:
        bool: True if enough memory is available, False otherwise.

    Prints:
        Information about required and available memory, and suggestions if memory is insufficient.
    """
    print("Checking available memory...")
    available: int = psutil.virtual_memory().available
    required: int = n_processes * file_size

    if required <= available:
        #print(f"Enough memory: need {required/(1024**2):.0f} MB, available {available/(1024**2):.0f} MB")
        print(f"Proceeding with {n_processes} CPUs.\n")
        return True
    else:
        print(f"Not enough memory: need {required/(1024**2):.0f} MB, available {available/(1024**2):.0f} MB.")
        print("Consider splitting the simulation into smaller chunks.")
        return False


def process_z_coords_extraction(
    args: tuple[
        int, int, str, str, str, str, float
    ]
) -> tuple[int, int, np.ndarray, np.ndarray, np.ndarray]:
    """
    Extracts the z-coordinates of membrane boundaries and selected molecules for a chunk of trajectory frames.

    Args:
        args (tuple): Contains:
            start_idx (int): Start frame index for this chunk.
            end_idx (int): End frame index for this chunk.
            traj_path (str): Path to trajectory file.
            top_path (str): Path to topology file.
            membrane_selection_str (str): Atom selection string for membrane.
            mols_selection_str (str): Atom selection string for molecules.
            MAX_DISPLACEMENT (float): Threshold for PBC jump detection.

    Returns:
        tuple: (start_idx, end_idx, lower_boundaries_chunk, upper_boundaries_chunk, mols_z_chunk)
            lower_boundaries_chunk (np.ndarray): Mean z of lower membrane boundary per frame.
            upper_boundaries_chunk (np.ndarray): Mean z of upper membrane boundary per frame.
            mols_z_chunk (np.ndarray): z-coordinates of selected molecules per frame.
    """
    # Get arguments
    start_idx, end_idx, traj_path, top_path, membrane_selection_str, mols_selection_str, MAX_DISPLACEMENT = args

    # Each process needs its own MDAnalysis Universe to avoid thread conflicts
    import MDAnalysis as mda
    u_local = mda.Universe(top_path, traj_path, in_memory=False)
    n_frames = end_idx - start_idx

    # Get box dimensions and membrane center
    box_x, box_y, box_z = u_local.coord.dimensions[:3]
    box_center = box_z / 2
    upper_memb_z_min, upper_memb_z_max = box_center, box_center + MAX_DISPLACEMENT
    lower_memb_z_max, lower_memb_z_min = box_center, box_center - MAX_DISPLACEMENT 

    # Select membrane atoms for upper and lower boundaries
    memb_sel = u_local.select_atoms(membrane_selection_str)
    upper_memb_sel = memb_sel.select_atoms(f"prop z > {upper_memb_z_min} and prop z < {upper_memb_z_max}", updating=True)
    lower_memb_sel = memb_sel.select_atoms(f"prop z < {lower_memb_z_max} and prop z > {lower_memb_z_min}", updating=True)
    # Select molecule atoms
    mols_sel = u_local.select_atoms(mols_selection_str)
    n_mols_atoms = mols_sel.n_atoms

    # Initialize arrays to store results for this chunk
    lower_boundaries_chunk = np.empty(n_frames, dtype=np.float64)
    upper_boundaries_chunk = np.empty(n_frames, dtype=np.float64)
    mols_z_chunk = np.empty((n_mols_atoms, n_frames), dtype=np.float64)

    # Loop over frames in this chunk and extract z positions
    for i, ts in enumerate(u_local.trajectory[start_idx:end_idx]):
        # Mean z position for lower and upper membrane boundaries
        lower_boundaries_chunk[i] = lower_memb_sel.positions[:, 2].mean()
        upper_boundaries_chunk[i] = upper_memb_sel.positions[:, 2].mean()
        # z positions for all selected molecules
        mols_z_chunk[:, i] = mols_sel.positions[:, 2]
    
    return start_idx, end_idx, lower_boundaries_chunk, upper_boundaries_chunk, mols_z_chunk


def process_permeations_search(
    args: tuple[
        int, int, np.ndarray, np.ndarray, np.ndarray, np.ndarray, float
    ]
) -> list:
    """
    Identifies permeation events for a chunk of molecules over a range of frames.

    Args:
        args (tuple): Contains:
            start_idx (int): Start molecule index for this chunk.
            end_idx (int): End molecule index for this chunk.
            upper_boundaries (np.ndarray): Array of upper membrane boundary z-coordinates per frame.
            lower_boundaries (np.ndarray): Array of lower membrane boundary z-coordinates per frame.
            mols_z_prev (np.ndarray): z-coordinates of molecules at previous frame.
            mols_z_curr (np.ndarray): z-coordinates of molecules at current frame.
            MAX_DISPLACEMENT (float): Threshold for PBC jump detection.

    Returns:
        list: List of permeation events for this chunk. Each event is a list:
            [molecule_index, start_frame, end_frame, duration_frames, direction]
    """
    # Get arguments
    start_idx, end_idx, upper_boundaries, lower_boundaries, mols_z_prev, mols_z_curr, MAX_DISPLACEMENT = args

    # Select molecule indices for this chunk
    mols_idx_chunk = np.arange(start_idx, end_idx)

    # Select z-coordinates for this chunk
    mols_z_prev_chunk = mols_z_prev[start_idx:end_idx, :]
    mols_z_curr_chunk = mols_z_curr[start_idx:end_idx, :]

    # Determine regions for each molecule at each frame
    in_upper_bulk = mols_z_curr_chunk > upper_boundaries
    in_lower_bulk = mols_z_curr_chunk < lower_boundaries
    in_membrane = np.invert(in_upper_bulk | in_lower_bulk)

    # Determine membrane leaflets
    half_memb = (upper_boundaries + lower_boundaries) / 2
    in_upper_memb = mols_z_curr_chunk > half_memb
    in_upper_leaflet = in_upper_memb & np.invert(in_upper_bulk)
    in_lower_memb = mols_z_curr_chunk < half_memb
    in_lower_leaflet = in_lower_memb & np.invert(in_lower_bulk)

    # Detect PBC jumps
    pbc_jump = abs(mols_z_curr_chunk - mols_z_prev_chunk) > MAX_DISPLACEMENT
    pbc_jump_prev = pbc_jump[:, :-1]
    pbc_jump_curr = pbc_jump[:, 1:]

    # Encode trajectory regions
    in_lower_bulk_prev = in_lower_bulk[:, :-1] * -1
    in_upper_bulk_prev = in_upper_bulk[:, :-1] * 1
    in_membrane_prev = in_membrane[:, :-1] * 0
    in_upper_leaflet_prev = in_upper_leaflet[:, :-1] * 1
    in_lower_leaflet_prev = in_lower_leaflet[:, :-1] * -1

    # Apply correction for PBC jumps in membrane
    pbc_correction = pbc_jump_prev * (in_upper_leaflet_prev + in_lower_leaflet_prev) + \
        pbc_jump_curr * (in_upper_leaflet_prev + in_lower_leaflet_prev)

    regions_chunk = in_lower_bulk_prev + \
        in_upper_bulk_prev + in_membrane_prev + pbc_correction

    # Search for permeation events
    permeations_chunk = []
    n_mols_chunk, n_frames = regions_chunk.shape

    # Track region history for each molecule
    history = np.zeros_like(regions_chunk, dtype=int)
    history[:, 0] = regions_chunk[:, 0]

    for i in range(1, n_frames):
        in_upper = (regions_chunk[:, i] == 1)
        in_lower = (regions_chunk[:, i] == -1)
        in_membrane = (regions_chunk[:, i] == 0)

        history[in_upper, i] = 1
        history[in_lower, i] = -1
        # For membrane: increment/decrement based on previous label
        prev_label = history[:, i-1]
        history[in_membrane, i] = prev_label[in_membrane] + np.sign(prev_label[in_membrane])

    # Identify jump types
    #     PBC jump: either -2 or 2
    #     Exit in the upper bulk: 1
    #     Exit in the lower bulk: -1
    #     Stays in the same region: 0    
    regions_prev_chunk = regions_chunk[:, :-1]
    regions_curr_chunk = regions_chunk[:, 1:]
    jump_type = regions_curr_chunk - regions_prev_chunk

    frame_indices = np.arange(n_frames - 1)
    departed_from_membrane = (regions_curr_chunk - jump_type == 0)
    opposite_direction = (history[:, :-1] * jump_type < 0)
    events_sparse = departed_from_membrane * opposite_direction * jump_type * frame_indices

    # Extract events for each molecule
    for idx in range(n_mols_chunk):
        events_bool = (events_sparse[idx] != 0)
        events = events_sparse[idx, events_bool]

        if events.size:
            duration_frames = np.abs(history[idx, :-1][events_bool])
            direction_perm = jump_type[idx, events_bool]
            mol_idx_chunk = mols_idx_chunk[idx]
            end_frame = np.abs(events)
            start_frame = end_frame - duration_frames + 1
            permeations_chunk.append([mol_idx_chunk, start_frame, end_frame, duration_frames, direction_perm])

    return permeations_chunk

# ---------------------------------------------------------------------------------------------------------

if __name__ == "__main__":

    welcome()

    # Check arguments
    args = parse_and_validate_args()

    top_path = args.topology
    traj_path = args.trajectory
    output_prefix = args.output
    # Check available memory
    n_processes = cpu_count() - 2 if cpu_count() > 2 else 1
    traj_file_size = os.path.getsize(traj_path)

    if not have_enough_memory(n_processes, traj_file_size):
        exit(1)

    # --------------------------------------------------------------------------
    # BEGIN PROGRAM
    # --------------------------------------------------------------------------

    # PARAMETERS
    MAX_DISPLACEMENT = 30.0         # Å, threshold for PBC jump detection
    mols_selection_str = "resname TIP3 SOL TP3 and name OW OH2"
    membrane_selection_str = "name P"
    
    start_time_prog = time.time() # Monitor the total execution time

    # SYSTEM INFO
    u = mda.Universe(top_path, traj_path, in_memory=False)
    n_frames = u.trajectory.n_frames
    mols_sel = u.select_atoms(mols_selection_str)
    n_mols_atoms = mols_sel.n_atoms


    # --------------------------------------------------------------------------
    # STEP 1: EXTRACT MEMBRANE BOUNDARIES (UPPER AND LOWER LEAFLETS) AND MOLECULES
    # Z POSITIONS
    # Multiprocess by chunks of frames

    chunk_size = n_frames // n_processes
    args_list = []
    for i in range(n_processes):
        start_idx = i * chunk_size
        end_idx = n_frames if i == n_processes - 1 else (i + 1) * chunk_size
        args_list.append((start_idx, end_idx, traj_path, top_path, membrane_selection_str, mols_selection_str, MAX_DISPLACEMENT))

    print(f"Extracting membrane and molecule coordinates on {n_frames} frames using {n_processes} processes...")
    start_time_step = time.time()

    with Pool(processes=n_processes) as pool:
        results = pool.map(process_z_coords_extraction, args_list)

    start_idx, end_idx, lower_boundaries_chunks, upper_boundaries_chunks, mols_z_chunks = zip(*results)
    upper_boundaries = np.concatenate(upper_boundaries_chunks)[1:]
    lower_boundaries = np.concatenate(lower_boundaries_chunks)[1:]
    mols_z = np.concatenate(mols_z_chunks, axis=1)
    mols_z_prev = mols_z[:, :-1]
    mols_z_curr = mols_z[:, 1:]

    elapsed_time = int(time.time() - start_time_step)
    print(f"Extraction completed in {convert_seconds_to_hours_min(elapsed_time)}\n")


    # --------------------------------------------------------------------------
    # STEP 2: FIND PERMEATION EVENTS
    # Multiprocess by chunks of molecules

    chunk_size = n_mols_atoms // n_processes
    args_list = []
    for i in range(n_processes):
        start_idx = i * chunk_size
        end_idx = n_mols_atoms if i == n_processes - 1 else (i + 1) * chunk_size
        args_list.append((start_idx, end_idx, upper_boundaries,
                        lower_boundaries, mols_z_prev, mols_z_curr, MAX_DISPLACEMENT))
        
    print(f"Finding permeation events of {n_mols_atoms} molecules using {n_processes} processes...")
    start_time_step = time.time()

    with Pool(processes=n_processes) as pool:
        results = pool.map(process_permeations_search, args_list)

    merged_results = [event for chunk in results if chunk for event in chunk]

    elapsed_time = int(time.time() - start_time_step)
    print(f"Search completed in {convert_seconds_to_hours_min(elapsed_time)}")
    print(f"Found {len(merged_results)} permeation events.\n")


    # --------------------------------------------------------------------------
    # STEP 3 GENERATE OUTPUTS

    print(f"Writing outputs...")

    # Save results to CSV
    output_file = f"{output_prefix}_permeation_events.csv"

    df = pd.DataFrame(merged_results, columns=["index_in_selection", "start_frame", "end_frame", "duration_frames", "direction"])
    df["resid_in_system"] = mols_sel.atoms.resids[df["index_in_selection"]]

    # index_in_selection: index of the molecule in the selection
    # resid_in_system: residue ID of the molecule in the original system
    # start_frame, end_frame: first and last frame of the permeation event
    # duration_frames: number of frames the permeation took
    # direction: -1 for lower to upper, +1 for upper to lower
    df = df[["index_in_selection", "resid_in_system", "start_frame", "end_frame", "duration_frames", "direction"]]
    df.to_csv(output_file, index=False)

    print(f"Permeations events written in {output_file}\n")

    # --------------------------------------------------------------------------
    # End program
    # --------------------------------------------------------------------------

    elapsed_time = int(time.time() - start_time_prog)
    print(f"Program terminated successfully in {convert_seconds_to_hours_min(elapsed_time)}")
    
