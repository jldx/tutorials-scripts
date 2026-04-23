# Perm MD Counter

**A Python module for detecting water permeation events in molecular dynamics (MD) trajectories of membrane transporters.**

This tool identifies and counts water permeation events through membrane transporters from MD simulation trajectories.

---

## Key Features

- **Efficient detection** of water permeation events in MD trajectories
- **Compatibility** with standard MD file formats (e.g., GROMACS, AMBER)
- **Preprocessing requirement**: Trajectories must be centered in the simulation box

---

## Installation

### Prerequisites
- [Conda](https://docs.conda.io/en/latest/) (recommended for environment management)
- Python 3.12

### Setup
1. Clone this repository:
   ```bash
   git clone https://github.com/ahardiag/perm-md-count.git
   cd perm-md-count
   ```

2. Create and activate the conda environment:
    ```bash
    conda env create -f perm_lip.yml
    conda activate perm_lip
    ```

## Usage

```bash

python3 perm_lip_water.py -r <topology_file> -f <trajectory_file> -o <output_folder> -x <output_prefix> -p -t --print_limits -s <number_of_water_subsets>   
```

## Output

* A log file with the list of processed files
* A csv file with the limits of the membrane upper and lower leaflets
* A csv files with the identified molecules, permeation description

## Licence

MIT Licence.

## Citation

Original paper:
> Hardiagon, A.; Murail, S.; Huang, L.-B.; van der Lee, A.; Sterpone, F.; Barboiu, M.; Baaden, M. Molecular Dynamics Simulations Reveal Statistics and Microscopic Mechanisms of Water Permeation in Membrane-Embedded Artificial Water Channel Nanoconstructs. J. Chem. Phys. 2021, 154 (18), 184102. https://doi.org/10.1063/5.0044360.



