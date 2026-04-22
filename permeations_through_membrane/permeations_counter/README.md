# Permeation Counter

**A Python module for detecting water permeation events in molecular dynamics (MD) trajectories of membrane transporters.**

This tool identifies and counts water permeation events through membrane transporters from MD simulation trajectories. It is optimized for efficiency and parallel processing, building on the original [perm-md-count](https://github.com/ahardiag/perm-md-count).

---

## Key Features

- **Efficient detection** of water permeation events in MD trajectories
- **Parallel processing** for faster analysis of large datasets
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
   git clone https://github.com/jldx/permeation-counter.git
   cd permeation-counter
   ```

2. Create and activate the conda environment:
    ```bash
    conda env create -f environment.yml
    conda activate permeations_counter
    ```

## Usage

```bash
python3 permeations_counter.py \
    --topology <topology_file> \
    --trajectory <trajectory_file> \
    --output <output_prefix>
```

## Output

* A log file with detected permeation events and timestamps
* A summary file with identified molecules, permeation description

## Licence

MIT Licence.

## Citation

Original paper:
> Hardiagon, A.; Murail, S.; Huang, L.-B.; van der Lee, A.; Sterpone, F.; Barboiu, M.; Baaden, M. Molecular Dynamics Simulations Reveal Statistics and Microscopic Mechanisms of Water Permeation in Membrane-Embedded Artificial Water Channel Nanoconstructs. J. Chem. Phys. 2021, 154 (18), 184102. https://doi.org/10.1063/5.0044360.



