# Code for "Quantifying Variability in Lagrangian Particle Dispersal in Ocean Ensemble Simulations: an Information Theory Approach"

This repository contains code for the Lagrangian simulations and the analysis for the manuscript submitted to Nonlinear Processes in Geophysics.

## Getting Started

To reproduce the simulations and the analysis, you will need to clone this repository to your machine. You can do this by running the following command in your terminal:

`git clone https://github.com/OceanParcels/NEMO_Ensemble_Lagrangian_Analysis.git`

To able to run the scripts we recommend creating a conda environment with the package versions required. You can create a new Conda environment from an [environment.yml](environment.yml) file using:

`conda env create -f my_environment.yml`

## File Structure

To reproduce the results it is necessary to execute the scripts in the order they are listed below:

#### [Simulations](simulations/)

First execute the simulations scripts. The simulations scripts consist of two types, spatial and temporal releases. Each type has its own script. We submit the scripts with `.sh` files. 

1. [`ensemble_Member_spatial.py`](simulation/ensemble_Member_spatial.py): Script for spatial release simualtions.
2. [`ensemble_Member_temporal.py`](simulation/ensemble_Member_temporal.py): Script for temporal release simualtions.
3. [`submit_spatial_jobs.sh`](submit_spatial_jobs.sh): submits the spatial simulations for the 50 members of the ensemble, for three cases: $0.1^\circ$, $1.0^\circ$, and $2.0^\circ$. It takes approx 8 hours per simulation, and there are 150 to be computed (be careful with the 120 hour limit).
4. [`submit_temporal_jobs.sh`](submit_temporal_jobs.sh): submits the temporal simulations for the 50 members of the ensemble, for three cases: 4, 12, and 20 weeks. It takes approx 8 hours per simulation, and there are 150 to be computed (be careful with the 120 hour limit).

#### [Analysis](analysis/)

Once the simulations are computed, we can perform the analysis. The first script that you need to run is `generate_hexbin_grids.py`. This generates the grid used to bin the particle positions. After this you can run the `single_memeber_...py` scripts, followed by `mixture_distributions_...py`, and then the `connectivity_...py` scripts.

1. [`single_member_statistics_spatial.py`](analysis\single_member_statistics_spatial.py):
2. [`single_member_statistics_temp.py`](analysis\single_member_statistics_temp.py): 
2. [`mixture_distributions_spatial.py`](analysis/mixture_distributions_spatial.py): 
3. [`mixture_distributions_temp.py`](analysis/mixture_distributions_spatial.py): 
4. [``](analysis/): 
5. [``](analysis/): 
6. [``](analysis/):

#### [Plot scripts](plots_scripts/)

#### [Functions](functions/)

#### [Data](data/)

#### [Figures](figs/)

## Overview Contents

```
NEMO_ensemble/
├── README.md
├── analysis/
├── check_corrupted_files/
├── data/
├── figs/
├── functions/
├── plot_scripts/
└── simulations/
```

The `NEMO_ensemble` folder contains the following files and directories:

- `README.md`: This file provides an overview of the repository and its contents.

- `analysis/`: Scripts to compute the analysis.

- `check_corrupted_files/`: Scripts that check if the downloaded fields are corrupted.

- `data/`: This directory stores the data files used in the NEMO Ensemble simulation.

- `figs/`: This directory contains the figures produced by the `analysis` or `notebooks`.

- `functions/`: Functions called as modules in several scripts.

- `plot_scripts`: Scripts to reproduce the plots in the manuscript.

- `simulations/`: This directory contains the simulations scripts. 