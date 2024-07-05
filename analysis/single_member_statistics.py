#!/usr/bin/env python
# coding: utf-8
import numpy as np
import xarray as xr
from tqdm import tqdm
import pickle
import hexbin_functions as hexfunc

def entropy(Pdf):
    # Pdf = Pdf / np.nansum(Pdf)  # Normalize Pdf to sum to 1, ignoring NaNs
    # Replace zeros with a very small number to avoid log(0)
    Pdf_safe = np.where(Pdf > 0, Pdf, np.finfo(float).eps)
    return -np.nansum(Pdf_safe * np.log(Pdf_safe))


def calculate_probability_and_entropy(pset, hexbin_grid, subgroups, entropy_function):
    """
    Calculates probability and entropy for particle sets over a hexagonal grid.

    Parameters
    ----------
    pset : xarray.Dataset
        Particle dataset containing longitude and latitude variables.
    hexbin_grid : hexGrid object
        An object representing the hexagonal grid, with a method count_2d to count particles within each hexbin.
    subgroups : dict
        Dictionary mapping t_gap values to indices of particles released at different times.
    entropy_function : function
        Function to calculate entropy given a probability distribution.

    Returns
    -------
    probability_sets : dict
        Dictionary of probability arrays for each t_gap, with dimensions (n_hex, obs_length).
    entropy_sets : dict
        Dictionary of entropy values for each t_gap, with dimension (obs_length).
    """
    obs_length = len(pset.obs)
    n_hex = hexbin_grid.n_hex
    probability_sets = {}
    entropy_sets = {}

    for t_gap in subgroups.keys():
        if t_gap == 0:
            sub_pset = pset
        else:
            sub_pset = pset.isel(trajectory=subgroups[t_gap])

        probability_subset = np.zeros((n_hex, obs_length))
        entropy_subset = np.zeros(obs_length)

        lons, lats = sub_pset['lon'][:, :].values, sub_pset['lat'][:, :].values

        for t in range(obs_length):
            probability_subset[:, t] = hexbin_grid.count_2d(lons[:, t], lats[:, t], normalize=True)
            entropy_subset[t] = entropy_function(probability_subset[:, t])

        probability_sets[t_gap] = probability_subset
        entropy_sets[t_gap] = entropy_subset

    return probability_sets, entropy_sets


def create_dataframe(probability_sets, entropy_sets, hexints):
    """
    Creates xarray Dataframe containing the probability and entropy data.

    Parameters
    ----------
    probability_sets : dict
        Dictionary containing probability data arrays for each delta_t.
    entropy_sets : dict
        Dictionary containing entropy data arrays for each delta_t.
    hexints : list
        List of hexagonal bin indices.
    obs_length : int
        The length of the observation period.
    filename : str
        The filename to save the NetCDF file.
    
    Returns
    -------
    ds : xarray.Dataset
        The dataset containing the probability and entropy data.
    """

    ds = xr.Dataset(
        {
            'probability': xr.DataArray(
                np.array([probability_sets[i] for i in range(0, 15)]),
                dims=['delta_t', 'hexint', 'time'],
                coords={
                    'delta_t': range(0, 15),
                    'hexint': hexints,
                    'time': range(probability_sets[0].shape[1])
                },
                attrs={
                    'description': 'Probability of occurrence for each time step, hexagonal bin, and observation time',
                    'units': 'probability'
                }
            ),
            'entropy': xr.DataArray(
                np.array([entropy_sets[i] for i in range(0, 15)]),
                dims=['delta_t', 'time'],
                coords={
                    'delta_t': range(0, 15),
                    'time': range(probability_sets[0].shape[1])
                },
                attrs={
                    'description': 'Entropy values for each time step and observation time',
                    'units': 'nats'
                }
            )
        }
    )
    
    return ds


location = 'Cape_Hatteras'
member = 50 # memeber
std = 0.1 # Standard deviation od initial dispersion

file_path = f"/storage/shared/oceanparcels/output_data/data_Claudio/NEMO_Ensemble/{location}/std_{std*100:03.0f}/{location}_std{std*100:03.0f}_m{member:03d}.zarr"
pset = xr.open_zarr(file_path)

obs_length = len(pset.obs)
# Create subgroups of particles released at different times
subgroups = {}

max_gap = 15
timesteps = np.linspace(1, 745, 745, dtype=int)
max_releases = timesteps[0::max_gap].shape[0]

set = np.linspace(0, 99, 100, dtype=int)

for gap in range(0, max_gap):
    
    if gap == 0:
        indexes = pset.trajectory.values # groups with time gap 0 are the full set of particles
    else:
        set_start = timesteps[0::gap][:max_releases]*100
        indexes = []
    
        for j in set_start:
            indexes = indexes + list(j+set)
    
    subgroups[gap] = np.array(indexes)


# Load the hexbin_grid for the domain
with open('../data/hexgrid_no_coast.pkl', 'rb') as f:
    hexbin_grid = pickle.load(f)
    
hexbin_grid = hexfunc.hexGrid(hexbin_grid, h3_res=3)


###### Calculate for all memebers and STDs ####
std_ranges = np.linspace(1, 20, 20)/100
location = 'Cape_Hatteras'
member = 1 # memeber
# std = 0.01 # Standard deviation od initial dispersion

for member in tqdm([1, 2, 47, 48, 49, 50]):
    for std in std_ranges:
        print(f"\U0001F914 Member: {member:03d},  std: {std:03.0f}")
        path = f"/storage/shared/oceanparcels/output_data/data_Claudio/NEMO_Ensemble/{location}/std_{std*100:03.0f}/{location}_std{std*100:03.0f}_m{member:03d}.zarr"

        pset = xr.open_zarr(path)
        P_m, Ent_m = calculate_probability_and_entropy(pset, hexbin_grid, subgroups, entropy)
        DF_m = create_dataframe(P_m, Ent_m, hexbin_grid.hexint)
        save_path = f"/storage/shared/oceanparcels/output_data/data_Claudio/NEMO_Ensemble/analysis/prob_distribution/{location}/P_std{std*100:03.0f}_m{member:03d}"
        DF_m.to_zarr(save_path)
        # print(f"Member {member:03d} saved at ../{location}/P_std{std*100:03.0f}_m{member:03d}")
        

