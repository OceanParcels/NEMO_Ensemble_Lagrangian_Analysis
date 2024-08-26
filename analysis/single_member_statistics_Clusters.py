#!/usr/bin/env python
# coding: utf-8
#%% 
import numpy as np
import xarray as xr
from tqdm import tqdm
import pickle
import sys
sys.path.append('../functions')
import hexbin_functions as hexfunc

def entropy(Pdf):
    # Shannon entropy
    # Pdf = Pdf / np.nansum(Pdf)  # Normalize Pdf to sum to 1, ignoring NaNs
    # Replace zeros with a very small number to avoid log(0)
    Pdf_safe = np.where(Pdf > 0, Pdf, np.finfo(float).eps)
    return -np.nansum(Pdf_safe * np.log2(Pdf_safe))


def calculate_probability_and_entropy(pset, hexbin_grid, entropy_function):
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

    probability_set = np.zeros((n_hex, obs_length))
    entropy_set = np.zeros(obs_length)

    lons, lats = pset['lon'][:, :].values, pset['lat'][:, :].values

    for t in range(obs_length):
        probability_set[:, t] = hexbin_grid.count_2d(lons[:, t], lats[:, t], normalize=True)
        entropy_set[t] = entropy_function(probability_set[:, t])


    return probability_set, entropy_set


def create_dataframe(probability_set, entropy_set, hexints, time_range):
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
                probability_set,
                dims=['hexint', 'time'],
                coords={
                    'hexint': hexints,
                    'time': time_range
                },
                attrs={
                    'description': 'Probability of occurrence for each time step, hexagonal bin, and observation time',
                    'units': 'probability'
                }
            ),
            'entropy': xr.DataArray(
                entropy_set,
                dims=['time'],
                coords={
                    'time': time_range
                },
                attrs={
                    'description': 'Entropy values for each time step and observation time',
                    'units': 'bits'
                }
            )
        }
    )
    
    return ds

#%%
member = 1 # memeber
cluster = 0 # cluster

path = f"/storage/shared/oceanparcels/output_data/data_Claudio/NEMO_Ensemble/Clusters/Cluster_{cluster}/cluster{cluster}_m{member:03d}.zarr"
pset = xr.open_zarr(path)
obs_range = pset.obs.values # Number of time steps in the observation period

# Load the hexbin_grid for the domain
hex_res = 3
with open(f'../data/hexgrid_no_coast_h{hex_res}.pkl', 'rb') as f:
    hexbin_grid = pickle.load(f)
    
hexbin_grid = hexfunc.hexGrid(hexbin_grid, h3_res=hex_res)
members = np.arange(1, 51)

#%%
for member in tqdm(members):

    print(f"Member: {member:03d}, Cluster: {cluster}")
    
    path = f"/storage/shared/oceanparcels/output_data/data_Claudio/NEMO_Ensemble/Clusters/Cluster_{cluster}/cluster{cluster}_m{member:03d}.zarr"
    pset = xr.open_zarr(path)
    print("--", len(pset.trajectory))
    
    P_m, Ent_m = calculate_probability_and_entropy(pset, hexbin_grid, entropy)
    DF_m = create_dataframe(P_m, Ent_m, hexbin_grid.hexint, obs_range)
    save_path = f"/storage/shared/oceanparcels/output_data/data_Claudio/NEMO_Ensemble/analysis/prob_distribution/Clusters/Cluster_{cluster}/P_cluster_{cluster}_m{member:03d}.nc"
    DF_m.to_netcdf(save_path)

# %%
