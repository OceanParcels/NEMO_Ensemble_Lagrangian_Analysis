#!/usr/bin/env python
# coding: utf-8
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

np.random.seed(43)
Subset = 20


# Load the hexbin_grid for the domain
with open('../data/hexgrid_no_coast.pkl', 'rb') as f:
    hexbin_grid = pickle.load(f)
    
hexbin_grid = hexfunc.hexGrid(hexbin_grid, h3_res=3)

###### Calculate for all memebers and delta_rs ####
# delta_r_ranges = np.linspace(0.1, 1, 10)
members = np.arange(2, 51)
N_subsets = 50

location = 'Cape_Hatteras'
delta_r = 0.1

for delta_r in [0.1, 1]:
    for k in range(1, N_subsets+1):
        member = 1
        
        path = f"/storage/shared/oceanparcels/output_data/data_Claudio/NEMO_Ensemble/{location}/spatial/dr_{delta_r*100:03.0f}/{location}_dr{delta_r*100:03.0f}_m{member:03d}.zarr"
        pset_members = xr.open_zarr(path)
        obs_range = range(len(pset_members.obs)) # Number of time steps in the observation period

        pset_members = pset_members.isel(trajectory=np.random.choice(pset_members.trajectory, 20, replace=False))

        print(f"Subset:{k} delta_r: {delta_r}")

        for member in tqdm(members):

            path = f"/storage/shared/oceanparcels/output_data/data_Claudio/NEMO_Ensemble/{location}/spatial/dr_{delta_r*100:03.0f}/{location}_dr{delta_r*100:03.0f}_m{member:03d}.zarr"
            pset = xr.open_zarr(path)
            
            
            pset = pset.isel(trajectory=np.random.choice(pset.trajectory, Subset, replace=False))
            
            pset_members = xr.concat([pset_members, pset], dim='trajectory')


        print("length pset_members: ", len(pset_members.trajectory))


        P_m, Ent_m = calculate_probability_and_entropy(pset_members, hexbin_grid, entropy)
        DF_m = create_dataframe(P_m, Ent_m, hexbin_grid.hexint, obs_range)
        save_path = f"/storage/shared/oceanparcels/output_data/data_Claudio/NEMO_Ensemble/analysis/prob_distribution/{location}_all/P_dr{delta_r*100:03.0f}_all_s{k}.nc"
        DF_m.to_netcdf(save_path)

