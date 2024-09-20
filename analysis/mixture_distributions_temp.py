#!/usr/bin/env python
# coding: utf-8
import numpy as np
import xarray as xr
from tqdm import tqdm
import pickle
import sys
sys.path.append('../functions')
import hexbin_functions as hexfunc
import concurrent.futures

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
    

    def process_time_step(t):
        prob = hexbin_grid.count_2d(pset['lon'][:, t], pset['lat'][:, t], normalize=True)
        ent = entropy_function(prob)
        return t, prob, ent

    with concurrent.futures.ThreadPoolExecutor() as executor:
        futures = [executor.submit(process_time_step, t) for t in range(obs_length)]
        for future in tqdm(concurrent.futures.as_completed(futures)):
            t, prob, ent = future.result()
            probability_set[:, t] = prob
            entropy_set[t] = ent

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


# Load the hexbin_grid for the domain
with open('../data/hexgrid_no_coast.pkl', 'rb') as f:
    hexbin_grid = pickle.load(f)
    
hexbin_grid = hexfunc.hexGrid(hexbin_grid, h3_res=3)

###### Calculate for all memebers and delta_rs ####
members = np.arange(2, 51)
N_subsets = 50

location = 'Cape_Hatteras'
week = 4 # Number of weeks
subset_particles = 148


def process_member(member, week, location, subset_particles):
    path = f"/storage/shared/oceanparcels/output_data/data_Claudio/NEMO_Ensemble/{location}/temporal_long/W_{week:01d}/{location}_W{week:01d}_m{member:03d}.zarr"
    pset = xr.open_zarr(path)
    pset = pset.isel(trajectory=np.random.choice(pset.trajectory, subset_particles, replace=False))
    return pset

for k in range(1, N_subsets+1):
    member = 1

    path = f"/storage/shared/oceanparcels/output_data/data_Claudio/NEMO_Ensemble/{location}/temporal_long/W_{week:01d}/{location}_W{week:01d}_m{member:03d}.zarr"
    pset_members = xr.open_zarr(path)
    obs_range = range(len(pset_members.obs)) # Number of time steps in the observation period

    pset_members = pset_members.isel(trajectory=np.random.choice(pset_members.trajectory, subset_particles, replace=False))

    print(f"Subset:{k} period: {week} weeks.")

    with concurrent.futures.ThreadPoolExecutor() as executor:
        futures = [executor.submit(process_member, member, week, location, subset_particles) for member in members]
        for future in tqdm(concurrent.futures.as_completed(futures), total=len(members)):
            pset = future.result()
            pset_members = xr.concat([pset_members, pset], dim='trajectory')

    print("length pset_members: ", len(pset_members.trajectory))

    P_m, Ent_m = calculate_probability_and_entropy(pset_members, hexbin_grid, entropy)
    DF_m = create_dataframe(P_m, Ent_m, hexbin_grid.hexint, obs_range)
    save_path = f"/storage/shared/oceanparcels/output_data/data_Claudio/NEMO_Ensemble/analysis/prob_distribution/{location}_all_long/P_W{week:02d}_all_s{k:03d}.nc"
    DF_m.to_netcdf(save_path)

