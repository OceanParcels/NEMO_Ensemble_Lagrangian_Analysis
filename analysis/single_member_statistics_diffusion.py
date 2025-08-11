#!/usr/bin/env python
# %% coding: utf-8
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
    number_particles_binned = np.zeros(obs_length)

    lons, lats = pset['lon'][:, :].values, pset['lat'][:, :].values

    for t in range(obs_length):
        _probability = hexbin_grid.count_2d(
            lons[:, t], lats[:, t], normalize=False)
        probability_set[:, t] = _probability
        entropy_set[t] = entropy_function(probability_set[:, t]/np.nansum(_probability))
        number_particles_binned[t] = np.nansum(_probability)

    return probability_set, entropy_set, number_particles_binned


def create_dataframe(probability_set, entropy_set, number_particles_set, hexints, time_range):
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
            ),
            'number_particles_binned': xr.DataArray(
                number_particles_set,
                dims=['time'],
                coords={
                    'time': time_range
                },
                attrs={
                    'description': 'Number of particles binned for each time step',
                    'units': 'count'
                }
            )
        }
    )

    return ds


# %%
location = 'Cape_Hatteras'
member = 1  # memeber
K_h = 1000  # Standard deviation od initial dispersion

path = f"/Volumes/Claudio SSD/Ensemble_article_data/simulations/diff_Kh_{K_h:01d}/{location}_diff_Kh_{K_h:01d}_m{member:03d}.zarr"
# path = f"/storage/shared/oceanparcels/output_data/data_Claudio/NEMO_Ensemble/{location}/diff_long/diff_Kh_{K_h:01d}/{location}_diff_Kh_{K_h:01d}_m{member:03d}.zarr"
pset = xr.open_zarr(path)

obs_range = pset.obs.values  # Number of time steps in the observation period

# Load the hexbin_grid for the domain
with open('../data/hexgrid_no_coast_h3.pkl', 'rb') as f:
    hexbin_grid = pickle.load(f)

hexbin_grid = hexfunc.hexGrid(hexbin_grid, h3_res=3)

mask_file = '/Volumes/Claudio SSD/Ensemble_article_data/NATL025-CJMenobs01_byte_mask.nc'
# mask_file = '/storage/shared/oceanparcels/input_data/NEMO_Ensemble/GRID/NATL025-CJMenobs01_byte_mask.nc'
mask = xr.open_dataset(mask_file, decode_times=False)

tmask = mask['tmask'][0,0].values
mask_lons = mask['nav_lon'][0, :].values
mask_lats = mask['nav_lat'][:, 0].values

# %%##### Calculate for all memebers and delta_rs ####
K_h_ranges = [1000]  # np.linspace(0.1, 1, 10)

keep_all_traj = False # If True, keep all trajectories, if False, subsample
subsample = 7500 # Number of particles to subsample if keep_all_traj is False

members = [4]  # np.arange(1, 51)

for member in members:
    for K_h in K_h_ranges:
        print(f"\U0001F914 Member: {member:03d},  K_h: {K_h}")
        path = f"/Volumes/Claudio SSD/Ensemble_article_data/simulations/diff_Kh_{K_h:01d}/{location}_diff_Kh_{K_h:01d}_m{member:03d}.zarr"
        # path = f"/storage/shared/oceanparcels/output_data/data_Claudio/NEMO_Ensemble/{location}/diff_long/diff_Kh_{K_h:01d}/{location}_diff_Kh_{K_h:01d}_m{member:03d}.zarr"
        pset = xr.open_zarr(path)
        # pset = xr.open_dataset(path)

        if K_h == 1000:

            full_trajectories = []

            # remove particles othat go inland
            lon_arr = pset['lon'].values
            lat_arr = pset['lat'].values

            for p in tqdm(range(pset.sizes['trajectory'])):

                lon_idx = np.digitize(
                    pset.lon[p, :].dropna(dim='obs'), mask_lons)

                lat_idx = np.digitize(
                    pset.lat[p, :].dropna(dim='obs'), mask_lats)
                
                print(lon_idx)

                if tmask.shape[0] in lat_idx:
                    # if lat_idx == shape[0], then find where lat_idx is == shape[0]. Remove those values from lat_idx and lon_idx
                    _lat_idx = lat_idx[lat_idx < tmask.shape[0]]
                    lon_idx = lon_idx[lat_idx < tmask.shape[0]]
                    lat_idx = _lat_idx
                if tmask.shape[1] in lon_idx:
                    print("AYE")
                    _lon_idx = lon_idx[lon_idx < tmask.shape[1]]
                    lat_idx = lat_idx[lon_idx < tmask.shape[1]]
                    lon_idx = _lon_idx

                print(lon_idx)

                tmask_values = tmask[lat_idx, lon_idx]
                idx = np.where(tmask_values == 0)[0] # Number of time steps outside the mask
                if len(idx) > 0:
                    cutoff = idx[0]  # First index outside the mask
    
                    lon_arr[p, cutoff:] = np.nan
                    lat_arr[p, cutoff:] = np.nan
                    
                elif len(idx) == 0:
                    full_trajectories.append(p)  # Store the last index if all are inside the mask

           
            full_trajectories = np.array(full_trajectories)
            #save the full trajectories to npz
            np.savez(f'../data/full_traj_K_h{K_h:01d}/full_trajectories_m{member:03d}_K_h{K_h:01d}.npz', full_trajectories=full_trajectories)

            pset = pset.isel(trajectory=full_trajectories)

        
            try :
                pset = pset.isel(trajectory=np.random.choice(
                    pset.sizes['trajectory'], subsample, replace=False))  # replace=False, no repeated trajectories
            except ValueError as e:
                print(f"Error: {e}. Subsampling {pset.sizes['trajectory']} trajectories, but requested {subsample} subsamples.")
                print("Setting subsample to the number of available trajectories.")
                # Subsample the trajectories again with the correct number

                pset = pset.isel(trajectory=np.random.choice(
                    pset.sizes['trajectory'], pset.sizes['trajectory'], replace=False)) # replace=False, no repeated trajectories
                

        # Calculate the probability and entropy
        P_m, Ent_m, Np_m = calculate_probability_and_entropy(
            pset, hexbin_grid, entropy)
        DF_m = create_dataframe(P_m, Ent_m, Np_m, hexbin_grid.hexint, obs_range)
        save_path = f"/Volumes/Claudio SSD/Ensemble_article_data/analysis/prob_distribution/{location}_diffusion_long/P_diff_Kh_{K_h:01d}_m{member:03d}.nc"
        # save_path = f"/storage/shared/oceanparcels/output_data/data_Claudio/NEMO_Ensemble/analysis/prob_distribution/{location}_diffusion_long/P_diff_Kh_{K_h:01d}_m{member:03d}{subsample_str}.nc"
        DF_m.to_netcdf(save_path)

# %%
