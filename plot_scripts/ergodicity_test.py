#%%!/usr/bin/env python
# coding: utf-8
import numpy as np
import xarray as xr
from tqdm import tqdm
import pickle
import sys
sys.path.append('../functions')
import hexbin_functions as hexfunc
import matplotlib.pyplot as plt

def entropy(Pdf):
    # Shannon entropy
    # Pdf = Pdf / np.nansum(Pdf)  # Normalize Pdf to sum to 1, ignoring NaNs
    # Replace zeros with a very small number to avoid log(0)
    Pdf_safe = np.where(Pdf > 0, Pdf, np.finfo(float).eps)
    return -np.nansum(Pdf_safe * np.log2(Pdf_safe))


def calculate_probability_one_trajectory(pset, hexbin_grid, entropy_function, traj=0):
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

    probability_set = np.zeros(n_hex)
    

    lons, lats = pset['lon'][traj, :].values, pset['lat'][traj, :].values

    probability_set[:] = hexbin_grid.count_2d(lons, lats, normalize=True)
    entropy = entropy_function(probability_set)

    return probability_set, entropy


def calculate_probability_and_entropy(pset, hexbin_grid, entropy_function, time, n_traj=100):
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

    probability_set = np.zeros(n_hex)
    n_particles = slice(n_traj)

    lons, lats = pset['lon'][n_particles, time].values, pset['lat'][n_particles, time].values
    
    probability_set[:] = hexbin_grid.count_2d(lons, lats, normalize=True)
    entropy = entropy_function(probability_set[:])

    return probability_set, entropy

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

location = 'Cape_Hatteras'
member = 1 # memeber
delta_r = 0.1 # Standard deviation od initial dispersion

path = f"/storage/shared/oceanparcels/output_data/data_Claudio/NEMO_Ensemble/{location}/spatial_long/dr_{delta_r*100:03.0f}/{location}_dr{delta_r*100:03.0f}_m{member:03d}.zarr"
pset = xr.open_zarr(path)

obs_range = pset.obs.values # Number of time steps in the observation period

# Load the hexbin_grid for the domain
with open('../data/hexgrid_no_coast.pkl', 'rb') as f:
    hexbin_grid = pickle.load(f)
    
hexbin_grid = hexfunc.hexGrid(hexbin_grid, h3_res=3)

particulos = 200
single_H = np.zeros(particulos)

for i in tqdm(range(particulos)):
    _, single_H[i] = calculate_probability_one_trajectory(pset, hexbin_grid, entropy, traj=i)

_, H_space = calculate_probability_and_entropy(pset, hexbin_grid, entropy, time=2000, n_traj=2000)

#%%
plt.hist(single_H, bins=20)
plt.axvline(H_space, color='r')
plt.axvline(np.mean(single_H), color='k', )
plt.text(H_space+0.1, 3, f'{H_space:0.1f} bits', rotation=90)
plt.text(np.mean(single_H)+0.1, 3, f'{np.mean(single_H):0.1f} bits', rotation=90)

print('Entropy one traj:', H)
print('Entropy 2189 particles as 2100 days:', H_space)

plt.xlabel('Entropy (bits)')
plt.ylabel('Frequency')
plt.title('Entropy distribution for 200 particles')

plt.savefig('../figs/ergodicity?nope.png', dpi=300)

# %%
