# %% Load the packages
import numpy as np
import xarray as xr
from tqdm import tqdm
import pandas as pd
import pickle
import os

location = "Cape_Hatteras"
path = "/Volumes/Claudio SSD/Ensemble_article_data/"

Latitude_limit = None
Longitude_limit = -40

mask_file = '/Volumes/Claudio SSD/Ensemble_article_data/NATL025-CJMenobs01_byte_mask.nc'
# mask_file = '/storage/shared/oceanparcels/input_data/NEMO_Ensemble/GRID/NATL025-CJMenobs01_byte_mask.nc'
mask = xr.open_dataset(mask_file, decode_times=False)

tmask = mask['tmask'][0,0].values
mask_lons = mask['nav_lon'][0, :].values
mask_lats = mask['nav_lat'][:, 0].values

# %% Spatial analysis
distributions = {}

total_members = 50
K_h = 1000
subsample = 7500

for member in tqdm(range(1, total_members + 1)):
    print(f"Member: {member:03d},  K_h: {K_h}")
    file_path = path + f"simulations/diff_Kh_{K_h:01d}/{location}_diff_Kh_{K_h:01d}_m{member:03d}.zarr"
    
    pset = xr.open_zarr(file_path)

        
    full_trajectories = np.load(f'../data/full_traj_K_h{K_h:01d}/full_trajectories_m{member:03d}_K_h{K_h:01d}.npz')['full_trajectories']
    pset = pset.isel(trajectory=full_trajectories)

    # Subsample the trajectories
    try :
        pset = pset.isel(trajectory=np.random.choice(
            pset.sizes['trajectory'], subsample, replace=False))  # replace=False, no repeated trajectories
    except ValueError as e:
        print(f"Error: {e}. Subsampling {pset.sizes['trajectory']} trajectories, but requested {subsample} subsamples.")
        print("Setting subsample to the number of available trajectories.")
        # Subsample the trajectories again with the correct number

        pset = pset.isel(trajectory=np.random.choice(
            pset.sizes['trajectory'], pset.sizes['trajectory'], replace=False)) # replace=False, no repeated trajectories


    N_particles = len(pset.trajectory)

    if Latitude_limit is not None:
        lats = pset.lat.load().values
        p_index, t_index = np.where(lats[:, :] > Latitude_limit)
    elif Longitude_limit is not None:
        lons = pset.lon.load().values
        p_index, t_index = np.where(lons[:, :] > Longitude_limit)
    
    
    subpolar_traj = np.unique(p_index) # it's no longer a subpolar but I keep the name
    drift_time = []

    if len(subpolar_traj) > 0:
        for i in subpolar_traj:
            idx_t = np.where(p_index == i)[0][0]
            drift_time.append(t_index[idx_t])
        
        drift_time = np.array(drift_time)
        
        depths = pset.z.load().values
        depths = depths[subpolar_traj, drift_time]

        distributions["member"] = member
        distributions["drift_time"] = drift_time
        distributions["depths"] = depths
        distributions["trajectory"] = np.unique(p_index)
        
        # SAVE DISTRIBUTIONS in a pickle file
        if Latitude_limit is not None:
            save_path = path + f"analysis/connectivity/Kh_{K_h:01d}_{Latitude_limit}N/Distributions_Kh_{K_h:01d}_m{member:03d}.pkl"
        elif Longitude_limit is not None:    
            save_path = path + f"analysis/connectivity/Kh_{K_h:01d}_{abs(Longitude_limit)}W/Distributions_Kh_{K_h:01d}_m{member:03d}.pkl"

        # Create the directory if it does not exist
        os.makedirs(os.path.dirname(save_path), exist_ok=True)
            
        with open(save_path, "wb") as f:
            pickle.dump(distributions, f)
        
    else:
        print(f"--EMPTY--")
    
#%% Build the Pandas Dataframes from the pickle files
    # ____________________Spatial__________________________
N_members = 50

stats = {}

n_members = np.arange(1, N_members + 1)
counts = np.zeros(N_members)
median_time = np.zeros(N_members)
mean_time = np.zeros(N_members)
min_time = np.zeros(N_members)
std_time = np.zeros(N_members)

mean_depth = np.zeros(N_members)
median_depth = np.zeros(N_members)
std_depth = np.zeros(N_members)


for member in range(1, N_members+1):
    
    if Latitude_limit is not None:
        pkl_path = path + f"analysis/connectivity/Kh_{K_h:01d}_{Latitude_limit}N/Distributions_Kh_{K_h:01d}_m{member:03d}.pkl"
    elif Longitude_limit is not None:    
        pkl_path = path + f"analysis/connectivity/Kh_{K_h:01d}_{abs(Longitude_limit)}W/Distributions_Kh_{K_h:01d}_m{member:03d}.pkl"
        
    
    if os.path.exists(pkl_path):
        with open(pkl_path, "rb") as f:
            distributions = pickle.load(f)
        
        drift_time = distributions["drift_time"]
        depths = distributions["depths"]
        trajectory = distributions["trajectory"]
        
        median_time[member - 1] = np.median(drift_time)
        mean_time[member - 1] = np.mean(drift_time)
        min_time[member - 1] = np.min(drift_time)
        std_time[member - 1] = np.std(drift_time)
        counts[member - 1] = len(drift_time)
                    
        mean_depth[member - 1] = np.mean(depths)
        median_depth[member - 1] = np.median(depths)
        std_depth[member - 1] = np.std(depths)
    else:
        print(f"File {pkl_path} does not exist. Skipping member {member}.")
        
        median_time[member - 1] = np.nan
        mean_time[member - 1] = np.nan
        min_time[member - 1] = np.nan
        std_time[member - 1] = np.nan
        counts[member - 1] = 0
                    
        mean_depth[member - 1] = np.nan
        median_depth[member - 1] = np.nan
        std_depth[member - 1] = np.nan

stats["subset"] = n_members
stats["counts"] = counts
stats["median_time"] = median_time
stats["mean_time"] = mean_time
stats["min_time"] = min_time
stats["std_time"] = std_time
stats["mean_depth"] = mean_depth
stats["median_depth"] = median_depth
stats["std_depth"] = std_depth

stats_df = pd.DataFrame(stats)

if Latitude_limit is not None:
    save_csv_path = path + f"analysis/connectivity/Stats/Stats_Kh_{K_h:01d}_{Latitude_limit}N.csv"
elif Longitude_limit is not None:    
    save_csv_path = path + f"analysis/connectivity/Stats/Stats_Kh_{K_h:01d}_{abs(Longitude_limit)}W.csv"

# Create the directory if it does not exist
os.makedirs(os.path.dirname(save_csv_path), exist_ok=True)

stats_df.to_csv(save_csv_path)
print(f"Saved {save_csv_path}")
    
# %%
