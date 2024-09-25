# %% Load the packages
import numpy as np
import xarray as xr
from tqdm import tqdm
import pandas as pd
import pickle
import concurrent.futures
import os

location = "Cape_Hatteras"
base_path = f"/storage/shared/oceanparcels/output_data/data_Claudio/NEMO_Ensemble/"

Latitude_limit = 44
Longitude_limit = None

print(f"Latitude limit: {Latitude_limit}")
print(f"Longitude limit: {Longitude_limit}")
# %% Spatial analysis

members = np.arange(2, 51)
N_subsets = 50

subset_particles = 150


def process_member(member, delta_r, location, subset_particles):
    path = base_path + f"{location}/spatial_long/dr_{delta_r*100:03.0f}/{location}_dr{delta_r*100:03.0f}_m{member:03d}.zarr"
    pset = xr.open_zarr(path)
    pset = pset.isel(trajectory=np.random.choice(pset.trajectory, subset_particles, replace=False))
    
    return pset


distributions = {}

for delta_r in [0.1, 1., 2.]:
    for k in range(1, N_subsets+1):
        
        member = 1
        
        path = base_path + f"{location}/spatial_long/dr_{delta_r*100:03.0f}/{location}_dr{delta_r*100:03.0f}_m{member:03d}.zarr"
        pset_members = xr.open_zarr(path)
        obs_length = len(pset_members.obs)
        obs_range = range(obs_length) # Number of time steps in the observation period

        pset_members = pset_members.isel(trajectory=np.random.choice(pset_members.trajectory, subset_particles, replace=False))

        with concurrent.futures.ThreadPoolExecutor() as executor:
            futures = [executor.submit(process_member, member, delta_r, location, subset_particles) for member in members]
            for future in tqdm(concurrent.futures.as_completed(futures), total=len(members)):
                pset = future.result()
            
                pset_members = xr.concat([pset_members, pset], dim='trajectory')
        
        print(f"Subset:{k} delta_r: {delta_r}. Number of particles: {len(pset_members.trajectory)}")

        N_particles = len(pset_members.trajectory)

        if Latitude_limit is not None:
            lats = pset.lat.load().values
            p_index, t_index = np.where(lats[:, :] > Latitude_limit)
            
        elif Longitude_limit is not None:
            lons = pset.lon.load().values
            p_index, t_index = np.where(lons[:, :] > Longitude_limit)
        
        
        subpolar_traj = np.unique(p_index)
        drift_time = []

        if len(subpolar_traj) > 0:
            for i in subpolar_traj:
                idx_t = np.where(p_index == i)[0][0]
                drift_time.append(t_index[idx_t])
            
            drift_time = np.array(drift_time)
            
            depths = pset_members.z.load().values
            depths = depths[subpolar_traj, drift_time]

            distributions["member"] = k
            distributions["drift_time"] = drift_time
            distributions["depths"] = depths
            distributions["trajectory"] = np.unique(p_index)
            
            # SAVE DISTRIBUTIONS in a pickle file
            if Latitude_limit is not None:
                save_path = base_path + f"analysis/connectivity/mix_dr_{delta_r*100:03.0f}_{Latitude_limit}N/Distributions_mix_dr{delta_r*100:03.0f}_s{k:03d}.pkl"
            elif Longitude_limit is not None:    
                save_path = base_path + f"analysis/connectivity/mix_dr_{delta_r*100:03.0f}_{Longitude_limit}W/Distributions_mix_dr{delta_r*100:03.0f}_s{k:03d}.pkl"
            
            
            with open(save_path, "wb") as f:
                pickle.dump(distributions, f)
            
        else:
            print(f"--EMPTY--")
           
# %% Make a dataframe with the statistics
N_subsets = 50
N_particles = 7500

stats = {}

n_members = np.arange(1, N_subsets + 1)
counts = np.zeros(N_subsets)
median_time = np.zeros(N_subsets)
mean_time = np.zeros(N_subsets)
min_time = np.zeros(N_subsets)
std_time = np.zeros(N_subsets)

mean_depth = np.zeros(N_subsets)
median_depth = np.zeros(N_subsets)
std_depth = np.zeros(N_subsets)

for delta_r in [0.1, 1., 2.]:
    for k in range(1, N_subsets+1):
        
        if Latitude_limit is not None:
            pkl_path = base_path + f"analysis/connectivity/mix_dr_{delta_r*100:03.0f}_{Latitude_limit}N/Distributions_mix_dr_{delta_r*100:03.0f}_s{k:03d}.pkl"
        elif Longitude_limit is not None:    
            pkl_path = base_path + f"analysis/connectivity/mix_dr_{delta_r*100:03.0f}_{Longitude_limit}W/Distributions_mix_dr_{delta_r*100:03.0f}_s{k:03d}.pkl"
            
        
        # pkl_path = f"/storage/shared/oceanparcels/output_data/data_Claudio/NEMO_Ensemble/analysis/connectivity/mix_dr_{delta_r*100:03.0f}/Distributions_mix_dr{delta_r*100:03.0f}_s{k:03d}.pkl"
        
        
        if os.path.exists(pkl_path):
            with open(pkl_path, "rb") as f:
                distributions = pickle.load(f)
            
            drift_time = distributions["drift_time"]
            depths = distributions["depths"]
            trajectory = distributions["trajectory"]
            
            median_time[k - 1] = np.median(drift_time)
            mean_time[k - 1] = np.mean(drift_time)
            min_time[k - 1] = np.min(drift_time)
            std_time[k - 1] = np.std(drift_time)
            counts[k - 1] = len(drift_time) #/ N_particles * 100
                        
            mean_depth[k - 1] = np.mean(depths)
            median_depth[k - 1] = np.median(depths)
            std_depth[k - 1] = np.std(depths)
        else:
            print(f"File {pkl_path} does not exist. Skipping subset {k}.")
            
            median_time[k - 1] = np.nan
            mean_time[k - 1] = np.nan
            min_time[k - 1] = np.nan
            std_time[k - 1] = np.nan
            counts[k - 1] = 0
                        
            mean_depth[k - 1] = np.nan
            median_depth[k - 1] = np.nan
            std_depth[k - 1] = np.nan

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
            save_csv_path = base_path + f"analysis/connectivity/Stats/Stats_mix_dr{delta_r*100:03.0f}_{Latitude_limit}N.csv"
        elif Longitude_limit is not None:    
            save_csv_path = base_path + f"analysis/connectivity/Stats/Stats_mix_dr{delta_r*100:03.0f}_{Longitude_limit}W.csv"
        
        stats_df.to_csv(save_csv_path)


# %%
