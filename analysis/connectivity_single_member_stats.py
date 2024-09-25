# %% Load the packages
import numpy as np
import xarray as xr
from tqdm import tqdm
import pandas as pd
import pickle
import os

location = "Cape_Hatteras"
path = f"/storage/shared/oceanparcels/output_data/data_Claudio/NEMO_Ensemble/"

Latitude_limit = None
Longitude_limit = -40

# %% Spatial analysis
distributions = {}

total_members = 50

for delta_r in [1, 2]:
    for member in tqdm(range(1, total_members + 1)):
        print(f"Member: {member:03d},  Delta_r: {delta_r}")
        file_path = path + f"{location}/spatial_long/dr_{delta_r*100:03.0f}/{location}_dr{delta_r*100:03.0f}_m{member:03d}.zarr"
        
        pset = xr.open_zarr(file_path)
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
                save_path = path + f"analysis/connectivity/dr_{delta_r*100:03.0f}_{Latitude_limit}N/Distributions_dr{delta_r*100:03.0f}_m{member:03d}.pkl"
            elif Longitude_limit is not None:    
                save_path = path + f"analysis/connectivity/dr_{delta_r*100:03.0f}_{abs(Longitude_limit)}W/Distributions_dr{delta_r*100:03.0f}_m{member:03d}.pkl"
                
            # save_path = f"/storage/shared/oceanparcels/output_data/data_Claudio/NEMO_Ensemble/analysis/connectivity/dr_{delta_r*100:03.0f}/Distributions_dr{delta_r*100:03.0f}_m{member:03d}.pkl"
            with open(save_path, "wb") as f:
                pickle.dump(distributions, f)
            
        else:
            print(f"--EMPTY--")

# %% Temporal analysis

distributions = {}
total_members = 50

for week in [4, 12, 20]:
    for member in tqdm(range(1, total_members + 1)):
        print(f"Member: {member:03d},  Week: {week}")
        file_path = path + f"{location}/temporal_long/W_{week:01d}/{location}_W{week:01d}_m{member:03d}.zarr"
        
        pset = xr.open_zarr(file_path)
        N_particles = len(pset.trajectory)

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
            
            depths = pset.z.load().values
            depths = depths[subpolar_traj, drift_time]

            distributions["member"] = member
            distributions["drift_time"] = drift_time
            distributions["depths"] = depths
            distributions["trajectory"] = np.unique(p_index)
            
            # SAVE DISTRIBUTIONS in a pickle file
            if Latitude_limit is not None:
                save_path = path + f"analysis/connectivity/W_{week:02d}_{Latitude_limit}N/Distributions_W_{week:02d}_m{member:03d}.pkl"
            elif Longitude_limit is not None:    
                save_path = path + f"analysis/connectivity/W_{week:02d}_{abs(Longitude_limit)}W/Distributions_W_{week:02d}_m{member:03d}.pkl"
            
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

for delta_r in [0.1, 1., 2.]:
    for member in range(1, N_members+1):
        
        if Latitude_limit is not None:
            pkl_path = path + f"analysis/connectivity/dr_{delta_r*100:03.0f}_{Latitude_limit}N/Distributions_dr{delta_r*100:03.0f}_m{member:03d}.pkl"
        elif Longitude_limit is not None:    
            pkl_path = path + f"analysis/connectivity/dr_{delta_r*100:03.0f}_{abs(Longitude_limit)}W/Distributions_dr{delta_r*100:03.0f}_m{member:03d}.pkl"
            
        
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
            save_csv_path = path + f"analysis/connectivity/Stats/Stats_dr{delta_r*100:03.0f}_{Latitude_limit}N.csv"
        elif Longitude_limit is not None:    
            save_csv_path = path + f"analysis/connectivity/Stats/Stats_dr{delta_r*100:03.0f}_{abs(Longitude_limit)}W.csv"
    
        stats_df.to_csv(save_csv_path)
    
#%% ___________________Temporal__________________________
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

for week in [4, 12, 20]:
    for member in range(1, N_members + 1):
        
        if Latitude_limit is not None:
            pkl_path = path + f"analysis/connectivity/W_{week:02d}_{Latitude_limit}N/Distributions_W{week:02d}_m{member:03d}.pkl"
        elif Longitude_limit is not None:    
            pkl_path = path + f"analysis/connectivity/W_{week:02d}_{abs(Longitude_limit)}W/Distributions_W{week:02d}_m{member:03d}.pkl"
            
        
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
            save_csv_path = path + f"analysis/connectivity/Stats/Stats_W{week:02d}_{Latitude_limit}N.csv"
        elif Longitude_limit is not None:    
            save_csv_path = path + f"analysis/connectivity/Stats/Stats_W{week:02d}_{abs(Longitude_limit)}W.csv"
        
        stats_df.to_csv(save_csv_path)
    
# %%
