# %% Load the packages
import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
import cartopy
from tqdm import tqdm
import pandas as pd
import pickle
import sys

sys.path.append("../functions")
import hexbin_functions as hexfunc

# %% Spatial analysis
location = "Cape_Hatteras"
stats = {}
distributions = {}

total_members = 50
Latitude_limit = 53

n_members = np.arange(1, total_members + 1)
counts = np.zeros(total_members)
median_time = np.zeros(total_members)
mean_time = np.zeros(total_members)
std_time = np.zeros(total_members)

mean_depth = np.zeros(total_members)
median_depth = np.zeros(total_members)
std_depth = np.zeros(total_members)

for delta_r in [0.1, 1, 2]:
    for member in tqdm(range(1, total_members + 1)):
        print(f"Member: {member:03d},  Delta_r: {delta_r}")
        file_path = f"/storage/shared/oceanparcels/output_data/data_Claudio/NEMO_Ensemble/{location}/spatial_long/dr_{delta_r*100:03.0f}/{location}_dr{delta_r*100:03.0f}_m{member:03d}.zarr"
        
        pset = xr.open_zarr(file_path)
        N_particles = len(pset.trajectory)

        lats = pset.lat.load().values
        p_index, t_index = np.where(lats[:, :] > Latitude_limit)
        subpolar_traj = np.unique(p_index)
        drift_time = []

        if len(subpolar_traj) > 0:
            for i in subpolar_traj:
                idx_t = np.where(p_index == i)[0][0]
                drift_time.append(t_index[idx_t])
            
            drift_time = np.array(drift_time)
            
            depths = pset.z.load().values
            depths = depths[subpolar_traj, drift_time]
            
            median_time[member - 1] = np.median(drift_time)
            mean_time[member - 1] = np.mean(drift_time)
            std_time[member - 1] = np.std(drift_time)
            counts[member - 1] = len(np.unique(p_index)) / N_particles * 100
            
            mean_depth[member - 1] = np.mean(depths)
            median_depth[member - 1] = np.median(depths)
            std_depth[member - 1] = np.std(depths)

            distributions["member"] = member
            distributions["drift_time"] = drift_time
            distributions["depths"] = depths
            distributions["trajectory"] = np.unique(p_index)
            
            # SAVE DISTRIBUTIONS in a pickle file
            save_path = f"/storage/shared/oceanparcels/output_data/data_Claudio/NEMO_Ensemble/analysis/connectivity/dr_{delta_r*100:03.0f}/Distributions_dr{delta_r*100:03.0f}_m{member:03d}.pkl"
            with open(save_path, "wb") as f:
                pickle.dump(distributions, f)
            
        else:
            median_time[member - 1] = np.nan
            mean_time[member - 1] = np.nan
            std_time[member - 1] = np.nan
            counts[member - 1] = np.nan
            
            mean_depth[member - 1] = np.nan
            median_depth[member - 1] = np.nan
            std_depth[member - 1] = np.nan


    stats["members"] = n_members
    stats["percentage"] = counts
    stats["median_time"] = median_time
    stats["mean_time"] = mean_time
    stats["std_time"] = std_time
    stats["mean_depth"] = mean_depth
    stats["median_depth"] = median_depth
    stats["std_depth"] = std_depth

    stats_df = pd.DataFrame(stats)

    save_csv_path = f"/storage/shared/oceanparcels/output_data/data_Claudio/NEMO_Ensemble/analysis/connectivity/Stats/Stats_dr{delta_r*100:03.0f}.csv"
    stats_df.to_csv(save_csv_path)

# %% Temporal analysis
location = "Cape_Hatteras"
stats = {}
distributions = {}

total_members = 50
Latitude_limit = 53

n_members = np.arange(1, total_members + 1)
counts = np.zeros(total_members)
median_time = np.zeros(total_members)
mean_time = np.zeros(total_members)
std_time = np.zeros(total_members)

mean_depth = np.zeros(total_members)
median_depth = np.zeros(total_members)
std_depth = np.zeros(total_members)

for week in [4, 12, 20]:
    for member in tqdm(range(1, total_members + 1)):
        print(f"Member: {member:03d},  Week: {week}")

        file_path = path = f"/storage/shared/oceanparcels/output_data/data_Claudio/NEMO_Ensemble/{location}/temporal_long/W_{week:01d}/{location}_W{week:01d}_m{member:03d}.zarr"
        
        pset = xr.open_zarr(file_path)
        N_particles = len(pset.trajectory)

        lats = pset.lat.load().values
        p_index, t_index = np.where(lats[:, :] > Latitude_limit)
        subpolar_traj = np.unique(p_index)
        drift_time = []
        
        if len(subpolar_traj) > 0:
            for i in subpolar_traj:
                idx_t = np.where(p_index == i)[0][0]
                drift_time.append(t_index[idx_t])
            
            drift_time = np.array(drift_time)
            
            depths = pset.z.load().values
            
            depths = depths[subpolar_traj, drift_time]
            
            median_time[member - 1] = np.median(drift_time)
            mean_time[member - 1] = np.mean(drift_time)
            std_time[member - 1] = np.std(drift_time)
            counts[member - 1] = len(np.unique(p_index)) / N_particles * 100
            
            mean_depth[member - 1] = np.mean(depths)
            median_depth[member - 1] = np.median(depths)
            std_depth[member - 1] = np.std(depths)

            distributions["member"] = member
            distributions["drift_time"] = drift_time
            distributions["depths"] = depths
            distributions["trajectory"] = np.unique(p_index)

            # SAVE DISTRIBUTIONS in a pickle file
            save_path = f"/storage/shared/oceanparcels/output_data/data_Claudio/NEMO_Ensemble/analysis/connectivity/W_{week:02d}/Distributions_W{week:02d}_m{member:03d}.pkl"
            with open(save_path, "wb") as f:
                pickle.dump(distributions, f)
                
        else:
            median_time[member - 1] = np.nan
            mean_time[member - 1] = np.nan
            std_time[member - 1] = np.nan
            counts[member - 1] = np.nan
            
            mean_depth[member - 1] = np.nan
            median_depth[member - 1] = np.nan
            std_depth[member - 1] = np.nan
        

    stats["members"] = n_members
    stats["percentage"] = counts
    stats["median_time"] = median_time
    stats["mean_time"] = mean_time
    stats["std_time"] = std_time
    stats["mean_depth"] = mean_depth
    stats["median_depth"] = median_depth
    stats["std_depth"] = std_depth

    stats_df = pd.DataFrame(stats)

    save_csv_path = f"/storage/shared/oceanparcels/output_data/data_Claudio/NEMO_Ensemble/analysis/connectivity/Stats/Stats_W{week:02d}_m{member:03d}.pkl.csv"
    stats_df.to_csv(save_csv_path)