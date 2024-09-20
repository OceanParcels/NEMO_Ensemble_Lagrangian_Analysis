# %% Load the packages
import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
import cartopy
from tqdm import tqdm
import pandas as pd
import pickle
import sys
import concurrent.futures

sys.path.append("../functions")
import hexbin_functions as hexfunc

# %%  Temporal analysis

members = np.arange(2, 51)
N_subsets = 50

location = 'Cape_Hatteras'
subset_particles = 150

def process_member(member, week, location, subset_particles):
    path = f"/storage/shared/oceanparcels/output_data/data_Claudio/NEMO_Ensemble/{location}/temporal_long/W_{week:01d}/{location}_W{week:01d}_m{member:03d}.zarr"
    pset = xr.open_zarr(path)
    pset = pset.isel(trajectory=np.random.choice(pset.trajectory, subset_particles, replace=False))
    return pset


stats = {}
distributions = {}

Latitude_limit = 53

n_members = np.arange(1, N_subsets + 1)
counts = np.zeros(N_subsets)
median_time = np.zeros(N_subsets)
mean_time = np.zeros(N_subsets)
std_time = np.zeros(N_subsets)

mean_depth = np.zeros(N_subsets)
median_depth = np.zeros(N_subsets)
std_depth = np.zeros(N_subsets)

for week in [20]:
    for k in tqdm(range(1, N_subsets+1)):
        
        member = 1
        
        path = f"/storage/shared/oceanparcels/output_data/data_Claudio/NEMO_Ensemble/{location}/temporal_long/W_{week:01d}/{location}_W{week:01d}_m{member:03d}.zarr"
        pset_members = xr.open_zarr(path)
        obs_length = len(pset_members.obs)
        obs_range = range(obs_length) # Number of time steps in the observation period

        pset_members = pset_members.isel(trajectory=np.random.choice(pset_members.trajectory, subset_particles, replace=False))

        with concurrent.futures.ThreadPoolExecutor() as executor:
            futures = [executor.submit(process_member, member, week, location, subset_particles) for member in members]
            for future in tqdm(concurrent.futures.as_completed(futures), total=len(members)):
                pset = future.result()
            
                pset_members = xr.concat([pset_members, pset], dim='trajectory')
        
        print(f"Subset:{k} week: {week}. Number of particles: {len(pset_members.trajectory)}")

        N_particles = len(pset_members.trajectory)

        lats = pset_members.lat.load().values
        p_index, t_index = np.where(lats[:, :] > Latitude_limit)
        subpolar_traj = np.unique(p_index)
        drift_time = []

        if len(subpolar_traj) > 0:
            for i in subpolar_traj:
                idx_t = np.where(p_index == i)[0][0]
                drift_time.append(t_index[idx_t])
            
            drift_time = np.array(drift_time)
            
            depths = pset_members.z.load().values
            depths = depths[subpolar_traj, drift_time]
            
            median_time[member - 1] = np.median(drift_time)
            mean_time[member - 1] = np.mean(drift_time)
            std_time[member - 1] = np.std(drift_time)
            counts[member - 1] = len(np.unique(p_index)) / N_particles * 100
            
            mean_depth[member - 1] = np.mean(depths)
            median_depth[member - 1] = np.median(depths)
            std_depth[member - 1] = np.std(depths)

            distributions["member"] = k
            distributions["drift_time"] = drift_time
            distributions["depths"] = depths
            distributions["trajectory"] = np.unique(p_index)
            
            # SAVE DISTRIBUTIONS in a pickle file
            save_path = f"/storage/shared/oceanparcels/output_data/data_Claudio/NEMO_Ensemble/analysis/connectivity/mix_W_{week:02d}/Distributions_mix_W{week:02d}_s{k:03d}.pkl"
            with open(save_path, "wb") as f:
                pickle.dump(distributions, f)
            
        else:
            print(f"--EMPTY--")
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

    save_csv_path = f"/storage/shared/oceanparcels/output_data/data_Claudio/NEMO_Ensemble/analysis/connectivity/Stats/Stats_mix_W{week:02d}.csv"
    stats_df.to_csv(save_csv_path)


# %%
