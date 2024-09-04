# %% Load the packages
import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
import cartopy
from tqdm import tqdm
import pandas as pd

import sys

sys.path.append("../functions")
import hexbin_functions as hexfunc

# %% Load the data

location = "Cape_Hatteras"
depth = 1
stats = {}
total_members = 50
Latitude_limit = 53

n_members = np.arange(1, total_members + 1)
counts = np.zeros(total_members)
median_time = np.zeros(total_members)
mean_time = np.zeros(total_members)
std_time = np.zeros(total_members)

median_depth = np.zeros(total_members)
std_depth = np.zeros(total_members)

for member in tqdm(range(1, total_members + 1)):

    file_path = f"/storage/shared/oceanparcels/output_data/data_Claudio/NEMO_Ensemble/{location}/temporal_connectivity/dep_{depth:01d}/{location}_dep{depth:01d}_m{member:03d}.zarr"
    pset = xr.open_zarr(file_path)
    N_particles = len(pset.trajectory)

    lats = pset.lat.load().values
    p_index, t_index = np.where(lats[:, :] > Latitude_limit)
    subpolar_traj = np.unique(p_index)
    drift_time = []

    for i in subpolar_traj:
        idx_t = np.where(p_index == i)[0][0]
        drift_time.append(t_index[idx_t])
    
    drift_time = np.array(drift_time)
    
    depths = pset.z.load().values
    depths = depths[subpolar_traj, drift_time]
    
    median_time[member - 1] = np.median(drift_time)
    mean_time[member - 1] = np.mean(drift_time)
    std_time[member - 1] = np.std(drift_time)
    counts[member - 1] = len(np.unique(p_index)) / N_particles
    
    median_depth[member - 1] = np.median(depths)
    std_depth[member - 1] = np.std(depths)

stats["members"] = n_members
stats["percentage"] = counts
stats["median_time"] = median_time
stats["mean_time"] = mean_time
stats["std_time"] = std_time
stats["median_depth"] = median_depth
stats["std_depth"] = std_depth


#%%
stats_df = pd.DataFrame(stats)

fig, ax = plt.subplots(3, 1, figsize=(10, 15))
ax[0].plot(stats_df["members"], stats_df["percentage"], '.k--')
ax[0].set_title("Percentage of Subpolar Trajectories")
ax[0].set_ylabel("Percentage")
ax[0].set_xlabel("Members")
ax[0].grid()

ax[1].plot(stats_df["members"], stats_df["median_time"]/365, 'sr--')
ax[1].set_title("Median time to reach $53^o$N")
ax[1].set_xlabel("Members")
ax[1].set_ylabel("Median drift time (years)")
ax[1].grid()

ax[2].plot(stats_df["members"], stats_df["median_depth"], 'og--')
ax[2].set_title("Median depth at $53^o$N")
ax[2].set_xlabel("Time (days)")
ax[2].set_ylabel("Median depth (m)")
ax[2].grid()
#%% Can you give me the ensembple mean and std of the drift time, percentage and depth at 53N? Just by printing the values
print(f"Ensemble mean percentage of subpolar trajectories: {np.mean(stats_df['percentage'])*100:.2f}% +/- {np.std(stats_df['percentage'])*100:.2f} %")
print(f"Ensemble mean median time to reach 53N: {np.mean(stats_df['median_time'])/365:.2f} years +/- {np.std(stats_df['median_time'])/365:.2f} years")
print(f"Ensemble mean median depth at 53N: {np.mean(stats_df['median_depth']):.2f} m +/- {np.std(stats_df['median_depth']):.2f} m")


# %%%%%%%%%%%%%%%%%%%%%%%%%%%STATS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
member = 2
file_path = f"/storage/shared/oceanparcels/output_data/data_Claudio/NEMO_Ensemble/{location}/temporal_connectivity/dep_{depth:01d}/{location}_dep{depth:01d}_m{member:03d}.zarr"
pset = xr.open_zarr(file_path)
N_particles = len(pset.trajectory)

lats = pset.lat.load().values
p_index, t_index = np.where(lats[:, :] > Latitude_limit)

subpolar_traj = np.unique(p_index)

drift_time = []

for i in subpolar_traj:
    idx_t = np.where(p_index == i)[0][0]
    drift_time.append(t_index[idx_t])

drift_time = np.array(drift_time)
pset = pset.isel(trajectory=subpolar_traj)
pset = pset.compute()
n_traj = len(pset.trajectory)

#%%index pset at the drift_time of the particles
pset_at_limit = pset.isel(obs=drift_time)

# %%%%%%%%%%%%%%%%%%%%%% STATS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fig = plt.figure()
ax = plt.axes(projection=cartopy.crs.PlateCarree())
# ax.set_extent([-80, -10, 25, 55], crs=cartopy.crs.PlateCarree())
# ax.set_extent([-80, -60, 28, 37], crs=cartopy.crs.PlateCarree())
ax.coastlines()
ax.gridlines(draw_labels=True, zorder=0, linestyle="--", linewidth=0.5)

for i in range(n_traj):
    # ax.scatter(pset.lon[i,0], pset.lat[i,0], s=1)
    ax.plot(pset.lon[i, :], pset.lat[i, :])


# %%%%%%%%%%%%%%%%%% MAPS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
###################################################################
path = f"/storage/shared/oceanparcels/output_data/data_Claudio/NEMO_Ensemble/analysis/prob_distribution/Connectivity/dep_{depth:01d}/Pconn_dep_{depth:01d}_m{member:03d}.nc"

location = "Cape_Hatteras"
member = 1
depth = 1
member_list = np.arange(1, 51)

# t = 730-12*7

file_path = f"/storage/shared/oceanparcels/output_data/data_Claudio/NEMO_Ensemble/analysis/prob_distribution/Connectivity/dep_{depth:01d}/Pconn_dep_{depth:01d}_m{member:03d}.nc"
P_m = xr.open_dataset(file_path)
P_m = P_m.sortby("hexint")

hex_grid = hexfunc.int_to_hex(P_m.hexint.values)
hexbin_grid = hexfunc.hexGrid(hex_grid, h3_res=3)

prior = 1/len(member_list) # 1/50 same for all members
posterior = np.zeros((len(member_list), len(P_m.hexint)))
counts = np.zeros((len(member_list), len(P_m.hexint)))

for i, member in enumerate(member_list):
    file_path = f"/storage/shared/oceanparcels/output_data/data_Claudio/NEMO_Ensemble/analysis/prob_distribution/Connectivity/dep_{depth:01d}/Pconn_dep_{depth:01d}_m{member:03d}.nc"
    P_m = xr.open_dataset(file_path)
    P_m = P_m.sortby("hexint")
    likelihood = P_m["probability"][:, :].values
    hypothesis = likelihood
    counts[i, :] = np.nansum(hypothesis, axis=1)

ensemble_mean = np.nanmean(counts, axis=0)
ensemble_std = np.nanstd(counts, axis=0)

# %%%%%%%%%%%%%%%%%%% Ensemble Mean %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fig = plt.figure(figsize=(10, 5))
ax = plt.axes(projection=cartopy.crs.PlateCarree())
ax.set_extent([-100, 30, 5, 75], crs=cartopy.crs.PlateCarree())

ax.add_feature(cartopy.feature.LAND, zorder=0, color="black")
gl = ax.gridlines(draw_labels=True, zorder=0, linestyle="--", linewidth=0.5, alpha=0.5)
gl.top_labels = False
gl.right_labels = False

im = hexbin_grid.pcolorhex(ensemble_mean, ax=ax, cmap="magma", draw_edges=False)
ax.plot([-30, -10], [Latitude_limit, Latitude_limit], color="red")
ax.text(-18, Latitude_limit + 1, "$53^o$N", color="red")
ax.scatter(-73.6, 35.6, color="deepskyblue", edgecolor='black', s=40, marker="o", label="Release Location")
ax.set_title(f"Release at ${depth}$ m")
# ax.legend()

# Adjust the layout to make space for the color bar
plt.subplots_adjust(bottom=0.2)


cbar_ax = fig.add_axes([0.2, 0.1, 0.6, 0.05])
cbar = plt.colorbar(im, ax=ax, cax=cbar_ax, orientation="horizontal", label=f"Ensemble Average Occurrences per Bin")
plt.savefig(f"../figs/Fig5-Ensemble_mean_Connectivity-depth{depth:01d}.png", dpi=300)

# %%%%%%%%%%%%%% Ensemble STD %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fig = plt.figure(figsize=(10, 5))
ax = plt.axes(projection=cartopy.crs.PlateCarree())
ax.set_extent([-100, 30, 5, 75], crs=cartopy.crs.PlateCarree())

ax.add_feature(cartopy.feature.LAND, zorder=0, color="black")
gl = ax.gridlines(draw_labels=True, zorder=0, linestyle="--", linewidth=0.5, alpha=0.5)
gl.top_labels = False
gl.right_labels = False

im = hexbin_grid.pcolorhex(ensemble_std, ax=ax, cmap="viridis", draw_edges=False)
ax.plot([-30, -10], [Latitude_limit, Latitude_limit], color="red")
ax.text(-18, Latitude_limit + 1, "$53^o$N", color="red")
ax.scatter(-73.6, 35.6, color="deepskyblue", edgecolor='black', s=40, marker="o")

ax.set_title(f"Ensemble Mean")

ax.set_title(f"Release at ${depth}$ m")

# Adjust the layout to make space for the color bar
plt.subplots_adjust(bottom=0.2)

cbar_ax = fig.add_axes([0.2, 0.1, 0.6, 0.05])
cbar = plt.colorbar(im, ax=ax, cax=cbar_ax, orientation="horizontal", label=f"Ensemble STD Occurrences per Bin")
plt.savefig(f"../figs/Fig6-Ensemble_std_Connectivity-depth{depth:01d}.png", dpi=300)
# %%
