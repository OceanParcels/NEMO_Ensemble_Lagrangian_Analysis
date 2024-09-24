# %% Load the packages
import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
import cartopy
from tqdm import tqdm
import pandas as pd
import pickle
import sys
import seaborn as sns

sys.path.append("../functions")
import hexbin_functions as hexfunc

# %% Load all dataframes

all_mix_temp = {}
all_mix_space = {}
all_temp = {}
all_space = {}
location = "Cape_Hatteras"

for week in [4, 12, 20]:
    
    df_path = f"/storage/shared/oceanparcels/output_data/data_Claudio/NEMO_Ensemble/analysis/connectivity/Stats/Stats_W{week:02d}.csv"
    _df = pd.read_csv(df_path)
    all_temp[week] = _df

    df_mix = f"/storage/shared/oceanparcels/output_data/data_Claudio/NEMO_Ensemble/analysis/connectivity/Stats/Stats_mix_W{week:02d}.csv"
    _df = pd.read_csv(df_mix)
    all_mix_temp[week] = _df
    
for delta_r in [0.1, 1., 2.]:
    df_path = f"/storage/shared/oceanparcels/output_data/data_Claudio/NEMO_Ensemble/analysis/connectivity/Stats/Stats_mix_dr{delta_r*100:03.0f}.csv"
    _df = pd.read_csv(df_path)
    all_mix_space[delta_r] = _df
    
    df_path = f"/storage/shared/oceanparcels/output_data/data_Claudio/NEMO_Ensemble/analysis/connectivity/Stats/Stats_dr{delta_r*100:03.0f}.csv"
    _df = pd.read_csv(df_path)
    all_space[delta_r] = _df


kde_space = {}
kde_temp = {}
kde_mix_space = {}
kde_mix_temp = {}

#%% Plot the percentage of subpolar trajectories for temporal and spatial members
fig, ax = plt.subplots(2, 3, figsize=(10, 6))

ax = ax.flatten()
# Plot distributions for percentage of subpolar trajectories

colors_space = ["midnightblue", "blueviolet", "teal"]
ls_space = [(0, (1, 1)), '--', '-.']
j = 0
for delta_r in [0.1, 1., 2.]:
    _counts = sns.kdeplot(all_space[delta_r]["counts"], ax=ax[0], label=f"$\delta_r = {delta_r}^o$", clip=(0, 500),
                fill=False, color=colors_space[j], linestyle=ls_space[j])
    kde_space[delta_r] = _counts.get_lines()[0].get_data() # save the KDE values for entropy computation
    
    sns.kdeplot(all_space[delta_r]["median_time"]/365, ax=ax[1], label=f"Spatial dr{delta_r}", clip=(0, 6),
                fill=False, color=colors_space[j], linestyle=ls_space[j])
    sns.kdeplot(all_space[delta_r]["min_time"]/365, ax=ax[2], label=f"Spatial dr{delta_r}", clip=(0, 6),
                fill=False, color=colors_space[j], linestyle=ls_space[j])
    j += 1

colors_temp = ["darkred", "orangered", "orange"]
ls_time = [(0, (1, 1)), '--', '-.', (0, (3, 1, 1, 1, 1, 1))]
j = 0
for week in [4, 12, 20]:
    _counts = sns.kdeplot(all_temp[week]["counts"], ax=ax[3], label=f"{week} weeks", clip=(0, 200),
                fill=False, color=colors_temp[j], linestyle=ls_time[j])
    kde_temp[week] = _counts.get_lines()[0].get_data() # save the KDE values for entropy computation
    
    sns.kdeplot(all_temp[week]["median_time"]/365, ax=ax[4], label=f"Temporal W{week}", clip=(0, 6),
                fill=False, color=colors_temp[j], linestyle=ls_time[j])
    sns.kdeplot(all_temp[week]["min_time"]/365, ax=ax[5], label=f"Temporal W{week}", clip=(0, 6),
                fill=False, color=colors_temp[j], linestyle=ls_time[j])
    j += 1

# Add labels 'A', 'B', ... for each subplot in the top left corner
labels = ['A', 'B', 'C', 'D', 'E', 'F']
for i, label in enumerate(labels):
    ax[i].text(0.05, 0.95, label, transform=ax[i].transAxes, fontsize=12, fontweight='bold', va='top', ha='left')

ax[0].set_xlabel("Counts")
ax[0].set_ylabel("Density")
ax[0].legend()

ax[1].set_xlabel("Median Drift Time (years)")
ax[1].set_ylabel("Density")

ax[2].set_xlabel("Minimum Drift Time (years)")
ax[2].set_ylabel("Density")

ax[3].set_xlabel("Counts")
ax[3].set_ylabel("Density")
ax[3].legend()

ax[4].set_xlabel("Median Drift Time (years)")
ax[4].set_ylabel("Density")

ax[5].set_xlabel("Minimum Drift Time (years)")
ax[5].set_ylabel("Density")

plt.tight_layout()
# save the figure
plt.savefig("../figs/Figx-Connect_tempNspace.png", dpi=300)

#%% Plot the percentage of subpolar trajectories for Mixture temporal and spatial members
fig, ax = plt.subplots(2, 3, figsize=(10, 6))

ax = ax.flatten()
# Plot distributions for percentage of subpolar trajectories

colors_space = ["midnightblue", "blueviolet", "teal"]
ls_space = [(0, (1, 1)), '--', '-.']
j=0
for delta_r in [0.1, 1., 2.]:
    _counts = sns.kdeplot(all_mix_space[delta_r]["counts"], ax=ax[0], label=f"Mix. $\delta_r = {delta_r}^o$", clip=(0, 100),
                fill=False, color=colors_space[j], linestyle=ls_space[j])
    kde_mix_space[delta_r] = _counts.get_lines()[0].get_data() # save the KDE values for entropy computation
    
    sns.kdeplot(all_mix_space[delta_r]["median_time"]/365, ax=ax[1], label=f"Mix. $\delta_r = {delta_r}^o$", clip=(0, 6),
                fill=False, color=colors_space[j], linestyle=ls_space[j])
    # sns.kdeplot(all_mix_space[delta_r]["std_time"]/365, ax=ax[5], label=f"Mix. $\delta_r = {delta_r}^o$", clip=(0, 6),
    #             fill=False, color=colors_space[j], linestyle=ls_space[j])
    sns.kdeplot(all_mix_space[delta_r]["min_time"]/365, ax=ax[2], label=f"Mix. $\delta_r = {delta_r}^o$", clip=(0, 6),
                fill=False, color=colors_space[j], linestyle=ls_space[j])
    j += 1

colors_temp = ["darkred", "orangered", "orange"]
ls_time = [(0, (1, 1)), '--', '-.', (0, (3, 1, 1, 1, 1, 1))]
j = 0
for week in [4, 12, 20]:
    _counts = sns.kdeplot(all_mix_temp[week]["counts"], ax=ax[3], label=f"Mix. {week} weeks", clip=(0, 100),
                fill=False, color=colors_temp[j], linestyle=ls_time[j])
    kde_mix_temp[week] = _counts.get_lines()[0].get_data() # save the KDE values for entropy computation
    
    sns.kdeplot(all_mix_temp[week]["median_time"]/365, ax=ax[4], label=f"Mix. {week} weeks", clip=(0, 6),
                fill=False, color=colors_temp[j], linestyle=ls_time[j])
    # sns.kdeplot(all_mix_temp[week]["std_time"]/365, ax=ax[2], label=f"Mix. {week} weeks", clip=(0, 6),
    #             fill=False, color=colors_temp[j], linestyle=ls_time[j])
    sns.kdeplot(all_mix_temp[week]["min_time"]/365, ax=ax[5], label=f"Mix. {week} weeks", clip=(0, 6),
               fill=False, color=colors_temp[j], linestyle=ls_time[j])
    j += 1

# Add labels 'A', 'B', ... for each subplot in the top left corner
labels = ['A', 'B', 'C', 'D', 'E', 'F']
for i, label in enumerate(labels):
    ax[i].text(0.05, 0.95, label, transform=ax[i].transAxes, fontsize=12, fontweight='bold', va='top', ha='left')
    

ax[0].set_xlabel("Counts")
ax[0].set_ylabel("Density")
ax[0].legend()
# ax[0].grid()

ax[1].set_xlabel("Median Drift Time (years)")
ax[1].set_ylabel("Density")
# ax[1].grid()

# ax[2].set_xlabel("STD Drift Time (years)")
ax[2].set_xlabel("Minimum Drift Time (years)")
ax[2].set_ylabel("Density")
# ax[2].grid()

ax[3].set_xlabel("Counts")
ax[3].set_ylabel("Density")
ax[3].legend(fontsize=7)
# ax[3].grid()

ax[4].set_xlabel("Median Drift Time (years)")
ax[4].set_ylabel("Density")
# ax[4].grid()

# ax[5].set_xlabel("STD Drift Time (years)")
ax[5].set_xlabel("Minimum Drift Time (years)")
ax[5].set_ylabel("Density")
# ax[5].grid()

plt.tight_layout()
plt.savefig("../figs/Figx-Connect_MIX_tempNspace.png", dpi=300)

###################################################################
#%% COMPARISON OF TEMPORAL AND SPATIAL CONNECTIVITY with Mixture
#%% Plot the percentage of subpolar trajectories for temporal and spatial members
fig, ax = plt.subplots(2, 3, figsize=(10, 6))

ax = ax.flatten()
# Plot distributions for percentage of subpolar trajectories

colors_space = ["midnightblue", "blueviolet", "teal"]
ls_space = [(0, (1, 1)), '--', '-.']
j = 0
for delta_r in [0.1, 1., 2.]:
    sns.kdeplot(all_space[delta_r]["counts"], ax=ax[0], label=f"$\delta_r = {delta_r}^o$", clip=(0, 500),
                fill=False, color=colors_space[j], linestyle=ls_space[j])
    sns.kdeplot(all_space[delta_r]["median_time"]/365, ax=ax[1], label=f"Spatial dr{delta_r}", clip=(0, 6),
                fill=False, color=colors_space[j], linestyle=ls_space[j])
    sns.kdeplot(all_space[delta_r]["min_time"]/365, ax=ax[2], label=f"Spatial dr{delta_r}", clip=(0, 6),
                fill=False, color=colors_space[j], linestyle=ls_space[j])
    j += 1

dr_ref = 0.1
sns.kdeplot(all_mix_space[dr_ref]["counts"], ax=ax[0], label=f"Mix. $\delta_r = {dr_ref}^o$", clip=(0, 100),
            fill=False, color='k', linestyle='-')
sns.kdeplot(all_mix_space[dr_ref]["median_time"]/365, ax=ax[1], label=f"Mix. $\delta_r = 1^o$", clip=(0, 6),
            fill=False, color='k', linestyle='-')
sns.kdeplot(all_mix_space[dr_ref]["min_time"]/365, ax=ax[2], label=f"Mix. $\delta_r = 1^o$", clip=(0, 6),
            fill=False, color='k', linestyle='-')

# sns.kdeplot(all_mix_temp[4]["counts"], ax=ax[0], label=f"Mix. 4 weeks", clip=(0, 100),
#             fill=False, color='k', linestyle=(0, (3, 1, 1, 3)))
# sns.kdeplot(all_mix_temp[4]["median_time"]/365, ax=ax[1], label=f"Mix. 4 weeks", clip=(0, 6),
#             fill=False, color='k', linestyle=(0, (3, 1, 1, 3)))
# sns.kdeplot(all_mix_temp[4]["min_time"]/365, ax=ax[2], label=f"Mix. 4 weeks", clip=(0, 6),
#             fill=False, color='k', linestyle=(0, (3, 1, 1, 3)))

colors_temp = ["darkred", "orangered", "orange"]
ls_time = [(0, (1, 1)), '--', '-.', (0, (3, 1, 1, 1, 1, 1))]
j = 0
for week in [4, 12, 20]:
    sns.kdeplot(all_temp[week]["counts"], ax=ax[3], label=f"{week} weeks", clip=(0, 200),
                fill=False, color=colors_temp[j], linestyle=ls_time[j])
    sns.kdeplot(all_temp[week]["median_time"]/365, ax=ax[4], label=f"Temporal W{week}", clip=(0, 6),
                fill=False, color=colors_temp[j], linestyle=ls_time[j])
    sns.kdeplot(all_temp[week]["min_time"]/365, ax=ax[5], label=f"Temporal W{week}", clip=(0, 6),
                fill=False, color=colors_temp[j], linestyle=ls_time[j])
    j += 1
    
dr_ref = 0.1
sns.kdeplot(all_mix_space[dr_ref]["counts"], ax=ax[3], label=f"Mix. $\delta_r = {dr_ref}^o$", clip=(0, 100),
            fill=False, color='k', linestyle='-')
sns.kdeplot(all_mix_space[dr_ref]["median_time"]/365, ax=ax[4], label=f"Mix. $\delta_r = 1^o$", clip=(0, 6),
            fill=False, color='k', linestyle='-')
sns.kdeplot(all_mix_space[dr_ref]["min_time"]/365, ax=ax[5], label=f"Mix. $\delta_r = 1^o$", clip=(0, 6),
            fill=False, color='k', linestyle='-')

sns.kdeplot(all_mix_temp[4]["counts"], ax=ax[3], label=f"Mix. 4 weeks", clip=(0, 100),
            fill=False, color='k', linestyle=(0, (3, 1, 1, 3)))
sns.kdeplot(all_mix_temp[4]["median_time"]/365, ax=ax[4], label=f"Mix. 4 weeks", clip=(0, 6),
            fill=False, color='k', linestyle=(0, (3, 1, 1, 3)))
sns.kdeplot(all_mix_temp[4]["min_time"]/365, ax=ax[5], label=f"Mix. 4 weeks", clip=(0, 6),
            fill=False, color='k', linestyle=(0, (3, 1, 1, 3)))

# Add labels 'A', 'B', ... for each subplot in the top left corner
labels = ['A', 'B', 'C', 'D', 'E', 'F']
for i, label in enumerate(labels):
    ax[i].text(0.05, 0.95, label, transform=ax[i].transAxes, fontsize=12, fontweight='bold', va='top', ha='left')

ax[0].set_xlabel("Counts")
ax[0].set_ylabel("Density")
ax[0].legend()

ax[1].set_xlabel("Median Drift Time (years)")
ax[1].set_ylabel("Density")

ax[2].set_xlabel("Minimum Drift Time (years)")
ax[2].set_ylabel("Density")

ax[3].set_xlabel("Counts")
ax[3].set_ylabel("Density")
ax[3].legend(fontsize=8)

ax[4].set_xlabel("Median Drift Time (years)")
ax[4].set_ylabel("Density")

ax[5].set_xlabel("Minimum Drift Time (years)")
ax[5].set_ylabel("Density")

plt.tight_layout()
# save the figure
plt.savefig("../figs/Figx-Connect_Comparisson.png", dpi=300)

#%% Average the temporal and spatial members and mixtures and put the values in a dataframe

zero_counts = np.zeros(3*4)
types = []

average_mean_time = np.zeros(3*4)
average_min_time = np.zeros(3*4)
average_mean_depth = np.zeros(3*4)

i = 0
for week in [4, 12, 20]:
    zero_counts[i] = (all_temp[week]['counts'] == 0).sum()
    types.append(f"{week} weeks")
    average_mean_time[i] = all_temp[week]["mean_time"].mean()
    average_min_time[i] = all_temp[week]["min_time"].mean()
    average_mean_depth[i] = all_temp[week]["mean_depth"].mean()
    i+=1

for week in [4, 12, 20]:
    zero_counts[i] = (all_mix_temp[week]['counts'] == 0).sum()
    types.append(f"Mixture {week} weeks")
    average_mean_time[i] = all_mix_temp[week]["mean_time"].mean()
    average_min_time[i] = all_mix_temp[week]["min_time"].mean()
    average_mean_depth[i] = all_mix_temp[week]["mean_depth"].mean()
    i+=1

for delta_r in [0.1, 1., 2.]:
    zero_counts[i] = (all_space[delta_r]['counts'] == 0).sum()
    types.append(f"dr{delta_r}")
    average_mean_time[i] = all_space[delta_r]["mean_time"].mean()
    average_min_time[i] = all_space[delta_r]["min_time"].mean()
    average_mean_depth[i] = all_space[delta_r]["mean_depth"].mean()
    i+=1
    
for delta_r in [0.1, 1., 2.]:
    zero_counts[i] = (all_mix_space[delta_r]['counts'] == 0).sum()
    types.append(f"Mixture dr{delta_r}")
    average_mean_time[i] = all_mix_space[delta_r]["mean_time"].mean()
    average_min_time[i] = all_mix_space[delta_r]["min_time"].mean()
    average_mean_depth[i] = all_mix_space[delta_r]["mean_depth"].mean()
    i+=1
    
# Create the dataframe
average_df = pd.DataFrame({"Type": types, "Zero_counts": zero_counts, "AVG_Mean_time": average_mean_time,
                            "AVG_Min_time": average_min_time, "AVG_Mean_depth": average_mean_depth})


#%% compute the KDE of the counts for the temporal and spatial members
# then compute the entropy of the KDE. No plotting here
    
# Compute the entropy of the KDE
from scipy.stats import entropy

types = ["Spatial 0.1", "Spatial 1.0", "Spatial 2.0", "Temporal 4", "Temporal 12", "Temporal 20"]

entropy_single = np.zeros(6)
entropy_mixture = np.zeros(6)

j = 0
for delta_r in [0.1, 1., 2.]:
    entropy_single[j] = entropy(kde_space[delta_r][0])
    entropy_mixture[j] = entropy(kde_mix_space[delta_r][1])
    j += 1
    
for j, week in enumerate([4, 12, 20]):
    entropy_single[j+3] = entropy(kde_temp[week][1])
    entropy_mixture[j+3] = entropy(kde_mix_temp[week][1])
    j += 1

# # Create the dataframe
# entropy_df = pd.DataFrame({"Type": types, "Entropy": [entropy_space[0.1], entropy_space[1.], entropy_space[2.],
#                                                             entropy_temp[4], entropy_temp[12], entropy_temp[20]]})

# entropy_mix_df = pd.DataFrame({"Type": types, "Entropy": [entropy_mix_space[0.1], entropy_mix_space[1.], entropy_mix_space[2.],
#                                                             entropy_mix_temp[4], entropy_mix_temp[12], entropy_mix_temp[20]]})


# %%%%%%%%%%%%%%%%%%%%%%%%%%%STATS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# member = 2
# file_path = f"/storage/shared/oceanparcels/output_data/data_Claudio/NEMO_Ensemble/{location}/temporal_connectivity/dep_{depth:01d}/{location}_dep{depth:01d}_m{member:03d}.zarr"
# pset = xr.open_zarr(file_path)
# N_particles = len(pset.trajectory)

# lats = pset.lat.load().values
# p_index, t_index = np.where(lats[:, :] > Latitude_limit)

# subpolar_traj = np.unique(p_index)

# drift_time = []

# for i in subpolar_traj:
#     idx_t = np.where(p_index == i)[0][0]
#     drift_time.append(t_index[idx_t])

# drift_time = np.array(drift_time)
# pset = pset.isel(trajectory=subpolar_traj)
# pset = pset.compute()
# n_traj = len(pset.trajectory)

# #%%index pset at the drift_time of the particles
# pset_at_limit = pset.isel(obs=drift_time)

# # %%%%%%%%%%%%%%%%%%%%%% STATS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# fig = plt.figure()
# ax = plt.axes(projection=cartopy.crs.PlateCarree())
# # ax.set_extent([-80, -10, 25, 55], crs=cartopy.crs.PlateCarree())
# # ax.set_extent([-80, -60, 28, 37], crs=cartopy.crs.PlateCarree())
# ax.coastlines()
# ax.gridlines(draw_labels=True, zorder=0, linestyle="--", linewidth=0.5)

# for i in range(n_traj):
#     # ax.scatter(pset.lon[i,0], pset.lat[i,0], s=1)
#     ax.plot(pset.lon[i, :], pset.lat[i, :])


# # %%%%%%%%%%%%%%%%%% MAPS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# ###################################################################
# path = f"/storage/shared/oceanparcels/output_data/data_Claudio/NEMO_Ensemble/analysis/prob_distribution/Connectivity/dep_{depth:01d}/Pconn_dep_{depth:01d}_m{member:03d}.nc"

# location = "Cape_Hatteras"
# member = 1
# depth = 1
# member_list = np.arange(1, 51)

# # t = 730-12*7

# file_path = f"/storage/shared/oceanparcels/output_data/data_Claudio/NEMO_Ensemble/analysis/prob_distribution/Connectivity/dep_{depth:01d}/Pconn_dep_{depth:01d}_m{member:03d}.nc"
# P_m = xr.open_dataset(file_path)
# P_m = P_m.sortby("hexint")

# hex_grid = hexfunc.int_to_hex(P_m.hexint.values)
# hexbin_grid = hexfunc.hexGrid(hex_grid, h3_res=3)

# prior = 1/len(member_list) # 1/50 same for all members
# posterior = np.zeros((len(member_list), len(P_m.hexint)))
# counts = np.zeros((len(member_list), len(P_m.hexint)))

# for i, member in enumerate(member_list):
#     file_path = f"/storage/shared/oceanparcels/output_data/data_Claudio/NEMO_Ensemble/analysis/prob_distribution/Connectivity/dep_{depth:01d}/Pconn_dep_{depth:01d}_m{member:03d}.nc"
#     P_m = xr.open_dataset(file_path)
#     P_m = P_m.sortby("hexint")
#     likelihood = P_m["probability"][:, :].values
#     hypothesis = likelihood
#     counts[i, :] = np.nansum(hypothesis, axis=1)

# ensemble_mean = np.nanmean(counts, axis=0)
# ensemble_std = np.nanstd(counts, axis=0)

# # %%%%%%%%%%%%%%%%%%% Ensemble Mean %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# fig = plt.figure(figsize=(10, 5))
# ax = plt.axes(projection=cartopy.crs.PlateCarree())
# ax.set_extent([-100, 30, 5, 75], crs=cartopy.crs.PlateCarree())

# ax.add_feature(cartopy.feature.LAND, zorder=0, color="black")
# gl = ax.gridlines(draw_labels=True, zorder=0, linestyle="--", linewidth=0.5, alpha=0.5)
# gl.top_labels = False
# gl.right_labels = False

# im = hexbin_grid.pcolorhex(ensemble_mean, ax=ax, cmap="magma", draw_edges=False)
# ax.plot([-30, -10], [Latitude_limit, Latitude_limit], color="red")
# ax.text(-18, Latitude_limit + 1, "$53^o$N", color="red")
# ax.scatter(-73.6, 35.6, color="deepskyblue", edgecolor='black', s=40, marker="o", label="Release Location")
# ax.set_title(f"Release at ${depth}$ m")
# # ax.legend()

# # Adjust the layout to make space for the color bar
# plt.subplots_adjust(bottom=0.2)


# cbar_ax = fig.add_axes([0.2, 0.1, 0.6, 0.05])
# cbar = plt.colorbar(im, ax=ax, cax=cbar_ax, orientation="horizontal", label=f"Ensemble Average Occurrences per Bin")
# plt.savefig(f"../figs/Fig5-Ensemble_mean_Connectivity-depth{depth:01d}.png", dpi=300)

# # %%%%%%%%%%%%%% Ensemble STD %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# fig = plt.figure(figsize=(10, 5))
# ax = plt.axes(projection=cartopy.crs.PlateCarree())
# ax.set_extent([-100, 30, 5, 75], crs=cartopy.crs.PlateCarree())

# ax.add_feature(cartopy.feature.LAND, zorder=0, color="black")
# gl = ax.gridlines(draw_labels=True, zorder=0, linestyle="--", linewidth=0.5, alpha=0.5)
# gl.top_labels = False
# gl.right_labels = False

# im = hexbin_grid.pcolorhex(ensemble_std, ax=ax, cmap="viridis", draw_edges=False)
# ax.plot([-30, -10], [Latitude_limit, Latitude_limit], color="red")
# ax.text(-18, Latitude_limit + 1, "$53^o$N", color="red")
# ax.scatter(-73.6, 35.6, color="deepskyblue", edgecolor='black', s=40, marker="o")

# ax.set_title(f"Ensemble Mean")

# ax.set_title(f"Release at ${depth}$ m")

# # Adjust the layout to make space for the color bar
# plt.subplots_adjust(bottom=0.2)

# cbar_ax = fig.add_axes([0.2, 0.1, 0.6, 0.05])
# cbar = plt.colorbar(im, ax=ax, cax=cbar_ax, orientation="horizontal", label=f"Ensemble STD Occurrences per Bin")
# plt.savefig(f"../figs/Fig6-Ensemble_std_Connectivity-depth{depth:01d}.png", dpi=300)
# # %%
