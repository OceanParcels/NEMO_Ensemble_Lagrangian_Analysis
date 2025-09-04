#%%
import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
import cartopy
from tqdm import tqdm
# import pickle
import cmocean.cm as cmo

import sys
sys.path.append('../functions')
import deprecated_hexbin_functions as hexfunc
from mpl_toolkits.axes_grid1 import make_axes_locatable

def Shannon_entropy(Pdf):
    # Pdf = Pdf / np.nansum(Pdf)  # Normalize Pdf to sum to 1, ignoring NaNs
    # Replace zeros with a very small number to avoid log(0)
    # Pdf_safe = np.where(Pdf > 0, Pdf, np.finfo(float).eps)
    Pdf_safe = Pdf
    return -np.nansum(Pdf_safe * np.log2(Pdf_safe))

#%%
location = 'Cape_Hatteras'
week = 20
member = 1
member_list = np.arange(1, 51)

# Define the file path for the NetCDF file containing probability distributions
file_path = f"/Volumes/Claudio SSD/Ensemble_article_data/analysis/prob_distribution/{location}_temporal_long/P_W{week:01d}_m{member:03d}.nc"

# Open the dataset and sort by 'hexint'
P_m = xr.open_dataset(file_path)
P_m = P_m.sortby('hexint')

# Convert hexint values to hex grid and create a hexbin grid with resolution 3
hex_grid = hexfunc.int_to_hex(P_m.hexint.values)
hexbin_grid = hexfunc.hexGrid(hex_grid, h3_res=3)


# Dictinary with all STDs maps
std_maps_temp = {}
std_maps_space = {}


for week in [4, 12, 20]:
    print(f"Week: {week}")
    time_average = np.zeros((len(member_list), len(P_m.hexint)))


    # Loop through each member to calculate the time_average probabilities
    for i, member in tqdm(enumerate(member_list)):
        file_path = f"/Volumes/Claudio SSD/Ensemble_article_data/analysis/prob_distribution/{location}_temporal_long/P_W{week:01d}_m{member:03d}.nc"
        P_m = xr.open_dataset(file_path)
        P_m = P_m.sortby('hexint')
        likelihood = P_m['probability'][:, :].values
        hypothesis = likelihood #* prior
        time_average[i, :] = np.nanmean(hypothesis, axis=1)


    #Calculate the Standard Deviation of the time_average probabilities
    std_maps_temp[week] = np.nanstd(time_average, axis=0)

    
for delta_r in[0.1, 1., 2.]:
    print(f"Delta_r: {delta_r}")
    time_average = np.zeros((len(member_list), len(P_m.hexint)))
    
    # Loop through each member to calculate the time_average probabilities
    for i, member in tqdm(enumerate(member_list)):
        file_path = f"/Volumes/Claudio SSD/Ensemble_article_data/analysis/prob_distribution/{location}_spatial_long/P_dr{delta_r*100:03.0f}_m{member:03d}.nc"
        P_m = xr.open_dataset(file_path)
        P_m = P_m.sortby('hexint')
        likelihood = P_m['probability'][:, :].values
        hypothesis = likelihood #* prior
        time_average[i, :] = np.nanmean(hypothesis, axis=1)


    #Calculate the Standard Deviation of the time_average probabilities
    std_maps_space[delta_r] = np.nanstd(time_average, axis=0)

#%%
std_maps_diff = {}

for K_h in [10, 1000]:
    print(f"K_h: {K_h}")
    time_average = np.zeros((len(member_list), len(P_m.hexint)))
    
    # Loop through each member to calculate the time_average probabilities
    for i, member in tqdm(enumerate(member_list)):
        file_path = f"/Volumes/Claudio SSD/Ensemble_article_data/analysis/prob_distribution/{location}_diffusion_long/P_diff_Kh_{K_h:01d}_m{member:03d}.nc"
        P_m = xr.open_dataset(file_path)
        P_m = P_m.sortby('hexint')
        likelihood = P_m['probability'][:, :].values
        hypothesis = likelihood #* prior
        time_average[i, :] = np.nanmean(hypothesis, axis=1)


    #Calculate the Standard Deviation of the time_average probabilities
    std_maps_diff[K_h] = np.nanstd(time_average, axis=0)
  

#%% Maximum value std_maps_space and std_maps_time for the colorbar
std_max = np.zeros(8)
std_max[0] = np.nanmax(std_maps_temp[4])
std_max[1] = np.nanmax(std_maps_temp[12])
std_max[2] = np.nanmax(std_maps_temp[20])
std_max[3] = np.nanmax(std_maps_space[0.1])
std_max[4] = np.nanmax(std_maps_space[1.])
std_max[5] = np.nanmax(std_maps_space[2.])
std_max[6] = np.nanmax(std_maps_diff[10])
std_max[7] = np.nanmax(std_maps_diff[1000])

max_STD = np.max(std_max)

# %%
Latitude_limit = 53
ncol = 2
nrow = 4  # Changed from 4 to 5 to add a new row
fig, axs = plt.subplots(ncols=ncol, nrows=nrow, figsize=(8.27, 14.0),
                        subplot_kw={'projection': cartopy.crs.PlateCarree()},
                        sharey=True, constrained_layout=True,
                        gridspec_kw={'height_ratios': [1, 1, 1, 1]})

axs = axs.reshape(ncol*nrow)

# formatting

for i in range(0, 8):  # Now 8 plots (6 + 2 for diffusion)
    axs[i].set_extent([-100, 20, -25, 77], crs=cartopy.crs.PlateCarree())
    axs[i].add_feature(cartopy.feature.LAND, zorder=10, color='black')
    axs[i].scatter(-73.6, 35.6, color="deepskyblue", edgecolor='black', s=30, 
                   marker="o", label="Release Location", zorder=10)
    
    axs[i].plot([-40, -40], [-5, 65], color="k", zorder=10, ls='--')
    axs[i].text(-40, 35, f"${40}^o$W", color="k", zorder=10, fontsize=8, rotation= -90)
    
    gl = axs[i].gridlines(crs=cartopy.crs.PlateCarree(), draw_labels=True,
                         linewidth=0.5, color='gray', alpha=0.3)

    if i % 2 == 1:
        gl.left_labels = False

    if i < 6:
        gl.bottom_labels = False

    gl.top_labels = False
    gl.right_labels = False

    if i % 2 == 1:
        gl.right_labels = True

    # Add labels like A, B, C, etc.
    label = chr(65 + i)  # 65 is the ASCII value for 'A'
    axs[i].text(0.02, 0.02, label, transform=axs[i].transAxes, fontsize=12, fontweight='bold', va='bottom', ha='left')

# for i in range(ncol*nrow - nrow, ncol*nrow):
#     axs[i].axis('off')

# Plot the maps
i = 0
for week in [4, 12, 20]:
    hexbin_grid.pcolorhex(std_maps_temp[week], ax=axs[i], cmap=cmo.amp, draw_edges=False, 
                               maxnorm=max_STD)
    axs[i].set_title(f'{week} weeks')
    i += 2

i = 1
for delta_r in [0.1, 1., 2.]:
    im = hexbin_grid.pcolorhex(std_maps_space[delta_r], ax=axs[i], cmap=cmo.amp, draw_edges=False, 
                               maxnorm=max_STD)
    axs[i].set_title(f'$\delta_r = {delta_r}^o$')
    i += 2

# Plot diffusion maps in the 5th row (positions 6 and 7)
axs[6].set_title(r'Diffusion $K_h=10$ m$^2$/s')
hexbin_grid.pcolorhex(std_maps_diff[10], ax=axs[6], cmap=cmo.amp, draw_edges=False, maxnorm=max_STD)

axs[7].set_title(r'Diffusion $K_h=1000$ m$^2$/s')
im = hexbin_grid.pcolorhex(std_maps_diff[1000], ax=axs[7], cmap=cmo.amp, draw_edges=False, maxnorm=max_STD)

# Place colorbar on the right, vertical
colorbar_axis = fig.add_axes([1.03, 0.25, 0.05, 0.5])  # [left, bottom, width, height]
colorbar = fig.colorbar(im, cax=colorbar_axis, orientation='vertical', label='Ensemble Standard Deviation of the Time Averaged Occurrence per Bin')

# Adjust right margin to prevent colorbar from being chopped
plt.subplots_adjust(right=0.85)

# save the plot
fig.savefig(f'../figs/FigS3_Ensemble_STD.png', dpi=300, pad_inches=0.1, bbox_inches='tight')


# %% ###############################################################
# ##### Repeat analysis for Mixture distributions #################
####################################################################
location = 'Cape_Hatteras'
week = 20
subset = 1
subset_list = np.arange(1, 29)

# Define the file path for the NetCDF file containing probability distributions
file_path = f"/Volumes/Claudio SSD/Ensemble_article_data/analysis/prob_distribution/{location}_all_long/P_W{week:02d}_all_s{subset:03d}.nc"

# Open the dataset and sort by 'hexint'
P_m = xr.open_dataset(file_path)
P_m = P_m.sortby('hexint')

# Convert hexint values to hex grid and create a hexbin grid with resolution 3
hex_grid = hexfunc.int_to_hex(P_m.hexint.values)
hexbin_grid = hexfunc.hexGrid(hex_grid, h3_res=3)

# Dictinary with all STDs maps
std_maps_temp_mix = {}
std_maps_space_mix = {}

for week in [4, 12, 20]:
    print(f"Week: {week}")
    time_average = np.zeros((len(subset_list), len(P_m.hexint)))

    # Loop through each subset to calculate the time_average probabilities
    for i, subset in tqdm(enumerate(subset_list)):
        file_path = f"/Volumes/Claudio SSD/Ensemble_article_data/analysis/prob_distribution/{location}_all_long/P_W{week:02d}_all_s{subset:03d}.nc"
        P_m = xr.open_dataset(file_path)
        P_m = P_m.sortby('hexint')
        likelihood = P_m['probability'][:, :].values
        hypothesis = likelihood #* prior
        time_average[i, :] = np.nanmean(hypothesis, axis=1)

    #Calculate the Standard Deviation of the time_average probabilities
    std_maps_temp_mix[week] = np.nanstd(time_average, axis=0)
    
for delta_r in[0.1, 1., 2.]:
    print(f"Delta_r: {delta_r}")
    time_average = np.zeros((len(subset_list), len(P_m.hexint)))

    # Loop through each subset to calculate the time_average probabilities
    for i, subset in tqdm(enumerate(subset_list)):
        file_path = f"/Volumes/Claudio SSD/Ensemble_article_data/analysis/prob_distribution/{location}_all_long/P_dr{delta_r*100:03.0f}_all_s{subset:03d}.nc"
        P_m = xr.open_dataset(file_path)
        P_m = P_m.sortby('hexint')
        likelihood = P_m['probability'][:, :].values
        hypothesis = likelihood #* prior
        time_average[i, :] = np.nanmean(hypothesis, axis=1)

    #Calculate the Standard Deviation of the time_average probabilities
    std_maps_space_mix[delta_r] = np.nanstd(time_average, axis=0)


# %% PLot the maps
Latitude_limit = 44
ncol = 2
nrow = 4
fig, axs = plt.subplots(ncols=ncol, nrows=nrow, figsize=(8.27, 11.69),
                        subplot_kw={'projection': cartopy.crs.PlateCarree()},
                        sharey=True, constrained_layout=True,
                        gridspec_kw={'height_ratios': [1, 1, 1, 0.4]})

axs = axs.reshape(ncol*nrow)

# formatiiing

for i in range(0, 6):
    axs[i].set_extent([-100, 20, -25, 77], crs=cartopy.crs.PlateCarree())
    axs[i].add_feature(cartopy.feature.LAND, zorder=10, color='black')
    axs[i].scatter(-73.6, 35.6, color="deepskyblue", edgecolor='black', s=30, 
                   marker="o", label="Release Location", zorder=10)
    # axs[i].plot([-50, -10], [Latitude_limit, Latitude_limit], color="cyan", zorder=10)
    # axs[i].text(-18, Latitude_limit + 1, f"${Latitude_limit}^o$N", color="cyan", zorder=10, fontsize=8)
    
    axs[i].plot([-40, -40], [-5, 65], color="k", zorder=10, ls='--')
    axs[i].text(-40, 35, f"${40}^o$W", color="k", zorder=10, fontsize=8, rotation= -90)
    
    # axs[i].set_title(f'All Members Subset {i+1}')
    gl = axs[i].gridlines(crs=cartopy.crs.PlateCarree(), draw_labels=True,
                         linewidth=0.5, color='gray', alpha=0.3)

    if i in [1, 3, 5]:
        gl.left_labels = False

    if i < 4:
        gl.bottom_labels = False

    gl.top_labels = False
    gl.right_labels = False

    if i in [1, 3, 5]:
        gl.right_labels = True

    # Add labels like A, B, C, etc.
    label = chr(65 + i)  # 65 is the ASCII value for 'A'
    axs[i].text(0.02, 0.02, label, transform=axs[i].transAxes, fontsize=12, fontweight='bold', va='bottom', ha='left')
        

for i in range(ncol*nrow - nrow, ncol*nrow):
    axs[i].axis('off')
    
# Plot the maps
i = 1
for delta_r in [0.1, 1., 2.]:
    im = hexbin_grid.pcolorhex(std_maps_space_mix[delta_r], ax=axs[i], cmap=cmo.amp, draw_edges=False, 
                               maxnorm=max_STD)
    axs[i].set_title(f'Mixture $\delta_r = {delta_r}^o$')
    i += 2
    
i = 0
for week in [4, 12, 20]:
    hexbin_grid.pcolorhex(std_maps_temp_mix[week], ax=axs[i], cmap=cmo.amp, draw_edges=False, 
                               maxnorm=max_STD)
    axs[i].set_title(f'Mixture {week} weeks')
    i += 2

colorbar_axis = fig.add_axes([0.1, 0.06, 0.8, 0.03])  # Adjust the height of the colorbar
colorbar = fig.colorbar(im, cax=colorbar_axis, orientation='horizontal', label=f'Ensemble Standard Deviation of the Time Averaged Occurrence per Bin')

# save the plot
plt.savefig(f'../figs/FigS4_Mixture_Ensemble_STD.png', dpi=300)

#%%