#%%
import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
import cartopy
from tqdm import tqdm
# import pickle

import sys
sys.path.append('../functions')
import deprecated_hexbin_functions as hexfunc

#%% Load the data
# MIXTURE
location = 'Cape_Hatteras'
delta_r = 0.1
subset = 43
member = 22

# Define the file path for the NetCDF file containing probability distributions
file_path = f"/storage/shared/oceanparcels/output_data/data_Claudio/NEMO_Ensemble/analysis/prob_distribution/{location}_all_long/P_dr{delta_r*100:03.0f}_all_s{subset:03d}.nc"

# Open the dataset and sort by 'hexint'
P_mix = xr.open_dataset(file_path)
P_mix = P_mix.sortby('hexint')

# Convert hexint values to hex grid and create a hexbin grid with resolution 3
hex_grid = hexfunc.int_to_hex(P_mix.hexint.values)
hexbin_grid = hexfunc.hexGrid(hex_grid, h3_res=3)

all_ps = {}

for delta_r in [2., 0.1]:
    file_path = f"/storage/shared/oceanparcels/output_data/data_Claudio/NEMO_Ensemble/analysis/prob_distribution/{location}_spatial_long/P_dr{delta_r*100:03.0f}_m{member:03d}.nc"
    P_m = xr.open_dataset(file_path)
    P_m = P_m.sortby('hexint')
    all_ps[delta_r] = P_m

    
for week in [20, 4]:
    file_path = f"/storage/shared/oceanparcels/output_data/data_Claudio/NEMO_Ensemble/analysis/prob_distribution/{location}_temporal_long/P_W{week:01d}_m{member:03d}.nc"
    P_m = xr.open_dataset(file_path)
    P_m = P_m.sortby('hexint')
    all_ps[week] = P_m

# %% Plot 100 days snapshots STD
Latitude_limit = 53
extent = [-85, -50, 20, 45]
ncol = 3
nrow = 5
fig, axs = plt.subplots(ncols=ncol, nrows=nrow, figsize=(8.27, 7.),
                        subplot_kw={'projection': cartopy.crs.PlateCarree()},
                        constrained_layout=True,
                        gridspec_kw={'height_ratios': [1, 1, 1, 1, 0.4]})

axs = axs.reshape(ncol * nrow)

colr_mapa = 'plasma'

hexbin_grid.pcolorhex(P_mix['probability'][:, 10], ax=axs[0], cmap=colr_mapa, draw_edges=False, maxnorm=0.1)
hexbin_grid.pcolorhex(P_mix['probability'][:, 100], ax=axs[1], cmap=colr_mapa, draw_edges=False, maxnorm=0.1)
im = hexbin_grid.pcolorhex(P_mix['probability'][:, 1000], ax=axs[2], cmap=colr_mapa, draw_edges=False, maxnorm=0.1)

hexbin_grid.pcolorhex(all_ps[20]['probability'][:, 10], ax=axs[3], cmap=colr_mapa, draw_edges=False, maxnorm=0.1)
hexbin_grid.pcolorhex(all_ps[20]['probability'][:, 100], ax=axs[4], cmap=colr_mapa, draw_edges=False, maxnorm=0.1)
hexbin_grid.pcolorhex(all_ps[20]['probability'][:, 1000], ax=axs[5], cmap=colr_mapa, draw_edges=False, maxnorm=0.1)

hexbin_grid.pcolorhex(all_ps[2.]['probability'][:, 10], ax=axs[6], cmap=colr_mapa, draw_edges=False, maxnorm=0.1)
hexbin_grid.pcolorhex(all_ps[2.]['probability'][:, 100], ax=axs[7], cmap=colr_mapa, draw_edges=False, maxnorm=0.1)
hexbin_grid.pcolorhex(all_ps[2.]['probability'][:, 1000], ax=axs[8], cmap=colr_mapa, draw_edges=False, maxnorm=0.1)

hexbin_grid.pcolorhex(all_ps[0.1]['probability'][:, 10], ax=axs[9], cmap=colr_mapa, draw_edges=False, maxnorm=0.1)
hexbin_grid.pcolorhex(all_ps[0.1]['probability'][:, 100], ax=axs[10], cmap=colr_mapa, draw_edges=False, maxnorm=0.1)
hexbin_grid.pcolorhex(all_ps[0.1]['probability'][:, 1000], ax=axs[11], cmap=colr_mapa, draw_edges=False, maxnorm=0.1)

# formating

for i in range(0, ncol * nrow - ncol):
    axs[i].set_extent([-90, -20, 10, 45], crs=cartopy.crs.PlateCarree())

    axs[i].add_feature(cartopy.feature.LAND, zorder=10, color='black')

    if i not in [0, 3, 6, 9]:
        axs[i].scatter(-73.6, 35.6, color="deepskyblue", edgecolor='black', s=30,
                       marker="o", label="Release Location", zorder=10)

    gl = axs[i].gridlines(crs=cartopy.crs.PlateCarree(), draw_labels=True,
                          linewidth=0.5, color='gray', alpha=0.3)

    gl.top_labels = False
    gl.left_labels = False
    gl.right_labels = False
    gl.bottom_labels = False

    if i in [2, 5, 8, 11]:
        gl.right_labels = True

    if i in [9, 10, 11]:
        gl.bottom_labels = True

    # Add labels like A, B, C, etc.
    label = chr(65 + i)  # 65 is the ASCII value for 'A'
    axs[i].text(0.02, 0.98, label, transform=axs[i].transAxes, fontsize=12, fontweight='bold', va='bottom', ha='left', color='w')

# Add labels to the left of axs 0, 3, 6, 9
axs[0].text(-0.01, 0.5, r'Mixture $\delta_r = 0.1^o$', transform=axs[0].transAxes, fontsize=12, fontweight='bold', va='center', ha='right', rotation='vertical')
axs[3].text(-0.01, 0.5, '20 weeks', transform=axs[3].transAxes, fontsize=12, va='center', ha='right', rotation='vertical')
axs[6].text(-0.01, 0.5, r'$\delta_r = 2.0^o$', transform=axs[6].transAxes, fontsize=12, va='center', ha='right', rotation='vertical')
axs[9].text(-0.01, 0.5, r'$\delta_r = 0.1^o$', transform=axs[9].transAxes, fontsize=12, va='center', ha='right', rotation='vertical')

for i in range(12, 15):
    axs[i].axis('off')

axs[0].set_title(f'Particle Age 10 days')
axs[1].set_title(f'Particle Age 100 days')
axs[2].set_title(f'Particle Age 1000 days')

colorbar_axis = fig.add_axes([0.1, 0.06, 0.8, 0.03])  # Adjust the height of the colorbar
colorbar = fig.colorbar(im, cax=colorbar_axis,
                        orientation='horizontal',
                        label=f'Probability', extend='max')

plt.savefig(f'../figs/FigS9_maps_snapshots_m{member:03d}_s{subset:03d}.png', dpi=300)
 # %%
