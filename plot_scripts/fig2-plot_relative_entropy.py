#%%
import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
import cartopy
from tqdm import tqdm
import pickle

import sys
sys.path.append('../functions')
import hexbin_functions as hexfunc

def entropy(Pdf):
    # Pdf = Pdf / np.nansum(Pdf)  # Normalize Pdf to sum to 1, ignoring NaNs
    # Replace zeros with a very small number to avoid log(0)
    Pdf_safe = np.where(Pdf > 0, Pdf, np.finfo(float).eps)
    return -np.nansum(Pdf_safe * np.log(Pdf_safe))

def information(Pdf):
    Pdf_safe = np.where(Pdf > 0, Pdf, np.finfo(float).eps)
    return - Pdf_safe * np.log2(Pdf_safe)

def relative_information(P, Q):
    P_safe = np.where(P > 0, P, np.finfo(float).eps)
    Q_safe = np.where(Q > 0, Q, np.finfo(float).eps)
    return Q_safe*np.log2(Q_safe/ P_safe)

def Shannon_entropy(Pdf):
    # Pdf = Pdf / np.nansum(Pdf)  # Normalize Pdf to sum to 1, ignoring NaNs
    # Replace zeros with a very small number to avoid log(0)
    Pdf_safe = np.where(Pdf > 0, Pdf, np.finfo(float).eps)
    return -np.nansum(Pdf_safe * np.log2(Pdf_safe))

location = 'Cape_Hatteras'
delta_r = 0.1
subset = 13

file_path_AX = f"/storage/shared/oceanparcels/output_data/data_Claudio/NEMO_Ensemble/analysis/prob_distribution/{location}_all_long/P_dr{delta_r*100:03.0f}_all_s{subset:03d}.nc"
P_AX = xr.open_dataset(file_path_AX)
P_AX = P_AX.sortby('hexint')

member = 13
delta_r = 0.1
file_path_M = f"/storage/shared/oceanparcels/output_data/data_Claudio/NEMO_Ensemble/analysis/prob_distribution/{location}_spatial_long/P_dr{delta_r*100:03.0f}_m{member:03d}.nc"
P_m = xr.open_dataset(file_path_M)
P_m = P_m.sortby('hexint')


P_AX['hexint'] == P_m['hexint']

hex_grid = hexfunc.int_to_hex(P_AX.hexint.values)
hexbin_grid = hexfunc.hexGrid(hex_grid, h3_res=3)

ncol = 4
nrow = 2
fig, axs = plt.subplots(ncols=ncol, nrows=nrow, figsize=(9, 3),
                       subplot_kw={'projection': cartopy.crs.PlateCarree()},
                       sharey=True, constrained_layout=True, gridspec_kw={'height_ratios': [0.8, 0.2]})

axs = axs.reshape(ncol*nrow)
t = 30
extent = [-85, -65, 30, 40]

# Function to set up each subplot
def setup_subplot(ax, extent, gridlines=True, left_labels=True, right_labels=False, bottom_labels=True):
    ax.set_extent(extent, crs=cartopy.crs.PlateCarree())
    ax.add_feature(cartopy.feature.LAND, zorder=0, color='gray')
    if gridlines:
        gl = ax.gridlines(crs=cartopy.crs.PlateCarree(), draw_labels=True,
                          linewidth=0.5, color='gray', alpha=0.)
        gl.top_labels = False
        gl.right_labels = right_labels
        gl.left_labels = left_labels
        gl.bottom_labels = bottom_labels

# Plot mixture probability
setup_subplot(axs[0], extent)
mixture = np.where(P_AX['probability'][:, t].values > 0, P_AX['probability'][:, t].values, np.finfo(float).eps)
single = np.where(P_m['probability'][:, t].values > 0, P_m['probability'][:, t].values, np.finfo(float).eps)
max_value = np.max(mixture)

im = hexbin_grid.pcolorhex(mixture, ax=axs[0], cmap='plasma', draw_edges=True, maxnorm=max_value)

# Plot single probability
setup_subplot(axs[1], extent, left_labels=False)
hexbin_grid.pcolorhex(single, ax=axs[1], cmap='plasma', draw_edges=True, maxnorm=max_value)

# Plot relative information
setup_subplot(axs[2], extent, left_labels=False)
rel_inf = relative_information(single, mixture)
rel_inf_r = relative_information(mixture, single)
max_rel_inf = np.max(rel_inf_r)
min_rel_inf = min(np.min(rel_inf), np.min(rel_inf_r))

im2 = hexbin_grid.pcolorhex(rel_inf, ax=axs[2], cmap='viridis', draw_edges=True, maxnorm=max_rel_inf, minnorm=min_rel_inf,
                            negative=False)

# Plot relative information of P_AX, P_m
setup_subplot(axs[3], extent, left_labels=False, right_labels=True)

hexbin_grid.pcolorhex(rel_inf_r, ax=axs[3], cmap='viridis', draw_edges=True, maxnorm=max_rel_inf, minnorm=min_rel_inf)

# Turn off unused subplots
for i in range(4, 8):
    axs[i].axis('off')

# Add colorbars
# Add labels to subplots
labels = ['A', 'B', 'C', 'D']
for i, label in enumerate(labels):
    axs[i].text(0.03, 0.95, label, transform=axs[i].transAxes, fontsize=12, fontweight='bold',
                verticalalignment='top', horizontalalignment='left', bbox=dict(facecolor='none', alpha=0))

# Add colorbars
bar_ax = fig.add_axes([0.1, 0.30, 0.35, 0.05])
cbar = fig.colorbar(im, cax=bar_ax, orientation='horizontal', label='Probability', extend='max')

bar_ax2 = fig.add_axes([0.55, 0.3, 0.35, 0.05])
cbar2 = fig.colorbar(im2, cax=bar_ax2, orientation='horizontal', label='Information Loss (bits)', extend='max')
# Add titles to each subplot
axs[0].set_title(r'$P_{mix}(X| \delta_r = 0.1^o)$', fontsize=10)
axs[1].set_title(r'$P_{m}(X| \delta_r = 0.1^o)$', fontsize=10)
axs[2].set_title(r'$D(P_{mix}||P_m) = ' + f'{np.sum(rel_inf):0.1f}$ bits', fontsize=10)
axs[3].set_title(r'$D(P_{m}||P_{mix}) = ' + f'{np.sum(rel_inf_r):0.1f}$ bits', fontsize=10)

# Save figure
plt.savefig(f'../figs/Fig2_explanations_relative_entropy.png', dpi=300)
#%%