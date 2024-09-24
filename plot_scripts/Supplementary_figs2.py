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

#%%%%%%%%%%%%%%%%%%%%%%%% SPATIAL ANALYSIS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
################################################################################
entropies_space_mean = {}
entropies_space_std = {}

for i, delta_r in enumerate([0.1, 1, 2]):
    _entropy = np.zeros((50, len(P_AX['entropy'].values)))
    for i, member in enumerate(range(1, 51)):
        file_path = f"/storage/shared/oceanparcels/output_data/data_Claudio/NEMO_Ensemble/analysis/prob_distribution/{location}_spatial_long/P_dr{delta_r*100:03.0f}_m{member:03d}.nc"
        P_m = xr.open_dataset(file_path)
        
        _entropy[i, :] = P_m['entropy'].values
        
    entropies_space_mean[delta_r] = np.mean(_entropy, axis=0)
    entropies_space_std[delta_r] = np.std(_entropy, axis=0)

#%% Plots Spatial entropy

member_list = range(1, 51)
ncol = 3
nrow = 1
fig, axs = plt.subplots(ncols=ncol, nrows=nrow, figsize=(8, 3),
                          sharey=True, constrained_layout=True)

axs = axs.reshape(ncol*nrow)

for i, delta_r in enumerate([0.1, 1, 2]):
        
    for member in member_list:
        file_path = f"/storage/shared/oceanparcels/output_data/data_Claudio/NEMO_Ensemble/analysis/prob_distribution/{location}_spatial_long/P_dr{delta_r*100:03.0f}_m{member:03d}.nc"
        P_m = xr.open_dataset(file_path)
        axs[i].plot(P_m['entropy'])
    
    axs[i].plot(entropies_space_mean[delta_r], ls='--', color='black', label='Mean')
    axs[i].grid()
    axs[i].set_title(f'$\delta r$ = {delta_r}$^o$')
    axs[i].set_xlabel('Time (days)')
    axs[i].set_xlim(0, 2189)

axs[0].set_ylabel('Shannon Entropy (bits)')
axs[0].legend()
plt.savefig('../figs/FigS1-Spatial-Representation_entropy_all.png', dpi=300)

# %%%%%%%%%%%%%%%%%%%%%% TEMPORAL Analysis %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
################################################################################
week_range = [4, 12, 20]


entropies_time_mean = {}
entropies_time_std = {}

for week in week_range:
    _entropy = np.zeros((50, len(P_AX['entropy'].values)))
    for i, member in enumerate(range(1, 51)):
        file_path = f"/storage/shared/oceanparcels/output_data/data_Claudio/NEMO_Ensemble/analysis/prob_distribution/{location}_temporal_long/P_W{week:01d}_m{member:03d}.nc"
        P_m = xr.open_dataset(file_path)
        _entropy[i, :] = P_m['entropy'].values
        
    entropies_time_mean[week] = np.mean(_entropy, axis=0)
    entropies_time_std[week] = np.std(_entropy, axis=0)
    
    
#%% Plots Temporal entropy

member_list = range(1, 51)
ncol = 3
nrow = 1
fig, axs = plt.subplots(ncols=ncol, nrows=nrow, figsize=(8, 3),
                       sharey=True, constrained_layout=True)

axs = axs.reshape(ncol*nrow)

for i, week in enumerate(week_range):
    
    for member in member_list:
        file_path = f"/storage/shared/oceanparcels/output_data/data_Claudio/NEMO_Ensemble/analysis/prob_distribution/{location}_temporal_long/P_W{week:01d}_m{member:03d}.nc"
        P_m = xr.open_dataset(file_path)
        chop_time = len(P_m['time'].values) - week*7
        axs[i].plot(P_m['entropy'][:chop_time])
    
    axs[i].plot(entropies_time_mean[week][:chop_time], ls='--', color='black', label='Mean')
    axs[i].grid()
    axs[i].set_title(f'Release span {week} weeks')
    axs[i].set_xlim(0, chop_time)
    # axs[i].set_xlabel('Time (days)')
    
axs[0].set_ylabel('Shannon Entropy (bits)')
axs[1].set_xlabel('Time (days)')
axs[2].set_xlabel('Time (days)')
axs[0].set_xlabel('Time (days)')
axs[0].legend()
plt.savefig('../figs/FigS2-Temporal-Representation_entropy_all.png', dpi=300)

#%% 