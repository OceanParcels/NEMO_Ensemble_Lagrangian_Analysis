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

location = 'Cape_Hatteras'
delta_r = 0.1
subset = 1

file_path_AX = f"/storage/shared/oceanparcels/output_data/data_Claudio/NEMO_Ensemble/analysis/prob_distribution/{location}_all/P_dr{delta_r*100:03.0f}_all_s{subset}.nc"
P_AX = xr.open_dataset(file_path_AX)

hex_grid = hexfunc.int_to_hex(P_AX.hexint.values)
hexbin_grid = hexfunc.hexGrid(hex_grid, h3_res=3)


# %% Average the subsets space

mixture_entropy_space = {}
std_mixture_entropy_space = {}

for delta_r in [0.1, 1, 2]:
    entropy_all = np.zeros((50, len(P_AX['entropy'].values)))
    for i in range(1, 51):
        
        file_path_AX = f"/storage/shared/oceanparcels/output_data/data_Claudio/NEMO_Ensemble/analysis/prob_distribution/{location}_all/P_dr{delta_r*100:03.0f}_all_s{i}.nc"
        P_AX = xr.open_dataset(file_path_AX)
        
        entropy_all[i-1, :] = P_AX['entropy'].values

    mixture_entropy_space[delta_r] = np.mean(entropy_all, axis=0)
    std_mixture_entropy_space[delta_r] = np.std(entropy_all, axis=0)

# %% Average the subsets TIME

mixture_entropy_time = {}
std_mixture_entropy_time = {}

for week in [4,16]:
    entropy_all = np.zeros((50, len(P_AX['entropy'].values)))
    
    for i in range(1, 51):
        
        file_path_AX = f"/storage/shared/oceanparcels/output_data/data_Claudio/NEMO_Ensemble/analysis/prob_distribution/{location}_all/P_W{week:01d}_all_s{i}.nc"
        P_AX = xr.open_dataset(file_path_AX)
        
        entropy_all[i-1, :] = P_AX['entropy'].values

    mixture_entropy_time[week] = np.mean(entropy_all, axis=0)
    std_mixture_entropy_time[week] = np.std(entropy_all, axis=0)


# %%
for delta_r in [0.1, 1, 2]:
    plt.plot(mixture_entropy_space[delta_r], label=f'$\delta_r={delta_r}^o$')

for week in [4, 16]:
    plt.plot(mixture_entropy_time[week], label=f'W={week}')

plt.semilogx()
plt.legend()
plt.ylabel('Shannon Entropy (bits)')
plt.xlabel('Time (days)')
plt.grid()
plt.title('Mean Shannon Entropy of the Subsets of All Members')


# %%%%%%%%%%%%%%%%%%%%%% TEMPORAL Analysis %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
################################################################################
week_range = [4,8,12, 16]

member_list = range(1, 51)
ncol = 2
nrow = 2
fig, axs = plt.subplots(ncols=ncol, nrows=nrow, figsize=(6, 6),
                       sharey=True, constrained_layout=True)

axs = axs.reshape(ncol*nrow)

for i, week in enumerate(week_range):
    
    for member in member_list:
        file_path = f"/storage/shared/oceanparcels/output_data/data_Claudio/NEMO_Ensemble/analysis/prob_distribution/{location}_temporal/P_W{week:01d}_m{member:03d}.nc"
        P_m = xr.open_dataset(file_path)
        chop_time = len(P_m['time'].values) - week*7
        axs[i].plot(P_m['entropy'][:chop_time])
    
    axs[i].plot(mixture_entropy_time[4][:702], ls='--', color='black', label='Mean All Members')
    axs[i].grid()
    axs[i].set_title(f'Release span {week} weeks')
    # axs[i].set_xlabel('Time (days)')
    
axs[0].set_ylabel('Shannon Entropy (bits)')
axs[2].set_ylabel('Shannon Entropy (bits)')
axs[2].set_xlabel('Time (days)')
axs[3].set_xlabel('Time (days)')
plt.savefig('../figs/FigS2-Temporal-Representation_entropy_all.png', dpi=300)

# %%%%%%%%%%%%%%%%%%%%%%% Average the subsets TIME %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

entropies_time_mean = {}
entropies_time_std = {}

for week in week_range:
    _entropy = np.zeros((50, len(P_AX['entropy'].values)))
    for i, member in enumerate(range(1, 51)):
        file_path = f"/storage/shared/oceanparcels/output_data/data_Claudio/NEMO_Ensemble/analysis/prob_distribution/{location}_temporal/P_W{week:01d}_m{member:03d}.nc"
        P_m = xr.open_dataset(file_path)
        _entropy[i, :] = P_m['entropy'].values
        
    entropies_time_mean[week] = np.mean(_entropy, axis=0)
    entropies_time_std[week] = np.std(_entropy, axis=0)

# %%%%%%%%%%%%%%%%%%%%%%%%% Temporal plot of the mean entropy %%%%%%%%%%%%%%%%%%%%%%%%%
time_range = np.arange(0, 730)
fig, ax = plt.subplots(figsize=(8, 6))

lss = [(0, (1, 1)), '--', '-.', (0, (3, 1, 1, 1, 1, 1))]
i=0
for week in week_range:
    
    chop_time = len(P_m['time'].values) - week*7
    ax.plot(time_range[:chop_time], entropies_time_mean[week][:chop_time],  
            label=f'{week} weeks release', linestyle=lss[i], linewidth=2)
    ax.fill_between(time_range[:chop_time], entropies_time_mean[week][:chop_time] - entropies_time_std[week][:chop_time], 
                    entropies_time_mean[week][:chop_time] + entropies_time_std[week][:chop_time], alpha=0.2)
    i+=1

ax.plot(time_range[:730-4*7], mixture_entropy_time[4][:730-4*7], ls=(0, (3, 1, 1, 3)), 
        color='black', label='Mixture: 4 weeks')

# ax.plot(time_range[:730-16*7], mixture_entropy_time[16][:730-16*7], ls='-', 
#         color='black', label='Mixture: 16 weeks')

ax.plot(time_range, mixture_entropy_space[0.1], ls='-', color='black', label=r'Mixture: $\delta_r = 0.1^o$')

ax.set_xlim(1, 1000)

ax.semilogx()
ax.legend()
ax.set_ylabel('Shannon Entropy (bits)')
ax.set_xlabel('Time (days)')
ax.grid()
plt.savefig('../figs/Fig4-Temporal-Representation_entropy.png', dpi=300)


#%%%%%%%%%%%%%%%%%%%%%%%% SPATIAL ANALYSIS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
################################################################################
entropies_space_mean = {}
entropies_space_std = {}

for i, delta_r in enumerate([0.1, 1, 2]):
    _entropy = np.zeros((50, len(P_AX['entropy'].values)))
    for i, member in enumerate(range(1, 51)):
        file_path = f"/storage/shared/oceanparcels/output_data/data_Claudio/NEMO_Ensemble/analysis/prob_distribution/{location}_spatial/P_dr{delta_r*100:03.0f}_m{member:03d}.nc"
        P_m = xr.open_dataset(file_path)
        
        _entropy[i, :] = P_m['entropy'].values
        
    entropies_space_mean[delta_r] = np.mean(_entropy, axis=0)
    entropies_space_std[delta_r] = np.std(_entropy, axis=0)


# %% spatial plot of the mean entropy
fig, ax = plt.subplots(figsize=(8, 6))



lss = [(0, (1, 1)), '--', '-.']
i=0
for delta_r in [0.1, 1., 2.]:
    
    ax.plot(time_range, entropies_space_mean[delta_r],  
            label=f'$\delta_r = {delta_r}^o$', linestyle=lss[i], linewidth=2)
    ax.fill_between(time_range, entropies_space_mean[delta_r] - entropies_space_std[delta_r], 
                    entropies_space_mean[delta_r] + entropies_space_std[delta_r], alpha=0.2)
    i+=1

ax.plot(time_range, mixture_entropy_space[0.1], ls='-', color='black', label=r'Mixture: $\delta_r = 0.1^o$')
ax.plot(time_range, mixture_entropy_space[2], ls=(0, (3, 1, 1, 1)), color='black', label=r'Mixture: $\delta_r = 2.0^o$')

ax.set_xlim(1, 1000)
ax.semilogx()
ax.legend()
ax.set_ylabel('Shannon Entropy (bits)')
ax.set_xlabel('Time (days)')
ax.grid()
plt.savefig('../figs/Fig3-Spatial-Representation_entropy.png', dpi=300)
# %%
