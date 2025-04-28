#%%
import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
# import cartopy
# from tqdm import tqdm
# import pickle

import sys
sys.path.append('../functions')
import hexbin_functions as hexfunc


location = 'Cape_Hatteras'
delta_r = 0.1
subset = 1

file_path_AX = f"/storage/shared/oceanparcels/output_data/data_Claudio/NEMO_Ensemble/analysis/prob_distribution/{location}_all_long/P_dr{delta_r*100:03.0f}_all_s{subset:03d}.nc"
P_AX = xr.open_dataset(file_path_AX)

hex_grid = hexfunc.int_to_hex(P_AX.hexint.values)
hexbin_grid = hexfunc.hexGrid(hex_grid, h3_res=3)


# %% Average the subsets space

mixture_entropy_space = {}
std_mixture_entropy_space = {}

for delta_r in [0.1, 2]:
    entropy_all = np.zeros((51, len(P_AX['entropy'].values)))
    for i in range(1, 51):
        
        file_path_AX = f"/storage/shared/oceanparcels/output_data/data_Claudio/NEMO_Ensemble/analysis/prob_distribution/{location}_all_long/P_dr{delta_r*100:03.0f}_all_s{i:03d}.nc"
        P_AX = xr.open_dataset(file_path_AX)
        
        entropy_all[i-1, :] = P_AX['entropy'].values

    mixture_entropy_space[delta_r] = np.mean(entropy_all, axis=0)
    std_mixture_entropy_space[delta_r] = np.std(entropy_all, axis=0)

# %% Average the subsets TIME

mixture_entropy_time = {}
std_mixture_entropy_time = {}

for week in [4, 20]:
    entropy_all = np.zeros((50, len(P_AX['entropy'].values)))
    
    for i in range(1, 50):
        
        file_path_AX = f"/storage/shared/oceanparcels/output_data/data_Claudio/NEMO_Ensemble/analysis/prob_distribution/{location}_all_long/P_W{week:02d}_all_s{i:03d}.nc"
        P_AX = xr.open_dataset(file_path_AX)
        
        entropy_all[i-1, :] = P_AX['entropy'].values

    mixture_entropy_time[week] = np.mean(entropy_all, axis=0)
    std_mixture_entropy_time[week] = np.std(entropy_all, axis=0)


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
    

# %%%%%%%%%%%%%%%%%%%%%% Diffusion Analysis %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
################################################################################
K_h_range = [10]

entropies_diff_mean = {}
entropies_diff_std = {}

for K_h in K_h_range:
    _entropy = np.zeros((50, len(P_AX['entropy'].values)))
    for i, member in enumerate(range(1, 51)):
        file_path = f"/storage/shared/oceanparcels/output_data/data_Claudio/NEMO_Ensemble/analysis/prob_distribution/{location}_diffusion_long/P_diff_Kh_{K_h:01d}_m{member:03d}.nc"
        P_m = xr.open_dataset(file_path)
        P_m = P_m.isel(time=slice(0, 2189))
        
        _entropy[i, :] = P_m['entropy']
        
    entropies_diff_mean[K_h] = np.mean(_entropy, axis=0)
    entropies_diff_std[K_h] = np.std(_entropy, axis=0)

# %% spatial plot of the mean entropy
time_range = np.arange(0, 2189)

fig, ax = plt.subplots()

lss = [(0, (1, 1)), '--', '-.']
i=0

colors_space = ['mediumblue', 'blueviolet', 'teal']
for delta_r in [0.1, 1., 2.]:
    
    ax.plot(time_range, entropies_space_mean[delta_r],  
            label=f'$\delta_r = {delta_r}^o$', linestyle=lss[i], linewidth=2, color=colors_space[i])
    ax.fill_between(time_range, entropies_space_mean[delta_r] - entropies_space_std[delta_r], 
                    entropies_space_mean[delta_r] + entropies_space_std[delta_r], 
                    alpha=0.2, color=colors_space[i])
    i+=1

ax.plot(time_range, mixture_entropy_space[0.1], ls='-', color='black', label=r'Mixture: $\delta_r = 0.1^o$')
ax.plot(time_range, mixture_entropy_space[2], ls=(0, (3, 1, 1, 1)), color='black', label=r'Mixture: $\delta_r = 2.0^o$')

# ax.semilogx(time_range[:10],  time_range[:10]**(3/8) + 3, ls='-', color='red')
ax.plot(time_range[10:1000],  1*np.log(time_range[10:1000]), ls='-', color='red')
ax.plot(time_range[1:1000],  2*np.log(time_range[1:1000]), ls='-', color='red')
ax.plot(time_range[10:1000],  (time_range[10:1000])**1, ls='-', color='green')
ax.set_xlim(1, 2189)
ax.set_ylim(0., 10.5)
ax.semilogx(base=10)
ax.legend(shadow=True)
ax.set_ylabel('Marginal Entropy (bits)')
ax.set_xlabel('Particle Age (days)')
ax.grid()
# plt.savefig('../figs/Fig3-Spatial-Representation_entropy.png', dpi=300)

# %%%%%%%%%%%%%%%%%%%%%%%%% Temporal plot of the mean entropy %%%%%%%%%%%%%%%%%%%%%%%%%

fig, ax = plt.subplots()

lss = [(0, (1, 1)), '--', '-.', (0, (3, 1, 1, 1, 1, 1))]
i=0
colors_time = ['darkred', 'orangered', 'orange']

for week in week_range:
    
    chop_time = len(P_m['time'].values) - week*7
    ax.plot(time_range[:chop_time], entropies_time_mean[week][:chop_time],  
            label=f'{week} weeks release', linestyle=lss[i], linewidth=2, color=colors_time[i])
    ax.fill_between(time_range[:chop_time], entropies_time_mean[week][:chop_time] - entropies_time_std[week][:chop_time], 
                    entropies_time_mean[week][:chop_time] + entropies_time_std[week][:chop_time],
                    alpha=0.2, color=colors_time[i])
    i+=1

ax.plot(time_range[:2189-4*7], mixture_entropy_time[4][:2189-4*7], ls=(0, (3, 1, 1, 3)), 
        color='black', label='Mixture: 4 weeks')

# ax.plot(time_range[:2189-20*7], mixture_entropy_time[20][:2189-20*7], ls=(0, (3, 1, 1, 5)), 
#         color='black', label='Mixture: 20 weeks')

# ax.plot(time_range[:730-16*7], mixture_entropy_time[16][:730-16*7], ls='-', 
#         color='black', label='Mixture: 16 weeks')

ax.plot(time_range, mixture_entropy_space[0.1], ls='-', color='black', label=r'Mixture: $\delta_r = 0.1^o$')

ax.set_xlim(1, 2189)
ax.set_ylim(0., 10.5)

ax.semilogx()
ax.legend(shadow=True)
ax.set_ylabel('Marginal Entropy (bits)')
ax.set_xlabel('Particle Age (days)')
ax.grid()
# plt.savefig('../figs/Fig4-Temporal-Representation_entropy.png', dpi=300)


#%%
# %% Combined spatial and temporal plot of the mean entropy
time_range = np.arange(0, 2189)

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4), sharey=True)

# Spatial plot
lss = [(0, (1, 1)), '--', '-.']
i = 0
colors_space = ['mediumblue', 'blueviolet', 'teal']
for delta_r in [0.1, 1., 2.]:
    ax1.plot(time_range, entropies_space_mean[delta_r],  
             label=f'$\delta_r = {delta_r}^o$', linestyle=lss[i], linewidth=2, color=colors_space[i])
    ax1.fill_between(time_range, entropies_space_mean[delta_r] - entropies_space_std[delta_r], 
                     entropies_space_mean[delta_r] + entropies_space_std[delta_r], 
                     alpha=0.2, color=colors_space[i])
    i += 1


ax1.plot(time_range, mixture_entropy_space[0.1], ls='-', color='black', label=r'Mixture: $\delta_r = 0.1^o$')
ax1.plot(time_range, mixture_entropy_space[2], ls=(0, (3, 1, 1, 1)), color='black', label=r'Mixture: $\delta_r = 2.0^o$')

ax1.set_xlim(1, 2189)
ax1.set_ylim(0., 10.5)
ax1.semilogx()
ax1.legend(shadow=True)
ax1.set_ylabel('Marginal Entropy (bits)')
ax1.set_xlabel('Particle Age (days)')
ax1.grid()
ax1.text(0.95, 0.05, 'A', transform=ax1.transAxes, fontsize=16, fontweight='bold', va='bottom', ha='right')

# Temporal plot
lss = [(0, (1, 1)), '--', '-.', (0, (3, 1, 1, 1, 1, 1))]
i = 0
colors_time = ['darkred', 'orangered', 'orange']

for week in week_range:
    chop_time = len(P_m['time'].values) - week*7
    ax2.plot(time_range[:chop_time], entropies_time_mean[week][:chop_time],  
             label=f'{week} weeks release', linestyle=lss[i], linewidth=2, color=colors_time[i])
    ax2.fill_between(time_range[:chop_time], entropies_time_mean[week][:chop_time] - entropies_time_std[week][:chop_time], 
                     entropies_time_mean[week][:chop_time] + entropies_time_std[week][:chop_time],
                     alpha=0.2, color=colors_time[i])
    i += 1

# diffusion
K_h = 10
# Create logarithmically spaced indices
# This will sample more points at the beginning and fewer at the end
n_samples = 100  # Choose the number of points you want
log_indices = np.unique(np.geomspace(1, len(time_range)-1, n_samples).astype(int))

# Use these indices to sample your data
time_sampled = time_range[log_indices]
y_sampled = entropies_diff_mean[K_h][log_indices]

# ax.plot(time_range[:], entropies_diff_mean[K_h][:],  
#         label=r'$K_h = 10 \ m^2 s^{-1}$', linestyle=(0, (2, 2, 10, 2)) , linewidth=2, color=color_diff)

ax2.fill_between(time_range[:], entropies_diff_mean[K_h][:] - entropies_diff_std[K_h][:], 
                entropies_diff_mean[K_h][:] + entropies_diff_std[K_h][:],
                alpha=0.2, color='grey')
ax2.scatter(time_sampled, y_sampled,
        label=r'$K_h = 10 \ m^2 s^{-1}$', color='k', marker='^', s=10)

ax2.plot(time_range, mixture_entropy_space[0.1], ls='-', color='black', label=r'Mixture: $\delta_r = 0.1^o$')
ax2.plot(time_range[:2189-4*7], mixture_entropy_time[4][:2189-4*7], ls=(0, (3, 1, 1, 3)), 
         color='black', label='Mixture: 4 weeks')
ax2.plot(time_range[:2189-20*7], mixture_entropy_time[20][:2189-20*7], ls=(0, (6, 4, 2, 2)), 
         color='black', label='Mixture: 20 weeks')



ax2.set_xlim(1, 2189)
ax2.set_ylim(0., 10.5)
ax2.semilogx()
ax2.legend(shadow=True)
ax2.set_xlabel('Particle Age (days)')
ax2.grid()
ax2.text(0.95, 0.05, 'B', transform=ax2.transAxes, fontsize=16, fontweight='bold', va='bottom', ha='right')

plt.tight_layout()
plt.savefig('../figs/Fig4_Combined-Representation_entropy.png', dpi=300)

# %% With diffusion

fig, ax = plt.subplots()

# Spatial plot
delta_r = 2.0
ax.plot(time_range, entropies_space_mean[delta_r],  
            label=f'$\delta_r = {delta_r}^o$', linestyle='-.', linewidth=2, color='teal')
ax.fill_between(time_range, entropies_space_mean[delta_r] - entropies_space_std[delta_r], 
                    entropies_space_mean[delta_r] + entropies_space_std[delta_r], 
                    alpha=0.2, color="teal")


# time
week = 20
chop_time = len(P_m['time'].values) - week*7
ax.plot(time_range[:chop_time], entropies_time_mean[week][:chop_time],  
        label=f'{week} weeks release', linestyle=(0, (3, 1, 1, 1, 1, 1)), linewidth=2, color='orange')
ax.fill_between(time_range[:chop_time], entropies_time_mean[week][:chop_time] - entropies_time_std[week][:chop_time], 
                entropies_time_mean[week][:chop_time] + entropies_time_std[week][:chop_time],
                alpha=0.2, color='orange')

# diffusion
K_h = 10
# Create logarithmically spaced indices
# This will sample more points at the beginning and fewer at the end
n_samples = 150  # Choose the number of points you want
log_indices = np.unique(np.geomspace(1, len(time_range)-1, n_samples).astype(int))

# Use these indices to sample your data
time_sampled = time_range[log_indices]
y_sampled = entropies_diff_mean[K_h][log_indices]

# ax.plot(time_range[:], entropies_diff_mean[K_h][:],  
#         label=r'$K_h = 10 \ m^2 s^{-1}$', linestyle=(0, (2, 2, 10, 2)) , linewidth=2, color=color_diff)

ax.fill_between(time_range[:], entropies_diff_mean[K_h][:] - entropies_diff_std[K_h][:], 
                entropies_diff_mean[K_h][:] + entropies_diff_std[K_h][:],
                alpha=0.5, color='grey')
ax.scatter(time_sampled, y_sampled,
        label=r'$K_h = 10 \ m^2 s^{-1}$', color='k', marker=',', s=1)


ax.plot(time_range, mixture_entropy_space[0.1], ls='-', color='black', label=r'Mixture: $\delta_r = 0.1^o$')

ax.set_xlim(1, 2189)
ax.set_ylim(0., 10.5)

ax.semilogx()
ax.legend(shadow=True)
ax.set_ylabel('Marginal Entropy (bits)')
ax.set_xlabel('Particle Age (days)')
ax.grid()
# plt.savefig('../figs/Fig4-Temporal-Represen
# %%
