#%%
import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
# import cartopy
# from tqdm import tqdm
# import pickle

# import sys
# sys.path.append('../functions')
# import hexbin_functions as hexfunc

location = "Cape_Hatteras"
time_length = 2189
t_range_space = range(0, time_length)
t_range_time = range(0, time_length-7*20)

#%%%%%%%%%%%%%%%%%%%%%%%% SPATIAL ANALYSIS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
################################################################################
entropies_space_mean = {}
entropies_space_std = {}

for i, delta_r in enumerate([0.1, 1, 2]):
    _entropy = np.zeros((50, time_length))
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
    _entropy = np.zeros((50, time_length))
    for i, member in enumerate(range(1, 51)):
        file_path = f"/storage/shared/oceanparcels/output_data/data_Claudio/NEMO_Ensemble/analysis/prob_distribution/{location}_temporal_long/P_W{week:01d}_m{member:03d}.nc"
        P_m = xr.open_dataset(file_path)
        _entropy[i, :] = P_m['entropy'].values
        
    entropies_time_mean[week] = np.mean(_entropy, axis=0)
    entropies_time_std[week] = np.std(_entropy, axis=0)

#%% Combined Plots Spatial and Temporal entropy

member_list = range(1, 51)
ncol = 3
nrow = 2
fig, axs = plt.subplots(ncols=ncol, nrows=nrow, figsize=(10, 6),
                        sharey=False, constrained_layout=True)

axs = axs.reshape(ncol*nrow)

# Plot Spatial entropy
for i, delta_r in enumerate([0.1, 1, 2]):
    for member in member_list:
        file_path = f"/storage/shared/oceanparcels/output_data/data_Claudio/NEMO_Ensemble/analysis/prob_distribution/{location}_spatial_long/P_dr{delta_r*100:03.0f}_m{member:03d}.nc"
        P_m = xr.open_dataset(file_path)
        axs[i].plot(t_range_space, P_m['entropy'])
    
    axs[i].plot(t_range_space, entropies_space_mean[delta_r], ls='--', color='black', label='Mean')
    axs[i].grid()
    axs[i].set_title(f'$\delta r$ = {delta_r}$^o$')
    axs[i].set_xlabel('Time (days)')
    axs[i].set_xlim(0, 2189)
    axs[i].set_ylim(0, 11)
    
axs[0].set_ylabel('Entropy (bits)')
axs[0].legend()

# Plot Temporal entropy
for i, week in enumerate(week_range):
    for member in member_list:
        file_path = f"/storage/shared/oceanparcels/output_data/data_Claudio/NEMO_Ensemble/analysis/prob_distribution/{location}_temporal_long/P_W{week:01d}_m{member:03d}.nc"
        P_m = xr.open_dataset(file_path)
        chop_time = len(P_m['time'].values) - week*7
        axs[i + 3].plot(P_m['entropy'][:chop_time])
    
    axs[i + 3].plot(entropies_time_mean[week][:chop_time], ls='--', color='black', label='Mean')
    axs[i + 3].grid()
    axs[i + 3].set_title(f'Release span {week} weeks')
    axs[i + 3].set_xlim(0, chop_time)
    
axs[3].set_ylabel('Marginal Entropy (bits)')
axs[4].set_xlabel('Particle Age (days)')
axs[5].set_xlabel('Particle Age (days)')
axs[3].set_xlabel('Particle Age (days)')
axs[3].legend()

plt.savefig('../figs/FigS7-Combined-Representation_entropy_all.png', dpi=300)

#%% Mixture spatial
entropies_space_mix_mean = {}
entropies_space_mix_std = {}

for i, delta_r in enumerate([0.1, 1, 2]):
    _entropy = np.zeros((50, time_length))
    for i, member in enumerate(range(1, 51)):
        file_path = f"/storage/shared/oceanparcels/output_data/data_Claudio/NEMO_Ensemble/analysis/prob_distribution/{location}_all_long/P_dr{delta_r*100:03.0f}_all_s{member:03d}.nc"
        P_m = xr.open_dataset(file_path)
        
        _entropy[i, :] = P_m['entropy'].values
        
    entropies_space_mix_mean[delta_r] = np.mean(_entropy, axis=0)
    entropies_space_mix_std[delta_r] = np.std(_entropy, axis=0)


#%% Mixture temporal
entropies_time_mix_mean = {}
entropies_time_mix_std = {}

for week in week_range:
    _entropy = np.zeros((50, time_length))
    for i, member in enumerate(range(1, 51)):
        file_path = f"/storage/shared/oceanparcels/output_data/data_Claudio/NEMO_Ensemble/analysis/prob_distribution/{location}_all_long/P_W{week:02d}_all_s{member:03d}.nc"
        P_m = xr.open_dataset(file_path)
        _entropy[i, :] = P_m['entropy'].values
        
    entropies_time_mix_mean[week] = np.mean(_entropy, axis=0)
    entropies_time_mix_std[week] = np.std(_entropy, axis=0)
# %% Combined Plots mixture spatial and temporal entropy

member_list = range(1, 51)
ncol = 3
nrow = 2
fig, axs = plt.subplots(ncols=ncol, nrows=nrow, figsize=(10, 6),
                        sharey=False, constrained_layout=True)

axs = axs.reshape(ncol*nrow)

# Plot Mixture Spatial entropy
for i, delta_r in enumerate([0.1, 1, 2]):
    for member in member_list:
        file_path = f"/storage/shared/oceanparcels/output_data/data_Claudio/NEMO_Ensemble/analysis/prob_distribution/{location}_all_long/P_dr{delta_r*100:03.0f}_all_s{member:03d}.nc"
        P_m = xr.open_dataset(file_path)
        axs[i].plot(t_range_space, P_m['entropy'])
    
    axs[i].plot(t_range_space, entropies_space_mix_mean[delta_r], ls='--', color='black', label='Mean')
    axs[i].grid()
    axs[i].set_title(f'$\delta r$ = {delta_r}$^o$ (Mixture)')
    axs[i].set_xlabel('Time (days)')
    axs[i].set_xlim(0, 2189)
    axs[i].set_ylim(0, 11)
    # axs[i].semilogx()
    
axs[0].set_ylabel('Entropy (bits)')
axs[0].legend()

# Plot Mixture Temporal entropy
for i, week in enumerate(week_range):
    for member in member_list:
        file_path = f"/storage/shared/oceanparcels/output_data/data_Claudio/NEMO_Ensemble/analysis/prob_distribution/{location}_all_long/P_W{week:02d}_all_s{member:03d}.nc"
        P_m = xr.open_dataset(file_path)
        chop_time = len(P_m['time'].values) - week*7
        axs[i + 3].plot(P_m['entropy'][:chop_time])
    
    axs[i + 3].plot(entropies_time_mix_mean[week][:chop_time], ls='--', color='black', label='Mean')
    axs[i + 3].grid()
    axs[i + 3].set_title(f'Release span {week} weeks (Mixture)')
    axs[i + 3].set_xlim(0, chop_time)
    
axs[3].set_ylabel('Marginal Entropy (bits)')
axs[4].set_xlabel('Particle Age (days)')
axs[5].set_xlabel('Particle Age (days)')
axs[3].set_xlabel('Particle Age (days)')
axs[3].legend()

plt.savefig('../figs/FigS8-Combined-Representation_entropy_mixture.png', dpi=300)

# %%Print the min and max std of the mixtures 
for week in week_range:
    print(f"week {week} min std: {np.min(entropies_time_mix_std[week]):0.3f}, max std: {np.max(entropies_time_mix_std[week]):0.3f}")
    
for delta_r in [0.1, 1, 2]:
    print(f"delta_r {delta_r} min std: {np.min(entropies_space_mix_std[delta_r]):0.3f}, max std: {np.max(entropies_space_mix_std[delta_r]):0.3f}")

# %%
K_h = 1000
location = 'Cape_Hatteras'
member = 48  # member
path = f"/Volumes/Claudio SSD/Ensemble_article_data/analysis/prob_distribution/{location}_diffusion_long/P_diff_Kh_{K_h:01d}_m{member:03d}.nc"

P_m = xr.open_dataset(path)

K_h = 10
location = 'Cape_Hatteras'
member = 48  # member
path = f"/Volumes/Claudio SSD/Ensemble_article_data/analysis/prob_distribution/{location}_diffusion_long/P_diff_Kh_{K_h:01d}_m{member:03d}.nc"

P_m_10 = xr.open_dataset(path)

#%%
fig, ax1 = plt.subplots()

color1 = 'tab:blue'
ax1.set_xlabel('Time (days)')
ax1.set_ylabel('Entropy (bits)', color=color1)
l1, = ax1.plot(P_m['time'], P_m['entropy'], label=r'$K_h=1000 \ m^2s^{-1}$', color=color1)
l2, = ax1.plot(P_m_10['time'], P_m_10['entropy'], label=r'$K_h=10 \ m^2s^{-1}$', color=color1, ls='--')
ax1.tick_params(axis='y', labelcolor=color1)
# ax1.set_xlim(0, 2189)
ax1.semilogx()

ax2 = ax1.twinx()
color2 = 'tab:orange'
ax2.set_ylabel('Number of particles binned', color=color2)
ax2.plot(P_m['time'], P_m['number_particles_binned'], label='Number of particles binned', color=color2)
ax2.plot(P_m_10['time'], P_m_10['number_particles_binned'], label='Number of particles binned K_h=10', color=color2, ls='--')
ax2.tick_params(axis='y', labelcolor=color2)

ax3 = ax1.twinx()
color3 = 'tab:green'
ax3.spines['right'].set_position(('outward', 60))
ax3.set_ylabel('Entropy (bits)/ Number of particles', color=color3)
ax3.plot(P_m['time'], P_m['entropy']/P_m['number_particles_binned'], label='Entropy/Number of particles', color=color3)
ax3.plot(P_m_10['time'], P_m_10['entropy']/P_m_10['number_particles_binned'], label='Entropy/Number of particles K_h=10', color=color3, ls='--')
ax3.tick_params(axis='y', labelcolor=color3)


# # Combine legends
lines = [l1, l2]
labels = [line.get_label() for line in lines]
ax1.legend(lines, labels, loc='lower right')

plt.savefig('../figs/FigS10-diff1000vsdiff10Entropy_curves_semilogx.png', dpi=300)

# plt.tight_layout()
# plt.ylabel('Value')

# %%
K_h = 1000
location = 'Cape_Hatteras'
member = 48  # member
path = f"/Volumes/Claudio SSD/Ensemble_article_data/analysis/prob_distribution/{location}_diffusion_long/P_diff_Kh_{K_h:01d}_m{member:03d}_subsample_4000.nc"

P_m_sub = xr.open_dataset(path)

K_h = 10
location = 'Cape_Hatteras'
member = 48  # member
path = f"/Volumes/Claudio SSD/Ensemble_article_data/analysis/prob_distribution/{location}_diffusion_long/P_diff_Kh_{K_h:01d}_m{member:03d}_subsample_4000.nc"

P_m_10_sub = xr.open_dataset(path)

#%%
fig, ax1 = plt.subplots()

color1 = 'tab:blue'
color2 = 'tab:orange'
ax1.set_xlabel('Time (days)')
ax1.set_ylabel('Entropy (bits)')
l1, = ax1.plot(P_m['time'], P_m['entropy'], label=r'$K_h=1000 \ m^2s^{-1}$')
l2, = ax1.plot(P_m_10['time'], P_m_10['entropy'], label=r'$K_h=10 \ m^2s^{-1}$', ls='--')
ax1.plot(P_m_sub['time'], P_m_sub['entropy'], label=r'$K_h=1000 \ m^2s^{-1}$, subsampled 4000 particles')
ax1.plot(P_m_10_sub['time'], P_m_10_sub['entropy'], label=r'$K_h=10 \ m^2s^{-1}$, subsampled 4000 particles', ls='--')


# ax1.tick_params(axis='y', labelcolor=color1)
ax1.set_xlim(0, 2189)
# ax1.semilogx()

# ax2 = ax1.twinx()
# color2 = 'tab:orange'
# ax2.set_ylabel('Number of particles binned', color=color2)
# ax2.plot(P_m['time'], P_m['number_particles_binned'], label='Number of particles binned', color=color2)
# ax2.plot(P_m_10['time'], P_m_10['number_particles_binned'], label='Number of particles binned K_h=10', color=color2, ls='--')
# ax2.tick_params(axis='y', labelcolor=color2)

# ax3 = ax1.twinx()
# color3 = 'tab:green'
# ax3.spines['right'].set_position(('outward', 60))
# ax3.set_ylabel('Entropy (bits)/ Number of particles', color=color3)
# ax3.plot(P_m['time'], P_m['entropy']/P_m['number_particles_binned'], label='Entropy/Number of particles', color=color3)
# ax3.plot(P_m_10['time'], P_m_10['entropy']/P_m_10['number_particles_binned'], label='Entropy/Number of particles K_h=10', color=color3, ls='--')
# ax3.tick_params(axis='y', labelcolor=color3)


# # Combine legends
lines = [l1, l2]
labels = [line.get_label() for line in lines]
ax1.legend(loc='lower right')

plt.savefig('../figs/FigS10-diff1000vsdiff10Entropy_curves_subsample_4000.png', dpi=300)
# %%
