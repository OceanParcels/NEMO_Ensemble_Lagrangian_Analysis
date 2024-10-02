#%%
import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
from tqdm import tqdm
import pickle
import pandas as pd
import seaborn as sns
import cmocean.cm as cmo
import sys
sys.path.append('../functions')
import hexbin_functions as hexfunc


def cross_entropy(P, Q, axis=0):
    # Cross entropy H_p(q) = - sum_x p(x) log q(x)
    # P is the true distribution, Q is the estimated distribution.
    # Pdf = Pdf / np.nansum(Pdf)  # Normalize Pdf to sum to 1, ignoring NaNs
    # Replace zeros with a very small number to avoid log(0)
    P_safe = np.where(P > 0, P, np.finfo(float).eps)
    return -np.nansum(Q * np.log2(P_safe), axis=axis)

def marginal_entropy(P, axis=0):
    # Shannon entropy
    # Pdf = Pdf / np.nansum(Pdf)  # Normalize Pdf to sum to 1, ignoring NaNs
    # Replace zeros with a very small number to avoid log(0)
    P_safe = np.where(P > 0, P, np.finfo(float).eps)
    return -np.nansum(P_safe * np.log2(P_safe), axis=axis)

def KLDivergence(P, Q, axis=0):
    # Kullback-Leibler divergence D_KL(P||Q)
    # P is the true distribution, Q is the estimated distribution.
    return cross_entropy(P, Q, axis=axis) - marginal_entropy(Q, axis=axis)

location = 'Cape_Hatteras'
delta_r = 0.1
subset = 1

base_path = "/storage/shared/oceanparcels/output_data/data_Claudio/NEMO_Ensemble/analysis/prob_distribution/"

file_path_AX = base_path + f"{location}_all_long/P_dr{delta_r*100:03.0f}_all_s{subset:03d}.nc"
P_AX = xr.open_dataset(file_path_AX)

hex_grid = hexfunc.int_to_hex(P_AX.hexint.values)
hexbin_grid = hexfunc.hexGrid(hex_grid, h3_res=3)


#%% KLD
N_members = 50

time_length = 2189 - 7*20
time_range = np.arange(0, time_length)

KLD_ALL_mean = {}
KLD_ALL_std = {}

for delta_ref in [0.1, 1., 2., 4, 12, 20]:
    
    KLDivergence_mean = {}
    KLDivergence_std = {}
    
    for set in [0.1, 2., 4, 20]:
        print(f'Processing set P: {delta_ref}, Q:{set}')

        _KLD = np.zeros((N_members**2, time_length))

        l = 0
        for member in tqdm(range(1, N_members+1)):
            for subset_ref in range(1, N_members+1):
                # print(f'member {member} and subset {subset_ref}')

                if delta_ref in [0.1, 1., 2.]:
                    file_path_ref = base_path + f"{location}_all_long/P_dr{delta_ref*100:03.0f}_all_s{subset_ref:03d}.nc"
                elif delta_ref in [4, 12, 20]:
                    file_path_ref = base_path + f"{location}_all_long/P_W{delta_ref:02d}_all_s{subset_ref:03d}.nc"
                
                P_ref = xr.open_dataset(file_path_ref)
                P_ref = P_ref.sortby('hexint')

                if set in [0.1, 1., 2.]:
                    file_path_Q = base_path + f"{location}_spatial_long/P_dr{set*100:03.0f}_m{member:03d}.nc"
                elif set in [4, 12, 20]:
                    file_path_Q = base_path + f"{location}_temporal_long/P_W{set:01d}_m{member:03d}.nc"
                
                P_Q = xr.open_dataset(file_path_Q)
                P_Q = P_Q.sortby('hexint')

                _KLD[l, :] = KLDivergence(P_ref['probability'], P_Q['probability'], axis=0)[:time_length]
                l+=1
                
        avg_KLD = np.mean(_KLD, axis=0)
        std_KLD = np.std(_KLD, axis=0)
        
        KLDivergence_mean[set] = avg_KLD
        KLDivergence_std[set] = std_KLD

    KLD_ALL_mean[delta_ref] = KLDivergence_mean
    KLD_ALL_std[delta_ref] = KLDivergence_std

#%% # save KLDivergence_mean and KLDivergence_std to file
with open(f'/storage/shared/oceanparcels/output_data/data_Claudio/NEMO_Ensemble/analysis/KLD_time/KLD_time_ALL_mean.pkl', 'wb') as f:
    pickle.dump(KLD_ALL_mean, f)
        
with open(f'/storage/shared/oceanparcels/output_data/data_Claudio/NEMO_Ensemble/analysis/KLD_time/KLD_time_ALL_std.pkl', 'wb') as f:
    pickle.dump(KLD_ALL_std, f)
    

#%%
fig, axs = plt.subplots(2, 2, figsize=(10, 8), sharex=True, sharey=True)
# fig.subplots_adjust(hspace=0.4, wspace=0.4)

# keys = list(KLD_ALL_mean[0.1].keys())
labels = ['A', 'B', 'C', 'D', 'E', 'F']

axs = axs.ravel()

for i, key in enumerate(KLD_ALL_mean.keys()):
    
    lss = [(0, (1, 1)), '-.', '--', (0, (3, 1, 1, 1, 1, 1))]
    colors = ['mediumblue', 'teal', 'darkred', 'black']
    
    ax = axs[i]
    for k, set in enumerate([0.1, 2., 4, 20]):
        if set > 2:
            labelz = f'{set:01d} weeks'
        else:
            labelz = f'$\delta_r = {set}^o$'
        
        ax.fill_between(time_range, KLD_ALL_mean[key][set] - KLD_ALL_std[key][set], 
                        KLD_ALL_mean[key][set] + KLD_ALL_std[key][set], alpha=0.2, 
                        color=colors[k])
        ax.plot(time_range, KLD_ALL_mean[key][set], label=labelz, color=colors[k], 
                linestyle=lss[k], linewidth=2)
    
    if key > 2:
        labeltt = f'Mixture {key:01d} weeks'
    else:
        labeltt = f'Mixture $\delta_r = {key}^o$'
        
    ax.text(0.05, 0.9, f'$\mathbf{{{labels[i]}}}$  Ref.: {labeltt}', fontsize=14, transform=ax.transAxes)

    ax.semilogx()
    ax.legend()
    # ax.grid()
    ax.set_xlim(0, time_length)
    ax.set_ylim(0, 30)

axs[2].set_xlabel('Time (days)')
axs[3].set_xlabel('Time (days)')
axs[0].set_ylabel('KL Divergence (bits)')
axs[2].set_ylabel('KL Divergence (bits)')

plt.tight_layout()

plt.savefig('../figs/KLDivergence_subplots.png', dpi=300)
# %%
DF = {}

for i, key in enumerate(KLD_ALL_mean.keys()):
    _kld = np.zeros(4)
    
    for k, set in enumerate([0.1, 2., 4, 20]):
        
        _kld[k] = np.mean(KLD_ALL_mean[key][set])
        
    if key > 2:
        labeltt = f'Mix. {key:01d} weeks'
    else:
        labeltt = f'Mix. $\delta_r = {key}^o$'
        
    DF[labeltt] = _kld

DF = pd.DataFrame(DF)
DF.set_index([[r'$\delta_r = 0.1^o$', r'$\delta_r = 2^o$', 'Mix 4 weeks', 'Mix 20 weeks']], inplace=True)

# %% Kullback-Leibler divergence plot 
fig, ax = plt.subplots(figsize=(6, 5))

cmmap = "Greens_r"

sns.heatmap(DF, annot=True, fmt=".3f", cmap=cmmap, ax=ax, cbar=True)

# Rotate y tick labels 90 degrees
ax.set_yticklabels(ax.get_yticklabels(), rotation=90, ha='center', fontsize=9)
# Rotate x tick labels 15 degrees
ax.set_xticklabels(ax.get_xticklabels(), rotation=15, ha='center', fontsize=9)

# Add colorbar label
cbar = ax.collections[0].colorbar
cbar.set_label('Time Average KL Divergence, $D_{P_{Mix}}(P_i)$ (bits)')

plt.tight_layout()
plt.show()

fig.savefig("../figs/FigX_Time_Average_KLDivergence.png", dpi=300)

# %%
