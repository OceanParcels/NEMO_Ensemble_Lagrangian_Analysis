#%%% LOAD PACKAGES
import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
import cartopy
from tqdm import tqdm
import pickle

import sys
sys.path.append('../functions')
import hexbin_functions as hexfunc

#%%% LOAD DATA
Particles = np.logspace(1.5, 3, 20, dtype=int)
Particles = np.concatenate(([1], Particles))

location = 'Cape_Hatteras'
delta_r = 0.1
subset = 1

file_path_AX = f"/storage/shared/oceanparcels/output_data/data_Claudio/NEMO_Ensemble/analysis/prob_distribution/{location}_all/P_dr{delta_r*100:03.0f}_all_s{subset}.nc"
P_AX = xr.open_dataset(file_path_AX)


redundancy_grids = {}
Hrate_grids = {}
H_grids = {}
H_std_grids = {}

 #%% Load data for each resolution %%%%%
for h_res in [2,3,4]:
   
    with open(f'../data/hexgrid_no_coast_h{h_res}.pkl', 'rb') as f:
        hexbin_grid = pickle.load(f)
    
    hexbin_grid = hexfunc.hexGrid(hexbin_grid, h3_res=h_res)
    p_max = 1/hexbin_grid.n_hex
    H_max = hexbin_grid.n_hex*(-p_max*np.log2(p_max))
    
    
    entropy_Nparticles = np.zeros((len(Particles), len(P_AX['entropy'].values)))
    redundancy = np.zeros(len(Particles))
    entropy_rate = np.zeros(len(Particles))
    entropy_rate_std = np.zeros(len(Particles))
    
    for i, p in enumerate(Particles):
        N = p*50
    
        file_path_AX = f"/storage/shared/oceanparcels/output_data/data_Claudio/NEMO_Ensemble/analysis/prob_distribution/{location}_all_Nparticles_h{h_res}/P_all_p{p:04d}_h{h_res}.nc"

        P_AX = xr.open_dataset(file_path_AX)
        
        entropy_Nparticles[i, :] = P_AX['entropy'].values
        
        # redundancy[i] = np.mean(P_AX['entropy'][600:].values)/H_max
        entropy_rate[i] = np.mean(P_AX['entropy'][600:].values)
        entropy_rate_std[i] = np.std(P_AX['entropy'][600:].values)
    
    redundancy_grids[h_res] = redundancy
    Hrate_grids[h_res] = entropy_rate
    H_grids[h_res] = entropy_Nparticles
    H_std_grids[h_res] = entropy_rate_std
    

#%%%%%%%% Plotting %%%%%%%%%%%
fig, ax1 = plt.subplots()

ax1.plot(Particles, Hrate_grids[2], 'ob-',  label='h = 2')
ax1.fill_between(Particles, Hrate_grids[2] - H_std_grids[2], 
                Hrate_grids[2] + H_std_grids[2], color='b', alpha=0.2)
ax1.plot(Particles, Hrate_grids[3], 'sk--',  label='h = 3')
ax1.fill_between(Particles, Hrate_grids[3] - H_std_grids[3], 
                Hrate_grids[3] + H_std_grids[3], color='k', alpha=0.2)
ax1.plot(Particles, Hrate_grids[4], 'xr-.',  label='h = 4')
ax1.fill_between(Particles, Hrate_grids[4] - H_std_grids[4], 
                Hrate_grids[4] + H_std_grids[4], color='r', alpha=0.2) 
# ax1.loglog()
ax1.set_xlabel('Particles sub-sampled per member')
ax1.set_ylabel('Entropy (bits)')
ax1.tick_params(axis='y')
ax1.legend()
ax1.grid()
plt.savefig(f'../figs/Fig2_entropy_hres.png', dpi=300)
# %%
