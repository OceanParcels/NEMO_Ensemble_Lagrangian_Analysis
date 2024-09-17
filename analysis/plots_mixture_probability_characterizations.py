import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
import pickle
import sys
sys.path.append('../functions')
import hexbin_functions as hexfunc

Particles = np.logspace(1.5, 3, 20, dtype=int)
Particles = np.array(list(Particles)) # + [25, 30, 35, 100, 500])

redundancy_grids = {}
Hrate_grids = {}
H_grids = {}
H_std_grids = {}


for h_res in [2,3]:
    #%% Load the hexbin_grid for the domain
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
        
        redundancy[i] = np.mean(P_AX['entropy'][600:].values)/H_max
        entropy_rate[i] = np.mean(P_AX['entropy'][600:].values)
        entropy_rate_std[i] = np.std(P_AX['entropy'][600:].values)
    
    redundancy_grids[h_res] = redundancy
    Hrate_grids[h_res] = entropy_rate
    H_grids[h_res] = entropy_Nparticles
    H_std_grids[h_res] = entropy_rate_std