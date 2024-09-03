#%% Load the packages
import numpy as np
import matplotlib.pyplot as plt
import xarray as xr
import cartopy
from tqdm import tqdm
# import pickle

import sys
sys.path.append('../functions')
import hexbin_functions as hexfunc

#%% Load the data

location = 'Cape_Hatteras'
depth = 1
member = 1

file_path = f"/storage/shared/oceanparcels/output_data/data_Claudio/NEMO_Ensemble/{location}/temporal_connectivity/dep_{depth:01d}/{location}_dep{depth:01d}_m{member:03d}.zarr"
pset = xr.open_zarr(file_path)
pset = pset.compute()


# %%
lats = pset.lat.values

p_index, t_index = np.where(lats[:,:]>53)
p_index = np.unique(p_index)

print(p_index)

# %%
fig = plt.figure()
ax = plt.axes(projection=cartopy.crs.PlateCarree())
# ax.set_extent([-80, -10, 25, 55], crs=cartopy.crs.PlateCarree())
# ax.set_extent([-80, -60, 28, 37], crs=cartopy.crs.PlateCarree())
ax.coastlines()
ax.gridlines(draw_labels=True, zorder=0, linestyle='--', linewidth=0.5)

for i in p_index:
    # ax.scatter(pset.lon[i,0], pset.lat[i,0], s=1)
    ax.plot(pset.lon[i,:], pset.lat[i,:])
# %%
