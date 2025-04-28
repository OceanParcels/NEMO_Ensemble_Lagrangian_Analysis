#%% 
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import cartopy 
#%%


#Some SIMULATION parameter
location = 'Cape_Hatteras'
start_time = np.datetime64('2010-01-02')
end_time =  np.datetime64('2015-12-31')
#% 
K_h = 10
# member = 1

for member in range(1, 50):
    outfile = f"/storage/shared/oceanparcels/output_data/data_Claudio/NEMO_Ensemble/{location}/diff_long/diff_Kh_{K_h:01d}/{location}_diff_Kh_{K_h:01d}_m{member:03d}.zarr"
    # outfile = f"/storage/shared/oceanparcels/output_data/data_Claudio/NEMO_Ensemble/{location}/diff_long/diff_Kh_{K_h}/{location}_diff_Kh_{K_h}_m{member:03d}.zarr"

    # print("Output file: ", outfile)
    mem_x = xr.load_dataset(outfile);
    
    # Are particles been deleted?
    # print(mem_x.lon.shape)
    # print(mem_x.lon[:, -12].isnull().sum()) #

    # which particles are being deleted?
    deleted_particles = np.where(mem_x.lon[:, -12].isnull())[0]
    print("Member: ", member, "particles deleted", deleted_particles)
    
    if len(deleted_particles) > 0:
        
        fig = plt.figure()
        ax = plt.axes(projection=cartopy.crs.PlateCarree())
        ax.add_feature(cartopy.feature.LAND, zorder=0, edgecolor='black')
        # Add gridlines
        gl = ax.gridlines(draw_labels=True, zorder=0, linestyle='--', linewidth=0.4)
        gl.xlabels_top = False
        gl.ylabels_right = False
        dt = -12
        ax.set_title(f"Member {member}")
        for i in deleted_particles:
            ax.plot(mem_x.lon[i, :dt], mem_x.lat[i, :dt])

# %% plot lat and lon of particles
plt.figure()
for i in range(0, 7500):
    plt.plot(mem_x.lon[i, :], -mem_x.z[i, :])

#%% 

fig = plt.figure()
ax = plt.axes(projection=cartopy.crs.PlateCarree())

# Set the extent of the map
# ax.set_extent([-77, -68, 31.5, 39.5], crs=cartopy.crs.PlateCarree())

# Add land feature to the map
ax.add_feature(cartopy.feature.LAND, zorder=0, edgecolor='black')

# Add gridlines
gl = ax.gridlines(draw_labels=True, zorder=0, linestyle='--', linewidth=0.4)
gl.xlabels_top = False
gl.ylabels_right = False

dt = -12
for i in range(0, 7500):
    ax.plot(mem_x.lon[i, :dt], mem_x.lat[i, :dt])

# %% Are particles been deleted?
print(mem_x.lon.shape)
print(mem_x.lon[:, -12].isnull().sum()) #

# which particles are being deleted?
deleted_particles = np.where(mem_x.lon[:, -12].isnull())[0]
print(deleted_particles)

# %%
for i in deleted_particles:
    plt.plot(mem_x.lon[i, :dt], -mem_x.z[i, :dt])

#%% 

fig = plt.figure()
ax = plt.axes(projection=cartopy.crs.PlateCarree())

# Set the extent of the map
# ax.set_extent([-77, -68, 31.5, 39.5], crs=cartopy.crs.PlateCarree())

# Add land feature to the map
ax.add_feature(cartopy.feature.LAND, zorder=0, edgecolor='black')

# Add gridlines
gl = ax.gridlines(draw_labels=True, zorder=0, linestyle='--', linewidth=0.4)
gl.xlabels_top = False
gl.ylabels_right = False

dt = -12
for i in deleted_particles:
    ax.plot(mem_x.lon[i, :dt], mem_x.lat[i, :dt])

# %%
