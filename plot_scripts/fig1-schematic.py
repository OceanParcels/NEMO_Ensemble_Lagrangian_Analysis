#%%
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import cartopy 


#%%

# Define the location, member, and delta_r values
location = 'Cape_Hatteras'
member = 3
delta_r = 0.5

# Define the file path for the spatial data
file_path = f"/storage/shared/oceanparcels/output_data/data_Claudio/NEMO_Ensemble/{location}/spatial/dr_{delta_r*100:03.0f}/{location}_dr{delta_r*100:03.0f}_m{member:03d}.zarr"

# Open and compute the spatial dataset
pset_space = xr.open_zarr(file_path)
pset_space.compute()

# Define the number of weeks
week = 4

# Define the file path for the temporal data
path = f"/storage/shared/oceanparcels/output_data/data_Claudio/NEMO_Ensemble/{location}/temporal/W_{week:01d}/{location}_W{week:01d}_m{member:03d}.zarr"

# Open and compute the temporal dataset
pset_temp = xr.open_zarr(path)
pset_temp.compute()

#%% Mixture set of particles
# Define the number of particles
N_particles = 50

# Initialize arrays to store particle positions
mix_lons = np.zeros((N_particles, len(pset_space.obs)))
mix_lats = np.zeros((N_particles, len(pset_space.obs)))

# Loop over each particle
for l, member in enumerate(range(1, N_particles+1)):
        # Define the file path for the particle data
        file_path = f"/storage/shared/oceanparcels/output_data/data_Claudio/NEMO_Ensemble/{location}/spatial/dr_{delta_r*100:03.0f}/{location}_dr{delta_r*100:03.0f}_m{member:03d}.zarr"
        
        # Open the particle dataset
        pset = xr.open_zarr(file_path)
        
        # Store the particle positions in the arrays
        mix_lons[l, :] = pset.lon[0,:].values
        mix_lats[l, :] = pset.lat[0,:].values


#%% Plot NA_domain on a map
depth = 0
t = 30
indexes = np.random.randint(0, 1000, N_particles)

indexes_space = np.concatenate([np.arange(1, 41), np.arange(40*6+1, 40*7+1)])
indexes_space = indexes_space[::2]


#%%%%%%%%%%%%% PLOT %%%%%%%%%%%%%%
# Create a figure and axes with PlateCarree projection
fig = plt.figure()
ax = plt.axes(projection=cartopy.crs.PlateCarree())

# Set the extent of the map
ax.set_extent([-77, -68, 31, 39], crs=cartopy.crs.PlateCarree())

# Add land feature to the map
ax.add_feature(cartopy.feature.LAND, zorder=0, edgecolor='black')

# Add gridlines
gl = ax.gridlines(draw_labels=True, zorder=0, linestyle='--', linewidth=0.4)
gl.xlabels_top = False
gl.ylabels_right = False

# Scatter plot for varying space particles
ax.scatter(pset_space.lon[indexes_space, 0], pset_space.lat[indexes_space, 0],
                   s=10, color='blueviolet', alpha=1, label='Varying Space', zorder=10)

# Scatter plot for varying time particle
ax.scatter(pset_temp.lon[0, 0], pset_temp.lat[0, 0],
                   s=20, color='orangered', alpha=1, marker='s', label='Varying Time', zorder=10)

# Plot trajectories for varying space particles
ax.plot(pset_space.lon[indexes_space, :t].T, pset_space.lat[indexes_space, :t].T, c='blueviolet', 
                ls='-', alpha=0.5)

# Plot trajectories for varying time particle
ax.plot(pset_temp.lon[indexes, :t].T, pset_temp.lat[indexes, :t].T, c='orangered',
                ls='-', alpha=0.5)

# Plot trajectories for varying members particles
ax.plot(mix_lons[:, :t].T, mix_lats[:, :t].T, c='k',
                ls='-', alpha=0.5, label='Varying Members', zorder=0)

# ax.annotate('', xy=(mix_lons[t, 0], mix_lats[t, 0]), xytext=(mix_lons[t-1, 0], mix_lats[t-1,0]),
#             arrowprops=dict(arrowstyle="->", color='green', lw=1.5))


# Add legend
handles, labels = ax.get_legend_handles_labels()
handles = [handles[0], handles[1], handles[-1]]
labels = [labels[0], labels[1], labels[-1]]
ax.legend(handles, labels, shadow=True, fontsize='small')

# Save the figure
plt.savefig(f'../figs/Fig1_schematic.png', dpi=300)
    
#%%