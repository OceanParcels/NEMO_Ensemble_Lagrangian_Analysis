#%%
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import cartopy 

#%%

# Define the location, member, and delta_r values
location = 'Cape_Hatteras'
member = 3
delta_r = 1.

# Define the file path for the spatial data
file_path = f"/storage/shared/oceanparcels/output_data/data_Claudio/NEMO_Ensemble/{location}/spatial_long/dr_{delta_r*100:03.0f}/{location}_dr{delta_r*100:03.0f}_m{member:03d}.zarr"

# Open and compute the spatial dataset
pset_space = xr.open_zarr(file_path)
pset_space.compute()

# Define the number of weeks
week = 20

# Define the file path for the temporal data
path = f"/storage/shared/oceanparcels/output_data/data_Claudio/NEMO_Ensemble/{location}/temporal_long/W_{week:01d}/{location}_W{week:01d}_m{member:03d}.zarr"

# Open and compute the temporal dataset
pset_temp = xr.open_zarr(path)
pset_temp.compute()

#%% Mixture set of particles
# Define the number of particles
N_particles = 50
delta_r = 0.1

# Initialize arrays to store particle positions
mix_lons = np.zeros((N_particles, len(pset_space.obs)))
mix_lats = np.zeros((N_particles, len(pset_space.obs)))

# Loop over each particle
for l, member in enumerate(range(1, N_particles+1)):
        print(f"Member {member}")
        # Define the file path for the particle data
        file_path = f"/storage/shared/oceanparcels/output_data/data_Claudio/NEMO_Ensemble/{location}/spatial_long/dr_{delta_r*100:03.0f}/{location}_dr{delta_r*100:03.0f}_m{member:03d}.zarr"
        
        # Open the particle dataset
        pset = xr.open_zarr(file_path)
        
        # Store the particle positions in the arrays
        mix_lons[l, :] = pset.lon[0,:].values
        mix_lats[l, :] = pset.lat[0,:].values


#%% Plot NA_domain on a map
np.random.seed(38)
depth = 0

indexes = np.random.randint(0, 7000, N_particles)

indexes_space = np.arange(1000, 1000 + 40*4) # np.concatenate([np.arange(1, 41), np.arange(40*2+1, 40*3+1)])
indexes_space = indexes_space[::2]

indexes_space = np.random.randint(1000, 4000, N_particles)
indexes_space_points = np.random.randint(0, 7000, 50)


#%%%%%%%%%%%%% PLOT %%%%%%%%%%%%%%
# Create a figure and axes with PlateCarree projection
t = 14

fig = plt.figure()
ax = plt.axes(projection=cartopy.crs.PlateCarree())

# Set the extent of the map
ax.set_extent([-77, -68, 31.5, 39.5], crs=cartopy.crs.PlateCarree())

# Add land feature to the map
ax.add_feature(cartopy.feature.LAND, zorder=0, edgecolor='black')

# Add gridlines
gl = ax.gridlines(draw_labels=True, zorder=0, linestyle='--', linewidth=0.4)
gl.xlabels_top = False
gl.ylabels_right = False

# Scatter plot for varying space particles
ax.scatter(pset_space.lon[indexes_space, 0], pset_space.lat[indexes_space, 0],
                   s=30, color='blueviolet', alpha=1, label='Varying Space', 
                   zorder=10, edgecolor='black')

# Scatter plot for varying time particle
ax.scatter(pset_temp.lon[0, 0], pset_temp.lat[0, 0],
                   s=50, color='orangered', alpha=1, marker='s', 
                   label='Varying Time', zorder=10, edgecolor='black')

# Plot trajectories for varying space particles
ax.plot(pset_space.lon[indexes_space, :t].T, pset_space.lat[indexes_space, :t].T, c='blueviolet', 
                ls='-', alpha=0.5)

# Plot trajectories for varying time particle
ax.plot(pset_temp.lon[indexes, :t].T, pset_temp.lat[indexes, :t].T, c='orangered',
                ls='-', alpha=0.5)

# Plot trajectories for varying members particles
ax.plot(mix_lons[:, :t].T, mix_lats[:, :t].T, c='k',
                ls='-', alpha=0.5, label='Varying Members', zorder=0)

for i in range(N_particles):
        ax.annotate('', xy=(mix_lons[i, t], mix_lats[i, t]), 
                xytext=(mix_lons[i, t-1], mix_lats[i, t-1]),
                arrowprops=dict(arrowstyle="-|>", color='black', lw=1.5, alpha=0.5), 
                zorder=10) 

for i in indexes:
        ax.annotate('', xy=(pset_temp.lon[i, t], pset_temp.lat[i, t]), 
                xytext=(pset_temp.lon[i, t-1], pset_temp.lat[i, t-1]),
                arrowprops=dict(arrowstyle="-|>", color='orangered', lw=1.5, alpha=0.5), 
                zorder=10)
        
        
for i in indexes_space:
        ax.annotate('', xy=(pset_space.lon[i, t], pset_space.lat[i, t]), 
                xytext=(pset_space.lon[i, t-1], pset_space.lat[i, t-1]),
                arrowprops=dict(arrowstyle="-|>", color='blueviolet', lw=1.5, alpha=0.5), 
                zorder=10)

# Add legend
handles, labels = ax.get_legend_handles_labels()
handles = [handles[0], handles[1], handles[-1]]
labels = [labels[0], labels[1], labels[-1]]
ax.legend(handles, labels, shadow=True, fontsize='small')

# Save the figure
plt.savefig(f'../figs/Fig1_schematic.png', dpi=300)

#%%