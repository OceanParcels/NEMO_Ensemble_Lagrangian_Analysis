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
file_path = f"/Volumes/Claudio SSD/Ensemble_article_data/simulations/dr_{delta_r*100:03.0f}/{location}_dr{delta_r*100:03.0f}_m{member:03d}.nc"

# Open and compute the spatial dataset
pset_space = xr.load_dataset(file_path)
pset_space.compute()

# Define the number of weeks
week = 20

# Define the file path for the temporal data
path = f"/Volumes/Claudio SSD/Ensemble_article_data/simulations/W_{week:01d}/{location}_W{week:01d}_m{member:03d}.nc"

# Open and compute the temporal dataset
pset_temp = xr.load_dataset(path)
pset_temp.compute()

#%% Define the file path for the diffusion data
member = 3
path = f"/Volumes/Claudio SSD/Ensemble_article_data/simulations/diff_Kh_10/Cape_Hatteras_diff_Kh_10_m{member:03d}.nc"
print(path)
# Open and compute the temporal dataset
pset_diff = xr.load_dataset(path)
pset_diff.compute()

member = 3
path = f"/Volumes/Claudio SSD/Ensemble_article_data/simulations/diff_Kh_1000/Cape_Hatteras_diff_Kh_1000_m{member:03d}.zarr"
print(path)
# Open and compute the temporal dataset
pset_diff2 = xr.open_zarr(path)
pset_diff2.compute()


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
        file_path = f"/Volumes/Claudio SSD/Ensemble_article_data/simulations/dr_{delta_r*100:03.0f}/{location}_dr{delta_r*100:03.0f}_m{member:03d}.nc"
        
        # Open the particle dataset
        pset = xr.load_dataset(file_path)
        
        # Store the particle positions in the arrays
        mix_lons[l, :] = pset.lon[0,:].values
        mix_lats[l, :] = pset.lat[0,:].values


#%% Plot NA_domain on a map
np.random.seed(38)
depth = 0

indexes = np.random.randint(0, 7000, N_particles)

indexes_space = np.arange(1000, 1000 + 40*4) # np.concatenate([np.arange(1, 41), np.arange(40*2+1, 40*3+1)])
indexes_space = indexes_space[::2]

indexes_space = np.random.randint(1000, 4000, 150)



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
                   zorder=12, edgecolor='black')

# Scatter plot for varying time particle
ax.scatter(pset_temp.lon[0, 0], pset_temp.lat[0, 0],
                   s=50, color='gold', alpha=1, marker='s', 
                   label='Varying Time', zorder=12, edgecolor='black')

# Plot trajectories for varying space particles
ax.plot(pset_space.lon[indexes_space, :t].T, pset_space.lat[indexes_space, :t].T, c='blueviolet', 
                ls='-', alpha=0.5, zorder=10)

# Plot trajectories for varying time particle
ax.plot(pset_temp.lon[indexes, :t].T, pset_temp.lat[indexes, :t].T, c='orangered',
                ls='-', alpha=0.5, zorder=11)

# Plot trajectories for varying members particles
ax.plot(mix_lons[:, :t].T, mix_lats[:, :t].T, c='k',
                ls='-', alpha=0.5, label='Varying Members', zorder=9)

for i in range(N_particles):
        ax.annotate('', xy=(mix_lons[i, t], mix_lats[i, t]), 
                xytext=(mix_lons[i, t-1], mix_lats[i, t-1]),
                arrowprops=dict(arrowstyle="-|>", color='black', lw=1.5, alpha=0.5), 
                zorder=9) 

for i in indexes:
        ax.annotate('', xy=(pset_temp.lon[i, t], pset_temp.lat[i, t]), 
                xytext=(pset_temp.lon[i, t-1], pset_temp.lat[i, t-1]),
                arrowprops=dict(arrowstyle="-|>", color='orangered', lw=1.5, alpha=0.5), 
                zorder=11)
        
for i in indexes_space:
        ax.annotate('', xy=(pset_space.lon[i, t], pset_space.lat[i, t]), 
                xytext=(pset_space.lon[i, t-1], pset_space.lat[i, t-1]),
                arrowprops=dict(arrowstyle="-|>", color='blueviolet', lw=1.5, alpha=0.5), 
                zorder=10)

# Add legend
handles, labels = ax.get_legend_handles_labels()
handles = [handles[-1], handles[0], handles[1]]
labels = [labels[-1], labels[0], labels[1]]
ax.legend(handles, labels, shadow=True, fontsize='small')

# Save the figure
plt.savefig(f'../figs/Fig1_schematic.png', dpi=300)

#%% plot one trajectory 

plt.scatter(0, 0, s=30, color='blueviolet', alpha=1, 
                   zorder=10, edgecolor='black')
plt.plot([0, 0.005], [0,0], c='blueviolet', 
                ls='-', alpha=1)
plt.annotate('', xy=(0.01, 0), 
                xytext=(0.005, 0),
                arrowprops=dict(arrowstyle="-|>", color='blueviolet', lw=1.5, alpha=1))

plt.scatter(0, 0.025, s=40, color='gold', alpha=1, marker='s',
                   zorder=10, edgecolor='black')
plt.plot([0, 0.005], [0.025,0.025], c='orangered', 
                ls='-', alpha=1)
plt.annotate('', xy=(0.01, 0.025), 
                xytext=(0.005, 0.025),
                arrowprops=dict(arrowstyle="-|>", color='orangered', lw=1.5, alpha=1))

plt.scatter(0, 0.05, s=40, color='gold', alpha=1, marker='s',
                   zorder=10, edgecolor='black')
plt.plot([0, 0.005], [0.05,0.05], c='black', 
                ls='-', alpha=1)
plt.annotate('', xy=(0.01, 0.05), 
                xytext=(0.005, 0.05),
                arrowprops=dict(arrowstyle="-|>", color='black', lw=1.5, alpha=1))

plt.scatter(0, 0.075, s=40, color='gold', alpha=1, marker='s',
                         zorder=10, edgecolor='black')
plt.plot([0, 0.005], [0.075,0.075], c='green',
                ls='-', alpha=1)
plt.annotate('', xy=(0.01, 0.075),
                xytext=(0.005, 0.075),
                arrowprops=dict(arrowstyle="-|>", color='green', lw=1.5, alpha=1))

plt.xlim(-0.1, 0.1)
plt.ylim(-0.1, 0.1)

plt.savefig(f'../figs/patch_for_Fig1_schematic.png', dpi=300)

# %% plot 2x2 subplots ploting varyin members, varying space, varying time and diffusion in idivudual subplots

t = 30
fig, axs = plt.subplots(2, 2, figsize=(7, 7), subplot_kw={'projection': cartopy.crs.PlateCarree()},
                          gridspec_kw={'wspace': 0.05, 'hspace': -0.15})
# Set the extent of the map
for ax in axs.flat:
    ax.set_extent([-77, -68, 31.5, 39.5], crs=cartopy.crs.PlateCarree())
    ax.add_feature(cartopy.feature.LAND, zorder=0, edgecolor='black')
    gl = ax.gridlines(draw_labels=True, zorder=0, linestyle='--', linewidth=0.4)
    gl.right_labels = False
    gl.top_labels = False

    if ax == axs[0, 0]:
        gl.bottom_labels = False
    elif ax == axs[0, 1]:
        gl.left_labels = False
        gl.bottom_labels = False
    elif ax == axs[1, 1]:
        gl.left_labels = False

# Plot varying members  
axs[0, 0].scatter(mix_lons[:, 0], mix_lats[:, 0],
                   s=50, color='gold', alpha=1, label='Varying Members', 
                   zorder=12, edgecolor='black', marker='s')

axs[0, 0].plot(mix_lons[:, :t].T, mix_lats[:, :t].T, c='k',
                ls='-', alpha=0.5, zorder=10)
for i in range(N_particles):
    axs[0, 0].annotate('', xy=(mix_lons[i, t], mix_lats[i, t]), 
                xytext=(mix_lons[i, t-1], mix_lats[i, t-1]),
                arrowprops=dict(arrowstyle="-|>", color='black', lw=1.5, alpha=0.5), 
                zorder=9)

# Plot varying space


axs[0, 1].scatter(pset_space.lon[indexes_space, 0], pset_space.lat[indexes_space, 0],
                   s=30, color='blueviolet', alpha=1, label='Varying Space', 
                   zorder=12, edgecolor='black')

axs[0, 1].plot(pset_space.lon[indexes_space[::2], :t].T, pset_space.lat[indexes_space[::2], :t].T, c='blueviolet',
                ls='-', alpha=0.5, zorder=10)
for i in indexes_space[::2]:
        axs[0, 1].annotate('', xy=(pset_space.lon[i, t], pset_space.lat[i, t]), 
                        xytext=(pset_space.lon[i, t-1], pset_space.lat[i, t-1]),
                        arrowprops=dict(arrowstyle="-|>", color='blueviolet', lw=1.5, alpha=0.5), 
                        zorder=10)
# Plot varying time

axs[1, 0].scatter(pset_temp.lon[0, 0], pset_temp.lat[0, 0],
                   s=50, color='gold', alpha=1, marker='s', 
                   label='Varying Time', zorder=12, edgecolor='black')
axs[1, 0].plot(pset_temp.lon[indexes, :t].T, pset_temp.lat[indexes, :t].T, c='orangered',
                ls='-', alpha=0.5, zorder=11)
for i in indexes:
        axs[1, 0].annotate('', xy=(pset_temp.lon[i, t], pset_temp.lat[i, t]), 
                        xytext=(pset_temp.lon[i, t-1], pset_temp.lat[i, t-1]),
                        arrowprops=dict(arrowstyle="-|>", color='orangered', lw=1.5, alpha=0.5), 
                        zorder=11)
# Plot diffusion
t = 30
axs[1, 1].plot(pset_diff2.lon[indexes, :t].T, pset_diff2.lat[indexes, :t].T, c='green',
                ls='-', alpha=0.5, zorder=11)
for i in indexes:
        axs[1, 1].annotate('', xy=(pset_diff2.lon[i, t], pset_diff2.lat[i, t]), 
                        xytext=(pset_diff2.lon[i, t-1], pset_diff2.lat[i, t-1]),
                        arrowprops=dict(arrowstyle="-|>", color='darkgreen', lw=1.5, alpha=0.5), 
                        zorder=11)


t = 50
axs[1, 1].scatter(pset_diff.lon[0, 0], pset_diff.lat[0, 0],
                   s=50, color='gold', alpha=1, marker='s', 
                   label='Added Diffusion', zorder=25, edgecolor='black')

axs[1, 1].plot(pset_diff.lon[indexes[:15], :t].T, pset_diff.lat[indexes[:15], :t].T, c='limegreen',
                ls='-', alpha=1, zorder=20)
for i in indexes[:15]:
        axs[1, 1].annotate('', xy=(pset_diff.lon[i, t], pset_diff.lat[i, t]), 
                        xytext=(pset_diff.lon[i, t-1], pset_diff.lat[i, t-1]),
                        arrowprops=dict(arrowstyle="-|>", color='limegreen', lw=1.5, alpha=0.5), 
                        zorder=25)


# Add legend
axs[0, 0].legend(loc='upper left', shadow=True, fontsize='small')
axs[0, 1].legend(loc='upper left', shadow=True, fontsize='small')
axs[1, 0].legend(loc='upper left', shadow=True, fontsize='small')
axs[1, 1].legend(loc='upper left', shadow=True, fontsize='small')

# add A, B, C, D labels to the subplots lower right corner, fontsize 14
axs[0, 0].text(-68.8, 31.8, 'A', fontsize=12, fontweight='bold', transform=cartopy.crs.PlateCarree())
axs[0, 1].text(-68.8, 31.8, 'B', fontsize=12, fontweight='bold', transform=cartopy.crs.PlateCarree())
axs[1, 0].text(-68.8, 31.8, 'C', fontsize=12, fontweight='bold', transform=cartopy.crs.PlateCarree())
axs[1, 1].text(-68.8, 31.8, 'D', fontsize=12, fontweight='bold', transform=cartopy.crs.PlateCarree())

# Save the figure
plt.savefig(f'../figs/Fig1_schematic_2x2.png', dpi=300)

# %%
