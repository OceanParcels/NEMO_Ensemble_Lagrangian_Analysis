#%%
from parcels import FieldSet, ParticleSet, JITParticle, Variable, ParticleFile
import numpy as np
import xarray as xr
from glob import glob
from datetime import timedelta as delta
import matplotlib.pyplot as plt
import h3
import cartopy

import pickle
import sys
sys.path.append('../functions')
import hexbin_functions as hexfunc
from scipy.optimize import curve_fit

mask_file = '/storage/shared/oceanparcels/input_data/NEMO_Ensemble/GRID/NATL025-CJMenobs01_byte_mask.nc'
mask = xr.open_dataset(mask_file, decode_times=False)

with open('../data/hexgrid_no_coast.pkl', 'rb') as f:
    hexagons_grid = pickle.load(f)
    
grid = hexfunc.hexGrid(hexagons_grid)
#%%
class EnsembleParticle(JITParticle):
    """
    Particle class definition with additional variables
    """
    # dynamic variables
    u = Variable('u', dtype=np.float32, initial=0)
    v = Variable('v', dtype=np.float32, initial=0)
    w = Variable('w', dtype=np.float32, initial=0)
    
    hexbin_id = Variable('hexbin_id', dtype=np.int16, initial=0)

def SampleField(particle, fieldset, time):
    """
    Sample the fieldset at the particle location and store it in the
    particle variable.
    """
    (ui, vi, wi) = fieldset.UVW.eval(time, particle.depth, particle.lat, particle.lon, 
                                     particle=particle, applyConversion=False)
    particle.u = ui
    particle.v = vi
    particle.w = wi


def compute_cosine_similarity(vec1, vec2):
    """Compute cosine similarity between two 2D vectors."""
    dot_product = np.dot(vec1, vec2)
    magnitude1 = np.linalg.norm(vec1)
    magnitude2 = np.linalg.norm(vec2)
    return dot_product / (magnitude1 * magnitude2)


def compute_correlation_function(x_components, y_components, max_lag):
    """
    Compute the correlation function based on lag for a set of 2D vectors.

    Parameters:
    - x_components (list): List of x-components of the 2D vectors.
    - y_components (list): List of y-components of the 2D vectors.
    - max_lag (int): Maximum lag to consider for computing the correlation function.

    Returns:
    - correlations (list): List of tuples containing the lag and the average correlation for that lag.

    The function computes the correlation function based on lag for a set of 2D vectors.
    It calculates the correlation between each pair of vectors separated by a lag value
    ranging from 1 to max_lag. The correlation is computed using the Pearson correlation
    coefficient. The average correlation for each lag is then calculated and stored in a list
    of tuples, where each tuple contains the lag and the average correlation for that lag.
    The list of tuples is returned as the output of the function.
    """
    n = len(x_components)
    correlations = []

    for lag in range(1, max_lag + 1):
        lag_correlations = []

        for i in range(n - lag):
            vec1 = (x_components[i], y_components[i])
            vec2 = (x_components[i + lag], y_components[i + lag])
            correlation = compute_cosine_similarity(vec1, vec2)
            lag_correlations.append(correlation)

        average_correlation = np.mean(lag_correlations)
        correlations.append(average_correlation)

    return correlations

#%% PArcels

start_time = np.datetime64('2010-01-02')

loc1_lon = -74.0
loc1_lat = 35.5

# Find the hexagon containing the location
loc1_hex = h3.geo_to_h3(loc1_lat, loc1_lon, 3)
loc1_lat, loc1_lon = h3.h3_to_geo(loc1_hex)
lon_0 = loc1_lon
lat_0 = loc1_lat

# Define the rings where we place the particles
span = 1.75
L_range = np.linspace(-span, span, 25)
theta_range = np.arange(0, 2*np.pi, np.pi/20)
lonp = [lon_0]
latp = [lat_0]

for r in L_range:
    for theta in theta_range:
    
        lonp.append(lon_0 + np.sin(theta)*r) 
        latp.append(lat_0 + np.cos(theta)*r)
        
times = [start_time]*len(lonp)
        
print("N particles: ", len(lonp))

loc1_lon = -74.0
loc1_lat = 35.5

# Find the hexagon containing the location
loc1_hex = h3.geo_to_h3(loc1_lat, loc1_lon, 3)
loc1_lat, loc1_lon = h3.h3_to_geo(loc1_hex)
lon_0 = loc1_lon
lat_0 = loc1_lat

# Define the rings where we place the particles
span = 1.6
L_range = np.arange(-span, span, 0.01)
theta_range = np.arange(0, 2*np.pi, np.pi/2)
lonp = [lon_0]
latp = [lat_0]

max_lag = 200
R_members = np.zeros((50, max_lag))
theta = np.pi/2

for member in range(1, 51):
    print(f"Member {member}")

    data_path = '/storage/shared/oceanparcels/input_data/NEMO_Ensemble/'
    ufiles = sorted(glob(f"{data_path}NATL025-CJMCYC3.{member:03d}-S/1d/2010/NATL025*U.nc"))
    vfiles = [f.replace('U.nc', 'V.nc') for f in ufiles]
    wfiles = [f.replace('U.nc', 'W.nc') for f in ufiles]
    mesh_mask = f"{data_path}GRID/coordinates_NATL025_v2.nc"
    maskfile = f"{data_path}GRID/NATL025-CJMenobs01_byte_mask.nc"

    filenames = {'U': {'lon': mesh_mask, 'lat': mesh_mask, 'depth': wfiles[0], 'data': ufiles},
                'V': {'lon': mesh_mask, 'lat': mesh_mask, 'depth': wfiles[0], 'data': vfiles},
                'W': {'lon': mesh_mask, 'lat': mesh_mask, 'depth': wfiles[0], 'data': wfiles},
                'mask': {'lon': mesh_mask, 'lat': mesh_mask, 'depth': wfiles[0], 'data': maskfile}}
    variables = {'U': 'vozocrtx', 'V': 'vomecrty', 'W': 'vovecrtz', 'mask': 'fmask'}
    dimensions = {'U': {'lon': 'glamf', 'lat': 'gphif', 'depth': 'depthw', 'time': 'time_counter'},
                'V': {'lon': 'glamf', 'lat': 'gphif', 'depth': 'depthw', 'time': 'time_counter'},
                'W': {'lon': 'glamf', 'lat': 'gphif', 'depth': 'depthw', 'time': 'time_counter'},
                'mask': {'lon': 'glamf', 'lat': 'gphif', 'depth': 'depthw'}}

    fieldset = FieldSet.from_nemo(filenames, variables, dimensions, netcdf_decodewarning=False)
        
    for r in L_range:
        
        # for theta in theta_range:
        lonp.append(lon_0 + np.sin(theta)*r) 
        latp.append(lat_0 + np.cos(theta)*r)

    end_time = np.datetime64('2010-01-02') + np.timedelta64(1, 'h')
    times = [np.datetime64('2010-01-02')]*len(lonp)
    depp = np.zeros(len(lonp))
    pset = ParticleSet(fieldset, EnsembleParticle, lon=lonp, lat=latp, depth=depp, time=times)
        
    pset.execute([SampleField], 
                dt=delta(hours=1), endtime=end_time)

    R_values = compute_correlation_function(pset.u, pset.v, max_lag)    
    R_members[member-1, :] = np.array(R_values)

#%% Save the results
np.save('../data/spatial_correlations_Cape_Hatteras.npy', R_members)

#%% COmpute stats
L_r = np.zeros(50)
criterium = 0.1

for i in range(0,50):
    if np.all(R_members[i,:] >= criterium, axis=0):
        # skip this iteration
        L_r[i] = np.nan
    elif np.all(R_members[i,:] == 0, axis=0):
        L_r[i] = np.nan
    else:
        L = np.where(R_members[i, :] < criterium)[0][0]
        L_r[i] = L
    
L_r_std = np.nanstd(L_r)*0.01
L_r_mean = np.nanmean(L_r)*0.01

print(f"Mean L_r: {L_r_mean}, std L_r: {L_r_std}")

lag_array = np.arange(0, max_lag)*0.01

#%%

# Define the exponential decay function
def exp_decay(x, b):
    return np.exp(-b * x)

# and you make lag array 2 dimensional so it matches R_members
lag_array_2d = np.zeros_like(R_members)

for i in range(50):
    lag_array_2d[i, :] = lag_array
    
# Perform the curve fitting
params, covariance = curve_fit(exp_decay, lag_array_2d.ravel(), R_members.ravel())

# Extract the fitted parameters
b_fitted = params

# Calculate the decorrelation length as 1/b
decorrelation_length = 1 / b_fitted

print(f"Decorrelation length: {decorrelation_length}")

#%% plot space correlations

for member in range(1, 51):
    plt.plot(lag_array, R_members[member - 1, :], color='black', linewidth=0.5, alpha=0.5)


plt.plot(lag_array, exp_decay(lag_array, b_fitted), color='red', label=r'$e^{-x/L_L}$', lw=1)
# plt.grid()

plt.axhline(y=0, color='k', linestyle='--', label='Zero correlation')

plt.axvline(x=decorrelation_length, color='r', linestyle='--', lw=1, label=r'$L_L = 0.41^o$')

plt.xlabel(r'$L$ (degrees)')
plt.ylabel(r'$\rho(L)$')
plt.legend()
plt.ylim([-0.75, 1])
plt.xlim([0, 2])
plt.title('Cape Hatteras Spatial Correlation')
plt.savefig('../figs/cape_hatteras_spatial_correlation.png', dpi=300)






#%% 
R_members_temp = np.load('../data/temporal_correlations_Cape_Hatteras.npy')

lag_tiempo = np.arange(0, 60)

# and you make lag array 2 dimensional so it matches R_members
lag_tiempo_2d = np.zeros_like(R_members_temp)

for i in range(50):
    lag_tiempo_2d[i, :] = lag_tiempo
    
# Perform the curve fitting
params_temp, _ = curve_fit(exp_decay, lag_tiempo_2d.ravel(), R_members_temp.ravel())

# Extract the fitted parameters
b_temp = params_temp

# Calculate the decorrelation length as 1/b
decorrelation_time = 1 / b_temp

print(f"Decorrelation time: {decorrelation_time} days")

# %% land MAsk

mask_file = '/storage/shared/oceanparcels/input_data/NEMO_Ensemble/GRID/NATL025-CJMenobs01_byte_mask.nc'
mask = xr.open_dataset(mask_file, decode_times=False)


#%% Combined Plot

fig, axs = plt.subplots(1, 3, figsize=(12, 4), subplot_kw={'projection': cartopy.crs.PlateCarree()})

# Map plot
ax = axs[0]
# ax = plt.axes(projection=cartopy.crs.PlateCarree())
ax.set_extent([-77, -71, 33, 38], crs=cartopy.crs.PlateCarree())
# ax.add_feature(cartopy.feature.LAND, zorder=0, color='gray')
ax.pcolormesh(mask['nav_lon'], mask['nav_lat'], mask['tmask'][0,0,:,:], cmap='Greens_r', alpha=0.5)
gl = ax.gridlines(draw_labels=True, zorder=0, linestyle='--', linewidth=0.5)

gl.xlabels_top = False
gl.ylabels_right = False

hexfunc.plot_hexagons(ax, hexagons_grid, colors='r', draw_edges=True, fill_polygons=False)

ax.scatter(pset.lon, pset.lat, 
           transform=cartopy.crs.PlateCarree(), s=0.2, c='blue', label='Spatial Correlations')

ax.scatter(loc1_lon, loc1_lat, 
           transform=cartopy.crs.PlateCarree(), s=15, c='darkred', label='Temporal Correlations')

ax.legend(shadow=True, loc='upper left')
ax.set_title(r'$\mathbf{A}$  Cape Hatteras', loc='left')

# Spatial correlations plot
ax = fig.add_subplot(1, 3, 2)  # Add subplot without projection
for member in range(1, 51):
    ax.plot(lag_array, R_members[member - 1, :], color='black', linewidth=0.5, alpha=0.5)

ax.plot(lag_array, exp_decay(lag_array, b_fitted), color='blue', label=r'$e^{-x/L_L}$', lw=2)
ax.axhline(y=0, color='k', linestyle='--', label='Zero correlation')
ax.axvline(x=decorrelation_length, color='green', linestyle='--', lw=2, label=r'$L_L = 0.41^o$')
ax.set_xlabel(r'$L$ (degrees)')
ax.set_ylabel(r'$\rho(L)$')
ax.legend(shadow=True)
ax.set_ylim([-1, 1])
ax.set_xlim([0, 2])
ax.set_title(r'$\mathbf{B}$  Spatial Correlation', loc='left')

# Temporal correlations plot
ax = fig.add_subplot(1, 3, 3)  # Add subplot without projection
for member in range(50):
    ax.plot(lag_tiempo, R_members_temp[member, :], color='black', linewidth=0.5, alpha=0.5)

ax.plot(lag_tiempo, exp_decay(lag_tiempo, b_temp), color='darkred', label=r'$e^{-x/\tau_L}$', lw=2)
ax.axvline(x=decorrelation_time, color='green', linestyle='--', lw=2, label=r'$\tau_L = 41$ days')
ax.axhline(y=0, color='k', linestyle='--', label='Zero correlation')

ax.set_xlabel(r'$\tau$ (days)')
ax.set_ylabel(r'$\rho(\tau)$')
ax.legend(shadow=True)
ax.set_xlim([0, 59])
ax.set_ylim([-1, 1])
ax.set_title(r'$\mathbf{C}$  Temporal Correlation', loc='left')

plt.tight_layout()
plt.savefig('../figs/figS0_cape_hatteras_correlations.png', dpi=300)
# %%
