from parcels import FieldSet, ParticleSet, JITParticle, Variable, ParticleFile
import numpy as np
import xarray as xr
from glob import glob
from datetime import timedelta as delta
import h3
import pickle
import sys
sys.path.append('../functions')
import hexbin_functions as hexfunc


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

list_locations = {'Cape_Hatteras': (-74.0, 35.5),
                  'Canary_Current': (-15, 31.5)}


with open('../data/hexgrid_no_coast.pkl', 'rb') as f:
    hexagons_grid = pickle.load(f)
    
grid = hexfunc.hexGrid(hexagons_grid)

max_lag = 60
start_time = np.datetime64('2010-01-02')
end_time = np.datetime64('2010-04-02')
location = "Canary_Current"
N_particles = 1

# define time range with 1 day intervals
time_range = np.arange(start_time, end_time, delta(days=1))

loc1_lon = list_locations[location][0] # -74.0
loc1_lat = list_locations[location][1] # 35.5

# Find the hexagon containing the location
loc1_hex = h3.geo_to_h3(loc1_lat, loc1_lon, 3)
loc1_lat, loc1_lon = h3.h3_to_geo(loc1_hex)
lon_0 = loc1_lon
lat_0 = loc1_lat

lonp = [lon_0]*N_particles
latp = [lat_0]*N_particles

times = [start_time]*N_particles
depp = np.zeros(N_particles)

R_members = np.zeros((50, max_lag))

for member in range(1, 51):
    print(f"Member {member}")
    
    outfile = f"/storage/shared/oceanparcels/output_data/data_Claudio/NEMO_Ensemble/analysis/correlations/temporal/temp_correlation{location}_m{member:03d}.zarr"

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
        
    
    depp = np.zeros(len(lonp))
    pset = ParticleSet(fieldset, EnsembleParticle, lon=lonp, lat=latp, depth=depp, time=times)
    
    
    pfile = ParticleFile(outfile, pset, outputdt=delta(days=1), chunks=(len(pset), 1))
    
    pset.execute([SampleField], 
                dt=delta(days=1), endtime=end_time, output_file=pfile)


    pset1 = xr.open_zarr(outfile)
    

    R_values = compute_correlation_function(pset1.u[0,:], pset1.v[0,:], max_lag)    
    R_members[member-1, :] = np.array(R_values)
    
# Save the correlation values to a file
np.save(f"../data/temporal_correlations_{location}.npy", R_members)
