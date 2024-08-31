#%%
from parcels import FieldSet, ParticleSet, AdvectionRK4_3D, ParticleFile, JITParticle, StatusCode, Variable
import numpy as np
from glob import glob
from datetime import timedelta as delta
from argparse import ArgumentParser

# Parse command line arguments
p = ArgumentParser()
p.add_argument('-m', '--member', type=int, default=1, help='Member number')
p.add_argument('-d', '--depth', type=int, default=1, help='Initialization depth')
p.add_argument('-dr', '--delta_r', type=float, help='delta r')
p.add_argument('-w', '--weeks', type=int, help='Weeks span')
args = p.parse_args()

# Extract member and weeks from command line arguments
member = args.member
depth = args.depth
delta_r = args.delta_r
weeks = args.weeks
# member = 1
# depth = 1
# delta_r = 1.75
# weeks = 12

# Set simulation parameters
location = 'Cape_Hatteras'
start_time = np.datetime64('2010-01-02')
end_time = start_time + np.timedelta64(weeks, 'W')
release_time_range = np.arange(start_time, end_time, delta(days=1))

# If start_time = 2010-01-02, sim_end_time = 2015-12-31
sim_end_time = start_time + np.timedelta64(2189, 'D')

outfile = f"/storage/shared/oceanparcels/output_data/data_Claudio/NEMO_Ensemble/{location}/temporal_connectivity/dep_{depth:01d}/{location}_dep{depth:01d}_m{member:03d}.zarr"
print("Output file: ", outfile)

#%% Import Fieldset
data_path = '/storage/shared/oceanparcels/input_data/NEMO_Ensemble/'
ufiles = sorted(glob(f"{data_path}NATL025-CJMCYC3.{member:03d}-S/1d/*/NATL025*U.nc")) #Load all years
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

#%% Declare the ParticleSet
# N_particles = 10000
# subset = int(N_particles/len(release_release_time_range)) # Number of particles to release at dt = 1 day
span = delta_r
L_range = np.linspace(-span, span, 25) # This L_range and theta_range makes that alway there are 1001 particles
theta_range = np.arange(0, 2*np.pi, np.pi/20)


lon_0, lat_0 = (-73.61184289610455, 35.60913368957989)
lonp = [lon_0]
latp = [lat_0]


for r in L_range:
    for theta in theta_range:
        lonp.append(lon_0 + np.sin(theta)*r) 
        latp.append(lat_0 + np.cos(theta)*r)

subset = len(lonp)
print(f"Number of particles to release at each time: {subset}")

times = []
lon_release = []
lat_release = []
for t in release_time_range:
    times += [t]*subset
    lon_release += lonp
    lat_release += latp

times = np.array(times)
lon_release = np.array(lon_release)
lat_release = np.array(lat_release)

N_particles = len(times)
print(f"Number of particles: {N_particles}")

depp = np.ones(N_particles)*depth
hex_ids = [590726022320619519]*N_particles
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
    

pset = ParticleSet(fieldset, EnsembleParticle, lon=lon_release, lat=lat_release, 
                   depth=depp, time=times, hexbin_id=hex_ids)

# %%Declare Kernels
def KeepInOcean(particle, fieldset, time):
    if particle.state == StatusCode.ErrorThroughSurface:
        particle_ddepth = 1.0
        particle.state = StatusCode.Success


def CheckOutOfBounds(particle, fieldset, time):
    if particle.state == StatusCode.ErrorOutOfBounds:
        particle.delete()


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


pfile = ParticleFile(outfile, pset, outputdt=delta(days=1), chunks=(len(pset), 1))

# Error handling kernels go at the end of the kernel list
pset.execute([AdvectionRK4_3D, SampleField, KeepInOcean, CheckOutOfBounds], 
             dt=delta(hours=1), endtime=sim_end_time,
             output_file=pfile)


# %%
