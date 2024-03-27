#%%
from parcels import FieldSet, ParticleSet, AdvectionRK4_3D, AdvectionRK4, ParticleFile, JITParticle, StatusCode
import numpy as np
from glob import glob
from datetime import timedelta as delta
from argparse import ArgumentParser

p = ArgumentParser()
p.add_argument('-m', '--member', type=int, default=1, help='Member number')
p.add_argument('-y', '--year', type=int, default=2010, help='Year')

args = p.parse_args()

member = args.member
year = args.year

print(year, member)

data_path = '/storage/shared/oceanparcels/input_data/NEMO_Ensemble/'
ufiles = sorted(glob(f"{data_path}NATL025-CJMCYC3.{member:03d}-S/1d/{year}/NATL025*U.nc"))
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

step=1/32.
min_lon = -74
max_lon = -72
min_lat = 35
max_lat = 37
min_depth = 0
max_depth = 1600
z_step = 100
X, Y, Z = np.meshgrid(np.arange(min_lon, max_lon, step), 
                         np.arange(min_lat, max_lat, step),
                         np.arange(min_depth, max_depth, z_step))
n_particles = len(X.flatten())
start_times = [np.datetime64(f'{year}-01-02')]*n_particles

pset = ParticleSet(fieldset, JITParticle, lon=X, lat=Y, depth=Z, time=start_times)

def KeepInOcean(particle, fieldset, time):
    if particle.state == StatusCode.ErrorThroughSurface:
        particle_ddepth = 1.0
        particle.state = StatusCode.Success

outfile = f"/storage/shared/oceanparcels/output_data/data_Claudio/NEMO_Ensemble/{year}/PGS_{year}_{member:03d}.zarr"
pfile = ParticleFile(outfile, pset, outputdt=delta(days=1), chunks=(len(pset), 1))

pset.execute([AdvectionRK4_3D, KeepInOcean], 
             dt=delta(hours=1), 
             output_file=pfile)

