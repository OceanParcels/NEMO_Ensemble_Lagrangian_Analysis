from parcels import FieldSet, ParticleSet, AdvectionRK4_3D, AdvectionRK4, ParticleFile, JITParticle, ErrorCode, StateCode
import numpy as np
import xarray as xr
from glob import glob
from datetime import timedelta as delta
from argparse import ArgumentParser


def run_expt(member):
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

    lonp = -75
    latp = 30
    depthp = 1

#     times = np.arange(np.datetime64('2010-01-02'), np.datetime64('2010-01-31'), np.timedelta64(1, 'h'))
    times = np.arange(0, 86400*30, 3600)
    pset = ParticleSet(fieldset, JITParticle, lon=lonp, lat=latp, 
                       depth=1)
    len(pset)

    def SubmergeParticle(particle, fieldset, time):
        particle.depth = 0
        AdvectionRK4(particle, fieldset, time)  # perform a 2D advection because vertical flow will always push up in this case
        particle.time = time + particle.dt  # to not trigger kernels again, otherwise infinite loop
        particle.set_state(StateCode.Success)

    outfile = f"Pspacetime_1p_{member:03d}.zarr"
    pfile = ParticleFile(outfile, pset, outputdt=delta(days=1), chunks=(len(pset), 1))
    pset.execute(AdvectionRK4_3D, dt=delta(hours=1), output_file=pfile, 
                 recovery={ErrorCode.ErrorThroughSurface: SubmergeParticle})


if __name__ == "__main__":
    for member in range(1, 51):
        run_expt(member)
