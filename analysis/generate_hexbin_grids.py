#%%
import numpy as np
import matplotlib.pyplot as plt
import xarray as xr

import h3
from scipy.interpolate import griddata
import pickle
from argparse import ArgumentParser

import sys
sys.path.append('../functions')
import hexbin_functions as hexfunc
#%%
# p = ArgumentParser()
# p.add_argument('-hr', '--hex_res', type=int, default=3, help='hexagon resolution')

# args = p.parse_args()

hex_res = 3 # args.hex_res


def get_coastal_nodes(landmask):
    """Function that detects the coastal nodes, i.e. the ocean nodes directly
    next to land. Computes the Laplacian of landmask.

    - landmask: the land mask built using `make_landmask`, where land cell = 1
                and ocean cell = 0.

    Output: 2D array array containing the coastal nodes, the coastal nodes are
            equal to one, and the rest is zero.
    """
    mask_lap = np.roll(landmask, -1, axis=0) + np.roll(landmask, 1, axis=0)
    mask_lap += np.roll(landmask, -1, axis=1) + np.roll(landmask, 1, axis=1)
    mask_lap -= 4*landmask
    coastal = np.ma.masked_array(landmask, mask_lap < 0)
    coastal = coastal.mask.astype('int')

    return coastal

analysis_domain = {"type": "Polygon",
                   "coordinates": [
                     [
  [73.43488716788852, -78.2549013096737],
  [54.57488797743113, -61.675459259113424],
  [52.453121315460066, -57.07005868951306],
  [50.95491032470406, -69.04410017047302],
  [30.75815419710834, -82.39976182231416],
  [30.162694581994558, -96.33109854535418],
  [27.8469723526252, -98.86406885863467],
  [22.536311461282153, -99.78514897255444],
  [18.001588273049236, -96.33109854535418],
  [13.347226126214665, -85.73867723527385],
  [8.601446991378992, -82.6300318507937],
  [8.829058069873284, -63.05707942999305],
  [-0.5768525742674768, -50.967902934792534],
  [-6.778286399290792, -37.151701225992326],
  [-13.573480027977041, -41.06629171015243],
  [-21.041840919672637, -42.447911881032894],
  [-21.792126579578806, 15.925540338649768],
  [-17.78471594728832, 13.277435011129484],
  [-12.226765214776748, 15.580135295929693],
  [8.829058069873284, 9.24770951272879],
  [8.259774873887878, -10.440377922311626],
  [13.906693603943978, -14.470103420711524],
  [22.6426133979418, -15.045778491911989],
  [31.251576239453954, -5.604707324230958],
  [39.40966407441758, -5.1441672672710865],
  [53.76615469626921, 9.477979541209123],
  [61.054602757311045, 8.672034441529178],
  [73.2035529895725, 27.09363671992969],
  [73.43488716788852, -78.2549013096737]
]]
}

grid_raw = h3.polyfill(analysis_domain, hex_res)

mask_file = '/storage/shared/oceanparcels/input_data/NEMO_Ensemble/GRID/NATL025-CJMenobs01_byte_mask.nc'
mask = xr.open_dataset(mask_file, decode_times=False)

# Extract the coordinates of the centers of the hexagons
hexgrid = list(grid_raw)
hexgrid_centres = np.array([h3.h3_to_geo(hex_) for hex_ in hexgrid])
hexgrid_centers_coordinates = (hexgrid_centres[:,1], hexgrid_centres[:,0])

# Declare arrays for lon and lat and mask values of model cells
depth = 0
mask_lons = mask['nav_lon'].values
mask_lats = mask['nav_lat'].values
mask_tland = mask['tmask'][0,depth,:,:].values
flipped_mask = (~mask_tland.astype(bool)).astype(int)

# Interpolate center of hexgrid to mask values

land_val_hexgrid = griddata((mask_lons.ravel(),mask_lats.ravel()), mask_tland.ravel(), hexgrid_centers_coordinates, method='nearest')

# discard hexagons where land_val_hexgrid is 0 (land)
grid_no_land = [hex_ for hex_, land in zip(hexgrid, land_val_hexgrid) if land == 1]

grid_no_land = set(grid_no_land) # making it a set to plot with hexfunc.plot_hexagons

#update the arrays with the centers of the grid witth no landcells
hexgrid_centres = np.array([h3.h3_to_geo(hex_) for hex_ in grid_no_land])
hexgrid_centers_coordinates = (hexgrid_centres[:,1], hexgrid_centres[:,0])

coastal_cells = get_coastal_nodes(mask['tmask'][0,0,:,:].values)
coastal_val_hexgrid = griddata((mask_lons.ravel(),mask_lats.ravel()), coastal_cells.ravel(), hexgrid_centers_coordinates, method='nearest')

grid_no_coast = [hex_ for hex_, coast in zip(grid_no_land, coastal_val_hexgrid) if coast == 0]

hexgrid_no_coast_centres = np.array([h3.h3_to_geo(hex_) for hex_ in grid_no_coast])

with open(f'../data/hexgrid_no_coast_h{hex_res}.pkl', 'wb') as f:
  pickle.dump(set(grid_no_coast), f)

with open(f'../data/hexgrid_land_h{hex_res}.pkl', 'wb') as f:
  pickle.dump(grid_land, f)
#%%
