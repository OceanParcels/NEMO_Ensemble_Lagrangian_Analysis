"""test code to get geojson strings for sample of 100 minimum paths between two stations"""

from geojson import Feature, MultiLineString, dump, FeatureCollection
import h3
from sourcecode.core import adjacencygraph as ag
from sourcecode.core import connectivityhelper as ch

hex_res = 3
data_folder = '/Users/dmanral/Desktop/Analysis/TARA/Task7D/'
depth = 100
min_accept_temp = 11.85
max_accept_temp = 28.85

# get stations lat-lon
source_code = '66SUR'
destination_code = 'OA002'

master_grids_list = ch.get_all_grids_hex_ids(data_folder + 'MasterHexList_Res3.npy')

s_hex, d_hex = ch.get_station_hexes_from_code(data_folder + 'AllStations_Tara.xls', hex_res, source_code,
                                              destination_code)
try:
    s_index, d_index = master_grids_list.index(s_hex), master_grids_list.index(d_hex)
except KeyError:
    print("Source/destination sourcecode not present in the domain. Recheck values")
    raise

atlantic_graph = ag.create_temp_min_max_graph(data_folder + 't{0}m/Annual_Avg_DomainAdjacency_csr.npz'.format(depth),
                                              data_folder + 't{0}m/Annual_Avg_MinTemperature_csr.npz'.format(depth),
                                              data_folder + 't{0}m/Annual_Avg_MaxTemperature_csr.npz'.format(depth),
                                              min_accept_temp, max_accept_temp, None)


def get_path_locations(path):
    centers = [h3.h3_to_geo(master_grids_list[ind]) for ind in path]
    new_c = list(map(lambda sub: (sub[1], sub[0]), centers))  # [(t[1], t[0]) for t in centers]
    return new_c


forward_path = ag.get_shortest_path(atlantic_graph, s_index, d_index)
forward_paths = ag.get_shortest_paths_subset(atlantic_graph, s_index, d_index, len(forward_path))
forward_locations = [get_path_locations(p) for p in forward_paths]

backward_path = ag.get_shortest_path(atlantic_graph, d_index, s_index)
backward_paths = ag.get_shortest_paths_subset(atlantic_graph, d_index, s_index, len(backward_path))
backward_locations = [get_path_locations(p) for p in backward_paths]

features = []
features.append(
    Feature(geometry=MultiLineString(forward_locations),
            properties={"source": source_code,
                        "destination": destination_code,
                        "depth": depth,
                        "time": len(forward_path) - 1}))
features.append(
    Feature(geometry=MultiLineString(backward_locations),
            properties={"source": destination_code,
                        "destination": source_code,
                        "depth": depth,
                        "time": len(backward_path) - 1}))

feature_collection = FeatureCollection(features)

with open(data_folder + 'testfile.geojson', 'w') as f:
    dump(feature_collection, f)
