import pandas as pd
import numpy as np
from sourcecode.core import adjacencygraph as ag
from sourcecode.core import connectivityhelper as ch
import h3
import time
from sklearn.preprocessing import normalize

home_folder = '/Users/dmanral/Desktop/Analysis/TARA/Task7D/'
# data_folder = home_folder + '2011_Lombard_Species/'

hex_res = 3
depth = 0


def fullconnectivity(species, mintemp, maxtemp, final_stations_master_indices, master_hex_ids):
    # to add constraints to the connectivity
    domain_adjacency_file = 't{0}m/Annual_Avg_DomainAdjacency_csr.npz'.format(depth)

    temp_constraint_range = np.NaN
    min_accept_temp = mintemp
    max_accept_temp = maxtemp

    # create graph
    if ~np.isnan(temp_constraint_range):
        print('Temp range: ', temp_constraint_range)
        atlantic_graph = ag.create_temp_range_graph(home_folder + domain_adjacency_file,
                                                    home_folder + 't{0}m/Annual_Avg_MinTemperature_csr.npz'.format(
                                                        depth),
                                                    home_folder + 't{0}m/Annual_Avg_MaxTemperature_csr.npz'.format(
                                                        depth),
                                                    temp_constraint_range,
                                                    None)
    elif ~np.isnan(max_accept_temp) and ~np.isnan(min_accept_temp):
        print('Min/Max Temp: ', min_accept_temp, max_accept_temp)
        atlantic_graph = ag.create_temp_min_max_graph(home_folder + domain_adjacency_file,
                                                      home_folder + 't{0}m/Annual_Avg_MinTemperature_csr.npz'.format(
                                                          depth),
                                                      home_folder + 't{0}m/Annual_Avg_MaxTemperature_csr.npz'.format(
                                                          depth),
                                                      min_accept_temp, max_accept_temp,
                                                      None)
    else:
        atlantic_graph = ag.create_simple_graph(home_folder + domain_adjacency_file,
                                                None)

    # get min_T paths for all pairs- forward and backward- 2 d matrix.
    North_con_density = np.zeros(len(master_hex_ids))
    N_enable = np.zeros(len(master_hex_ids), dtype=bool)
    South_con_density = np.zeros(len(master_hex_ids))
    S_enable = np.zeros(len(master_hex_ids), dtype=bool)

    station_count = len(final_stations_master_indices)
    st_time = time.time()
    for i in range(station_count):
        s_idx = final_stations_master_indices[i]
        for j in range(i, station_count):
            d_idx = final_stations_master_indices[j]
            if s_idx != d_idx:
                f_path = ag.get_shortest_path(atlantic_graph, s_idx, d_idx)
                if f_path:
                    # for ind in f_path:
                    #     North_con_density[ind] += 1
                    #     N_enable[ind] = True
                    forward_paths = ag.get_shortest_paths_subset(atlantic_graph, s_idx, d_idx)
                    for path in forward_paths:
                        for ind in path:
                            North_con_density[ind] += 1
                            N_enable[ind] = True

                b_path = ag.get_shortest_path(atlantic_graph, d_idx, s_idx)
                if b_path:
                    #     for ind in b_path:
                    #         South_con_density[ind] += 1
                    #         S_enable[ind] = True
                    backward_paths = ag.get_shortest_paths_subset(atlantic_graph, d_idx, s_idx)
                    for path in backward_paths:
                        for ind in path:
                            South_con_density[ind] += 1
                            S_enable[ind] = True
            # else:
            #     if ag.check_if_edge_exists(atlantic_graph, s_idx, d_idx):
            #         if i != j: #stations with same Hex ID
            #             min_T_matrix[j][i] = 0
            # min_T_matrix[i][j] = 0
    print("Full time:", time.time() - st_time)
    hex_centers = [h3.h3_to_geo(str(h)) for h in master_hex_ids]
    lats = np.array([x1[0] for x1 in hex_centers])
    lons = np.array([x1[1] for x1 in hex_centers])

    df = pd.DataFrame(data={'S_N_gridId': np.array(master_hex_ids)[N_enable],
                            'S_N_latitudes': lats[N_enable],
                            'S_N_longitudes': lons[N_enable],
                            'S_N_Density': North_con_density[N_enable]},
                      columns=['S_N_gridId', 'S_N_latitudes', 'S_N_longitudes', 'S_N_Density'])
    df.to_csv(home_folder + 'Min100_S_N_PathsDensity_d{0}.csv'.format(depth))

    df2 = pd.DataFrame(data={'N_S_gridId': np.array(master_hex_ids)[S_enable],
                             'N_S_latitudes': lats[S_enable],
                             'N_S_longitudes': lons[S_enable],
                             'N_S_Density': South_con_density[S_enable]},
                       columns=['N_S_gridId', 'N_S_latitudes', 'N_S_longitudes', 'N_S_Density'])
    df2.to_csv(home_folder + 'Min100_N_S_PathsDensity_d{0}.csv'.format(depth))


def main1():
    # species_info = pd.read_csv(data_folder + '2011_lombard_forams.csv',
    #                            delimiter=';|,', header=0)
    # Get all sorted stations- station codes between -100 and 20 Longitude
    stations = pd.read_csv(home_folder + 'Tara_Stations_hexId_sorted.csv', header=0, index_col=0)
    stations_code, stations_hex = stations.index.values, stations['H3Id_res3'].values

    master_hex_ids = ch.get_all_grids_hex_ids(home_folder + 'MasterHexList_Res3.npy')
    # map station to master hex (some stations lie in the same hex- same connectivity)
    # domain_mask = np.in1d(stations_hex, master_hex_ids)
    # final_stations_code = stations_code[domain_mask]
    # final_stations_hex = stations_hex[domain_mask]

    final_stations_master_indices = np.array([master_hex_ids.index(h) for h in stations_hex])

    # [fullconnectivity(entry['Species'], entry['MinTemp'], entry['MaxTemp'], len(final_stations_hex),
    #                   final_stations_master_indices, final_stations_code) for index, entry in species_info.iterrows()]
    fullconnectivity(None, np.nan, np.nan, final_stations_master_indices, master_hex_ids)


def main():
    file = home_folder + 'Min100_N_S_PathsDensity_d{0}.csv'.format(depth)
    data = pd.read_csv(file)
    norm_density = data['N_S_Density'] / np.sum(data['N_S_Density'])
    print(np.sum(norm_density))
    data['N_S_norm_dens'] = norm_density.apply(lambda x: '%.5f' % x)
    data.to_csv(file)

    file = home_folder + 'Min100_S_N_PathsDensity_d{0}.csv'.format(depth)
    data = pd.read_csv(file)
    norm_density = data['S_N_Density'] / np.sum(data['S_N_Density'])
    print(np.sum(norm_density))
    data['S_N_norm_dens'] = norm_density.apply(lambda x: '%.5f' % x)
    data.to_csv(file)


if __name__ == '__main__':
    main()
