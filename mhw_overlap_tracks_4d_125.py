import numpy as np
from datetime import date
from scipy import io as sio
import h5py
from scipy import ndimage
import time
import argparse

# input parameters
parser = argparse.ArgumentParser()
parser.add_argument("--start_year", type=int, default=1993, help="the start year of the data")
parser.add_argument("--end_year", type=int, default=2020, help="the end year of the data")
parser.add_argument('--data_path', type=str, default='./glorys_mean/glorys_mean/', help="data storage address")
parser.add_argument("--mask", type=str, default='map_mask.mat', help='land mask file')
parser.add_argument("--alpha", type=float, default=0.5, help='threshold of the connected domain')
parser.add_argument("--longitude_length", type=int, default=360, help='the number of longitude grids')
parser.add_argument("--latitude_length", type=int, default=120, help='the number of latitude grids')
parser.add_argument("--depth", type=int, default=20, help='The number of depth grids')
parser.add_argument("--cut_off", type=int, default=5, help='mhw_tracks shorter than cut_off days will be removed')
parser.add_argument("--judgment", type=str, default='region', help='the criteria for determining a heat wave include "volume" and "region"')

opt = parser.parse_args()
START_YEAR = opt.start_year
END_YEAR = opt.end_year
DATA_PATH = opt.data_path
LAND_MASK = opt.mask
ALPHA = opt.alpha
LON_LEN = opt.longitude_length
LAT_LEN = opt.latitude_length
DEPTH = opt.depth
CUT_OFF = opt.cut_off
JUDGMENT = opt.judgment
earth_radius = 6378.137
pi180 = np.pi / 180


def data_load(year, data_path, land_mask):
    # load one year's data
    with h5py.File(data_path + str(year) + '/mask_now.mat') as F:
        mask = F['mask_now'][:]
    mask = np.transpose(mask, (3, 2, 1, 0))
    map_mask = sio.loadmat(land_mask)
    sst = map_mask['sst']
    for each_day in range(0, (date(year, 12, 31) - date(year, 1, 1)).days + 1):
        for each_depth in range(0, 20):
            depth_day_slice = mask[:, :, each_depth, each_day]
            depth_day_slice[(depth_day_slice == 1) & (sst == 1)] = 0
            mask[:, :, each_depth, each_day] = depth_day_slice
    return mask


def make_MHWs(mask, date_count_f, judgment='volume'):
    MHWs = []
    for day in range(0, mask.shape[-1]):
        MHWs_dict = {'day': day + 1 + date_count_f, 'zloc': [], 'yloc': [], 'xloc': []}
        cal_mask = mask[:, :, :, day]
        labels, num_mhws = ndimage.label(cal_mask, structure=ndimage.generate_binary_structure(3, 3))
        for mhw in range(0, num_mhws):
            mhw_ind = np.where(labels == (mhw + 1))

            # use volume measures whether it is a Marine heat wave
            if judgment == 'volume':
                volume = 0
                for mhw_cell in range(0, mhw_ind[0].shape[0]):
                    volume = volume + (earth_radius ** 2 * 1 * pi180 * (
                                np.sin(((mhw_ind[1][mhw_cell] + 1) + 1 / 2 + 30 - 90) * pi180) - np.sin(
                            ((mhw_ind[1][mhw_cell] + 1) - 1 / 2 + 30 - 90) * pi180))) * 10 * 1e-3
                if volume >= 1.0732e+04 * 10 * 1e-3 * 125:
                    MHWs_dict['zloc'].append([int(z) for z in mhw_ind[2]])
                    MHWs_dict['yloc'].append([int(y) for y in mhw_ind[1]])
                    MHWs_dict['xloc'].append([int(x) for x in mhw_ind[0]])

            # use occupied mesh ragion measures whether it is a Marine heat wave
            elif judgment == 'region':
                if len(mhw_ind[0]) >= 125:
                    MHWs_dict['zloc'].append([int(z) for z in mhw_ind[2]])
                    MHWs_dict['yloc'].append([int(y) for y in mhw_ind[1]])
                    MHWs_dict['xloc'].append([int(x) for x in mhw_ind[0]])
            else:
                raise ValueError("Invalid mhws judgment method value provided.")
        MHWs.append(MHWs_dict)
    return MHWs


def put_mhw_in_mesh(dims: tuple, x, y, z):
    mesh_mhw = np.zeros(dims)
    mesh_mhw[x, y, z] = 1
    return mesh_mhw


def mapping(t1_loc, t2_loc, lon_len, lat_len, depth, alpha):
    mapping_t1_t2 = np.zeros((len(t1_loc), len(t2_loc)))
    for each_t1mhw in range(len(t1_loc)):
        # each t2mhw intersects with the current t1mhw
        overlap = np.zeros(len(t2_loc))
        judge_t1 = put_mhw_in_mesh((lon_len, lat_len, depth), t1_loc[each_t1mhw]['xloc'][0],
                                   t1_loc[each_t1mhw]['yloc'][0], t1_loc[each_t1mhw]['zloc'][0])
        for each_t2mhw in range(len(t2_loc)):
            judge_t2 = put_mhw_in_mesh((lon_len, lat_len, depth), t2_loc[each_t2mhw]['xloc'][0],
                                       t2_loc[each_t2mhw]['yloc'][0], t2_loc[each_t2mhw]['zloc'][0])
            # ------------------- Important Part ! ----------------
            # compute overlap with a fixed t1mhw when loop for t2mhw
            overlap[each_t2mhw] = np.sum(np.logical_and(judge_t1, judge_t2)) / min((np.sum(judge_t2)),
                                                                                   (np.sum(judge_t1)))
            # -----------------------------------------------------
        idx = np.where(overlap >= alpha)[0]
        mapping_t1_t2[each_t1mhw, idx] = mapping_t1_t2[each_t1mhw, idx] + 1
    return mapping_t1_t2


def simple_and_spliting(mapping_t1_t2, search, x, y, z, day):
    # calculate the case where one t1_mhw corresponds to one or more t2_mhws
    t1_idx = np.where(np.sum(mapping_t1_t2, axis=1) > 0)[0]
    for idx_1 in t1_idx:
        t2_idx = np.where(mapping_t1_t2[idx_1] > 0)[0]
        temp_mhw_x = []
        temp_mhw_y = []
        temp_mhw_z = []
        for idx_2 in t2_idx:
            temp_mhw_x = temp_mhw_x + x[idx_2]
            temp_mhw_y = temp_mhw_y + y[idx_2]
            temp_mhw_z = temp_mhw_z + z[idx_2]
        search[idx_1]['day'].append(day)
        search[idx_1]['xloc'].append(temp_mhw_x)
        search[idx_1]['yloc'].append(temp_mhw_y)
        search[idx_1]['zloc'].append(temp_mhw_z)
    return search


def merging(mapping_t1_t2, search, loc_now, lon_len, lat_len, depth):
    # calculate the case where one t2_mhw corresponds to one or more t1_mhws
    t2_idx = np.where(np.sum(mapping_t1_t2, axis=0) > 1)[0]
    if len(t2_idx) > 0:
        for idx_2 in t2_idx:
            t1_idx = np.where(mapping_t1_t2[:, idx_2] > 0)[0]

            new_fusion = []
            new_fusion_dict = {'zloc': [], 'yloc': [], 'xloc': []}
            new_fusion_loc = []
            include_x = loc_now[idx_2]['xloc'][0]
            include_y = loc_now[idx_2]['yloc'][0]
            include_z = loc_now[idx_2]['zloc'][0]
            mean_x, mean_y, mean_z = [], [], []

            for idx_1, value_1 in np.ndenumerate(t1_idx):
                # calculate the fraction of t1_mhw in t2_mhw
                past_2d = put_mhw_in_mesh((lon_len, lat_len, depth),
                                          search[value_1]['xloc'][-2],
                                          search[value_1]['yloc'][-2],
                                          search[value_1]['zloc'][-2])

                now_2d = put_mhw_in_mesh((lon_len, lat_len, depth),
                                         search[value_1]['xloc'][-1],
                                         search[value_1]['yloc'][-1],
                                         search[value_1]['zloc'][-1])

                past_now_idx = np.where(np.logical_and(past_2d, now_2d))

                nd_x = np.array(include_x)
                nd_y = np.array(include_y)
                nd_z = np.array(include_z)
                # use different negative values to represent different t1_mhws in t2_mhw
                for i in range(past_now_idx[0].shape[0]):
                    indices = np.nonzero(np.all(np.array([nd_x, nd_y, nd_z]) == (
                        [past_now_idx[0][i]], [past_now_idx[1][i]], [past_now_idx[2][i]]), axis=0))
                    nd_x[indices] = nd_y[indices] = nd_z[indices] = -(idx_1[0] + 1)
                include_x = nd_x.tolist()
                include_y = nd_y.tolist()
                include_z = nd_z.tolist()
                # calculate the mean central point of each t1_mhw
                mean_x.append(np.mean(search[value_1]['xloc'][-2]))
                mean_y.append(np.mean(search[value_1]['yloc'][-2]))
                mean_z.append(np.mean(search[value_1]['zloc'][-2]))

            # calculate other part of t2_mhw
            include_x_0 = np.array(include_x)[np.where(np.array(include_x) > -1)]
            include_y_0 = np.array(include_y)[np.where(np.array(include_y) > -1)]
            include_z_0 = np.array(include_z)[np.where(np.array(include_z) > -1)]
            others = np.where(np.array(include_x) >= 0)[0]

            if include_x_0.shape[0] > 0:
                distance = []
                # calculate the distance from each other point to the center of each t1_mhw
                for i in range(len(mean_x)):
                    distance.append(np.sqrt((include_x_0 - mean_x[i]) ** 2 + (include_y_0 - mean_y[i]) ** 2 + (
                            include_z_0 - mean_z[i]) ** 2))
                min_distance = np.argmin(distance, axis=0)
                # add each other point to the nearest t1_mhw
                for each_point in range(len(min_distance)):
                    include_x[others[each_point]] = -(min_distance[each_point] + 1)
                    include_y[others[each_point]] = -(min_distance[each_point] + 1)
                    include_z[others[each_point]] = -(min_distance[each_point] + 1)

            for idx_1, value_1 in np.ndenumerate(t1_idx):
                # in the case of a merge t1_mhw does not split
                if len(search[value_1]['xloc'][-1]) == len(include_x):
                    index_include_x = np.where(np.array(include_x) == -(idx_1[0] + 1))[0]
                    search[value_1]['xloc'][-1] = (np.array(search[value_1]['xloc'][-1])[index_include_x]).tolist()
                    index_include_y = np.where(np.array(include_y) == -(idx_1[0] + 1))[0]
                    search[value_1]['yloc'][-1] = (np.array(search[value_1]['yloc'][-1])[index_include_y]).tolist()
                    index_include_z = np.where(np.array(include_z) == -(idx_1[0] + 1))[0]
                    search[value_1]['zloc'][-1] = (np.array(search[value_1]['zloc'][-1])[index_include_z]).tolist()
                # in the case of a merge t1_mhw  split
                else:
                    # extract1 is t1mhw after spliting, so extract1 is spuerset of extract2
                    extract1 = put_mhw_in_mesh((lon_len, lat_len, depth),
                                               search[value_1]['xloc'][-1],
                                               search[value_1]['yloc'][-1],
                                               search[value_1]['zloc'][-1])
                    # extract2 is current t2mhw in loop
                    extract2 = put_mhw_in_mesh((lon_len, lat_len, depth),
                                               loc_now[idx_2]['xloc'][0],
                                               loc_now[idx_2]['yloc'][0],
                                               loc_now[idx_2]['zloc'][0])

                    indice1 = np.where((extract1 == 1) & (extract2 == 0))
                    indice2 = np.where((extract1 == 1) & (extract2 == 1))

                    index_include_x = [index for index, value in enumerate(include_x) if value == -(idx_1[0] + 1)]
                    search[value_1]['xloc'][-1] = np.concatenate(
                        (indice1[0], indice2[0][index_include_x])).tolist()
                    index_include_y = [index for index, value in enumerate(include_y) if value == -(idx_1[0] + 1)]
                    search[value_1]['yloc'][-1] = np.concatenate(
                        (indice1[1], indice2[1][index_include_y])).tolist()
                    index_include_z = [index for index, value in enumerate(include_z) if value == -(idx_1[0] + 1)]
                    search[value_1]['zloc'][-1] = np.concatenate(
                        (indice1[2], indice2[2][index_include_z])).tolist()

                    new_fusion_loc.append(idx_1[0])
                    new_fusion_dict['xloc'].append(indice2[0][index_include_x].tolist())
                    new_fusion_dict['yloc'].append(indice2[1][index_include_y].tolist())
                    new_fusion_dict['zloc'].append(indice2[2][index_include_z].tolist())
                    new_fusion.append(new_fusion_dict)
                    new_fusion_dict = {'zloc': [], 'yloc': [], 'xloc': []}

            # The following is for those unconnected according to distance
            if len(new_fusion) == 0:
                # if there is no new split, loop for all xloc and yloc
                # to find connectivity
                for idx_1, value_1 in np.ndenumerate(t1_idx):
                    split_fusion = put_mhw_in_mesh((lon_len, lat_len, depth),
                                                   search[value_1]['xloc'][-1],
                                                   search[value_1]['yloc'][-1],
                                                   search[value_1]['zloc'][-1])
                    # calculate current mhw connectivity
                    D_labels, D_num_features = ndimage.label(split_fusion,
                                                             structure=ndimage.generate_binary_structure(3, 2))
                    D_pixel_idx_list = ndimage.find_objects(D_labels)

                    # if find unconnected contours
                    if D_num_features > 1:
                        connect_region = []
                        for r in range(D_num_features):
                            connect_region.append(len(np.where(split_fusion[D_pixel_idx_list[r]] > 0)[0]))
                        max_mhw = max(connect_region)
                        max_idx = connect_region.index(max_mhw)
                        for r in range(D_num_features):
                            if r == max_idx:
                                continue
                            # get x, y, z for the current connectivity area
                            conn_indice = np.where(D_labels == (r + 1))
                            # check whether the current connected area belongs to other t1_mhws
                            for idx_11, value_11 in np.ndenumerate(t1_idx):
                                if idx_11 == idx_1:
                                    continue
                                # concatenate the current connected area with other t1_mhws
                                append_x = search[value_11]['xloc'][-1] + conn_indice[0].tolist()
                                append_y = search[value_11]['yloc'][-1] + conn_indice[1].tolist()
                                append_z = search[value_11]['zloc'][-1] + conn_indice[2].tolist()

                                other_t1mhw = put_mhw_in_mesh((lon_len, lat_len, depth),
                                                              search[value_11]['xloc'][-1],
                                                              search[value_11]['yloc'][-1],
                                                              search[value_11]['zloc'][-1])
                                D_other_labels, D_other_num_features = ndimage.label(other_t1mhw, structure=ndimage.generate_binary_structure(3, 2))
                                append_mhw = put_mhw_in_mesh((lon_len, lat_len, depth), append_x, append_y, append_z)
                                D_append_labels, D_append_num_features = ndimage.label(append_mhw, structure=ndimage.generate_binary_structure(3, 2))

                                if D_append_num_features <= D_other_num_features:
                                    # put the append_mhw to the other mhw
                                    search[value_11]['xloc'][-1] = append_x
                                    search[value_11]['yloc'][-1] = append_y
                                    search[value_11]['zloc'][-1] = append_z
                                    # extract3 is current connectivity area
                                    extract3 = put_mhw_in_mesh((lon_len, lat_len, depth), conn_indice[0],
                                                               conn_indice[1], conn_indice[2])

                                    extract4 = put_mhw_in_mesh((lon_len, lat_len, depth),
                                                               search[value_1]['xloc'][-1],
                                                               search[value_1]['yloc'][-1],
                                                               search[value_1]['zloc'][-1])
                                    # cut off the currently connected area from the original t1mhw
                                    ex_indice = np.where(np.logical_and(extract3 == 0, extract4 == 1))
                                    search[value_1]['xloc'][-1] = ex_indice[0].tolist()
                                    search[value_1]['yloc'][-1] = ex_indice[1].tolist()
                                    search[value_1]['zloc'][-1] = ex_indice[2].tolist()
                                    break

            else:
                # if there are new splits, loop for those old xloc and yloc
                # to find connectivity
                # ------------------------------------------------------------------------
                for idx_1, value_1 in np.ndenumerate(t1_idx):
                    # the current t1mhw does not split and merge at the same time
                    # the processing is the same as the previous code
                    if len(np.where(np.array(new_fusion_loc) == idx_1[0])[0]) == 0:
                        split_fusion = put_mhw_in_mesh((lon_len, lat_len, depth),
                                                       search[value_1]['xloc'][-1],
                                                       search[value_1]['yloc'][-1],
                                                       search[value_1]['zloc'][-1])
                        # calculate current mhw connectivity
                        D_labels, D_num_features = ndimage.label(split_fusion,
                                                                 structure=ndimage.generate_binary_structure(3, 2))
                        D_pixel_idx_list = ndimage.find_objects(D_labels)
                    # the current t1mhw split and merge at the same time
                    else:
                        split_fusion = put_mhw_in_mesh((lon_len, lat_len, depth),
                                                       new_fusion[np.where(np.array(new_fusion_loc) == idx_1[0])[0][0]][
                                                           'xloc'][0],
                                                       new_fusion[np.where(np.array(new_fusion_loc) == idx_1[0])[0][0]][
                                                           'yloc'][0],
                                                       new_fusion[np.where(np.array(new_fusion_loc) == idx_1[0])[0][0]][
                                                           'zloc'][0])

                        D_labels, D_num_features = ndimage.label(split_fusion, structure=ndimage.generate_binary_structure(3, 2))
                        D_pixel_idx_list = ndimage.find_objects(D_labels)

                    if D_num_features > 1:
                        connect_region = []
                        for r in range(D_num_features):
                            connect_region.append(len(np.where(split_fusion[D_pixel_idx_list[r]] > 0)[0]))
                        max_mhw = max(connect_region)
                        max_idx = connect_region.index(max_mhw)
                        for r in range(D_num_features):
                            if r == max_idx:
                                continue
                            # get x, y, z for the current connectivity area
                            conn_indice = np.where(D_labels == (r + 1))
                            # check whether the current connected area belongs to other t1_mhws
                            for idx_11, value_11 in np.ndenumerate(t1_idx):
                                if idx_11 == idx_1:
                                    continue
                                # concatenate the current connected area with other t1_mhws
                                append_x = search[value_11]['xloc'][-1] + conn_indice[0].tolist()
                                append_y = search[value_11]['yloc'][-1] + conn_indice[1].tolist()
                                append_z = search[value_11]['zloc'][-1] + conn_indice[2].tolist()

                                other_t1mhw = put_mhw_in_mesh((lon_len, lat_len, depth),
                                                              search[value_11]['xloc'][-1],
                                                              search[value_11]['yloc'][-1],
                                                              search[value_11]['zloc'][-1])
                                D_other_labels, D_other_num_features = ndimage.label(other_t1mhw,
                                                                                     structure=ndimage.generate_binary_structure(3, 2))
                                append_mhw = put_mhw_in_mesh((lon_len, lat_len, depth), append_x, append_y, append_z)
                                D_append_labels, D_append_num_features = ndimage.label(append_mhw,
                                                                                       structure=ndimage.generate_binary_structure(3, 2))

                                if D_append_num_features <= D_other_num_features:
                                    # put the append_mhw to the other mhw
                                    search[value_11]['xloc'][-1] = append_x
                                    search[value_11]['yloc'][-1] = append_y
                                    search[value_11]['zloc'][-1] = append_z
                                    # extract3 is current connectivity area
                                    extract3 = put_mhw_in_mesh((lon_len, lat_len, depth), conn_indice[0],
                                                               conn_indice[1], conn_indice[2])

                                    extract4 = put_mhw_in_mesh((lon_len, lat_len, depth),
                                                               search[value_1]['xloc'][-1],
                                                               search[value_1]['yloc'][-1],
                                                               search[value_1]['zloc'][-1])
                                    # cut off the currently connected area from the original t1mhw
                                    ex_indice = np.where(np.logical_and(extract3 == 0, extract4 == 1))
                                    search[value_1]['xloc'][-1] = ex_indice[0].tolist()
                                    search[value_1]['yloc'][-1] = ex_indice[1].tolist()
                                    search[value_1]['zloc'][-1] = ex_indice[2].tolist()
                                    break
    return search


def tracking(current_MHWs, search, lon_len, lat_len, depth, alpha, tracks):
    for day_in_MHWs in range(0, len(current_MHWs)):
        current_day = current_MHWs[day_in_MHWs]['day']
        mhw_xloc = current_MHWs[day_in_MHWs]['xloc']
        mhw_yloc = current_MHWs[day_in_MHWs]['yloc']
        mhw_zloc = current_MHWs[day_in_MHWs]['zloc']
        # deal with the first day
        if current_day == 1:
            # put each mhw into search,
            for each_mhw in range(len(mhw_xloc)):
                search_dict = {'day': [], 'zloc': [], 'yloc': [], 'xloc': []}
                search_dict['day'].append(current_day)
                search_dict['xloc'].append(mhw_xloc[each_mhw])
                search_dict['yloc'].append(mhw_yloc[each_mhw])
                search_dict['zloc'].append(mhw_zloc[each_mhw])
                search.append(search_dict)

        # deal with other days
        else:
            # put each mhw from current day into loc_now
            loc_now = []
            for each_t2mhw in range(len(mhw_xloc)):
                loc_now_dict = {'zloc': [], 'yloc': [], 'xloc': []}
                loc_now_dict['xloc'].append(mhw_xloc[each_t2mhw])
                loc_now_dict['yloc'].append(mhw_yloc[each_t2mhw])
                loc_now_dict['zloc'].append(mhw_zloc[each_t2mhw])
                loc_now.append(loc_now_dict)
            # put each mhw from the day before the current day into loc_old
            loc_old = []
            for each_t1mhw in range(len(search)):
                loc_old_dict = {'day': [], 'zloc': [], 'yloc': [], 'xloc': []}
                loc_old_dict['day'].append(search[each_t1mhw]['day'][-1])
                loc_old_dict['xloc'].append(search[each_t1mhw]['xloc'][-1])
                loc_old_dict['yloc'].append(search[each_t1mhw]['yloc'][-1])
                loc_old_dict['zloc'].append(search[each_t1mhw]['zloc'][-1])
                loc_old.append(loc_old_dict)

            mapping_t1_t2 = mapping(loc_old, loc_now, lon_len, lat_len, depth, alpha)
            search = simple_and_spliting(mapping_t1_t2, search, mhw_xloc, mhw_yloc, mhw_zloc, current_day)
            search = merging(mapping_t1_t2, search, loc_now, lon_len, lat_len, depth)

            # for new tracks
            mhw_xloc = [value for index, value in enumerate(mhw_xloc) if np.sum(mapping_t1_t2, axis=0)[index] == 0]
            mhw_yloc = [value for index, value in enumerate(mhw_yloc) if np.sum(mapping_t1_t2, axis=0)[index] == 0]
            mhw_zloc = [value for index, value in enumerate(mhw_zloc) if np.sum(mapping_t1_t2, axis=0)[index] == 0]
            if len(mhw_xloc) > 0:
                for new_track in range(len(mhw_xloc)):
                    search_dict = {'day': [], 'zloc': [], 'yloc': [], 'xloc': []}
                    search_dict['day'].append(current_day)
                    search_dict['xloc'].append(mhw_xloc[new_track])
                    search_dict['yloc'].append(mhw_yloc[new_track])
                    search_dict['zloc'].append(mhw_zloc[new_track])
                    search.append(search_dict)

            # for moved tracks
            # when the mhw no longer lasts, move to tracks
            moved = []
            for mhw in range(len(search)):
                moved.append(search[mhw]['day'][-1] <= current_day - 1)
                if moved[mhw]:
                    search[mhw]['xloc'] = [list(map(lambda x: x + 1, row)) for row in search[mhw]['xloc']]
                    search[mhw]['yloc'] = [list(map(lambda x: x + 1, row)) for row in search[mhw]['yloc']]
                    search[mhw]['zloc'] = [list(map(lambda x: x + 1, row)) for row in search[mhw]['zloc']]
                    tracks.append(search[mhw])

            # remove tracks from open track array
            search = [value for value, flag in zip(search, moved) if not flag]
    return search, tracks


def main(start_year, end_year, data_path, land_mask, lon_len, lat_len, depth, alpha, cut_off, judgment):
    tracks = []
    search = []
    date_count = 0
    for current_year in range(start_year, end_year + 1):
        print('current_year:', current_year)
        load_start = time.time()
        current_year_mask = data_load(current_year, data_path, land_mask)
        current_MHWs = make_MHWs(current_year_mask, date_count, judgment)
        load_end = time.time()
        print('load_time:', load_end - load_start)
        search, tracks = tracking(current_MHWs, search, lon_len, lat_len, depth, alpha, tracks)
        date_count = date_count + (date(current_year, 12, 31) - date(current_year, 1, 1)).days + 1
    # Add tracks that are still in the search array at the end of the
    # time-series
    for mhw in range(len(search)):
        search[mhw]['xloc'] = [list(map(lambda x: x + 1, row)) for row in search[mhw]['xloc']]
        search[mhw]['yloc'] = [list(map(lambda x: x + 1, row)) for row in search[mhw]['yloc']]
        search[mhw]['zloc'] = [list(map(lambda x: x + 1, row)) for row in search[mhw]['zloc']]
        tracks.append(search[mhw])

    # %%%%%%%%%%%%%%%%%%%% remove tracks shorter than cut_off days %%%%%%%%%%%%%%
    short = []
    for mhw_track in range(len(tracks)):
        short.append(len(tracks[mhw_track]['day']) < cut_off)
    tracks = [value for value, flag in zip(tracks, short) if not flag]
    np.savez_compressed('./test/MHW_tracks_3d_200m_1x1_60_125_coef_0.6_final', tracks)


if __name__ == '__main__':

    main(START_YEAR, END_YEAR, DATA_PATH, LAND_MASK, LON_LEN, LAT_LEN, DEPTH, ALPHA, CUT_OFF, JUDGMENT)

