import numpy as np
from datetime import date
from scipy import io as sio
import h5py
from scipy import ndimage
import os

depth = 20
# for k in range(1992, 2021):
#     mask_now = np.zeros((360, 120, depth, (date(k, 12, 31) - date(k, 1, 1)).days + 1))
#     for i in range(1, (date(k, 12, 31) - date(k, 1, 1)).days + 2):
#         load_filename = './' + str(k) + '/test', + "{:03d}".format(i) + '.mat'
#         mask = sio.loadmat(load_filename)
#         mask_now[:, :, :, i] = mask['mask']
#     save_filename = './' + str(k) + '/mask_now.mat'
#     file = h5py.File(save_filename, "w")
#     file.create_dataset('mask_now', data=mask_now)
#     file.close()

with h5py.File('./1994/mask_now.mat') as F:
    mask = F['mask_now'][:]

# for i in range(1993, 2021):
#     mask_file_name = './' + str(i) + '/mask_now.mat'
#     with h5py.File(mask_file_name) as F:
#         mask_now = F['mask_now'][:]
#         mask = np.concatenate((mask, mask_now), axis=3)
mask = np.transpose(mask,(3,2,1,0))
map_mask = sio.loadmat('map_mask.mat')
sst = map_mask['sst']
for k in range(0, (date(1993, 12, 31) - date(1993, 1, 1)).days + 1):
    for m in range(0, 20):
        extract = mask[:, :, m, k]
        extract[(extract == 1) & (sst == 1)] = 0
        mask[:, :, m, k] = extract


r = 6378.137
MHWs = []
mask_daily = mask.copy()
for i in range(0, mask_daily.shape[-1]):
    MHWs_dict = {'day': i + 1, 'xloc': [], 'yloc': [], 'zloc': []}
    cal_mask = mask_daily[:, :, :, i]
    labels, num_features = ndimage.label(cal_mask, structure=np.ones((3, 3, 3)))
    pixel_idx_list = ndimage.find_objects(labels)
    for j in range(0, num_features):
        region_slice = pixel_idx_list[j]
        mhw_ind = np.where(labels == (j+1))
        #region = cal_mask[region_slice[0], region_slice[1], region_slice[2]]
        # mhw_ind is a tuple
        #mhw_ind = np.nonzero(region)
        volume = 0
        for k in range(0, mhw_ind[0].shape[0]):
            volume = volume + (r ** 2 * 1 * np.pi / 180 * (np.sin(((mhw_ind[1][k] + 1) + 1 / 2 + 30 - 90) * np.pi / 180) - np.sin(((mhw_ind[1][k] + 1) - 1 / 2 + 30 - 90) * np.pi / 180))) * 10 *1e-3

        if volume >= 1.0732e+04*10*1e-3*125:

            MHWs_dict['zloc'].append([np.float32(z) for z in mhw_ind[2]])
            MHWs_dict['yloc'].append([np.float32(y) for y in mhw_ind[1]])
            MHWs_dict['xloc'].append([np.float32(x) for x in mhw_ind[0]])


    MHWs.append(MHWs_dict)

# ---------------------------------------------------------------------
# --------------------------- MHW tracking ----------------------------

# Input the parameter
alpha = 0.5

lon_len = 360
lat_len = 120
depth = 20

cut_off = 5
search = []

tracks = []

pi180 = np.pi/180
earth_radius = 6378.137
s_tb = np.ones((3, 3))
s_tb[[0, 0, -1, -1], [0, -1, 0, -1]] = 0
s_m = np.ones((3, 3))
structure = np.concatenate((np.expand_dims(s_tb, axis=0), np.expand_dims(s_m, axis=0), np.expand_dims(s_tb, axis=0)))

# ------------------------------ Beginning --------------------------------
for i in range(len(MHWs)):
    day = MHWs[i]['day']
    print('day', day)
    # 取出某一天的所有热浪
    mhw_xloc = MHWs[i]['xloc']
    mhw_yloc = MHWs[i]['yloc']
    mhw_zloc = MHWs[i]['zloc']
    # if the first day, open tracks
    if i == 0:
        for i2 in range(len(mhw_xloc)):
            search_dict = {'day': [], 'xloc': [], 'yloc': [], 'zloc': []}
            search_dict['day'].append(day)
            search_dict['xloc'].append(mhw_xloc[i2])
            search_dict['yloc'].append(mhw_yloc[i2])
            search_dict['zloc'].append(mhw_zloc[i2])
            search.append(search_dict)
    else:
        # if not the first day, put the xloc and yloc into loc_now
        loc_now = []
        for i2 in range(len(mhw_xloc)):
            loc_now_dict = {'xloc': [], 'yloc': [], 'zloc': []}
            loc_now_dict['xloc'].append(mhw_xloc[i2])
            loc_now_dict['yloc'].append(mhw_yloc[i2])
            loc_now_dict['zloc'].append(mhw_zloc[i2])
            loc_now.append(loc_now_dict)

        count = np.zeros(len(mhw_xloc))  # for counting times that idx is used

        for i2 in range(len(search)):
            # find xloc and yloc in open tracks from the previous day
            if (day - search[i2]['day'][-1]) == 1:
                loc_old = []
                loc_old.append(search[i2]['xloc'][-1])
                loc_old.append(search[i2]['yloc'][-1])
                loc_old.append(search[i2]['zloc'][-1])

                overlap = np.zeros(len(loc_now))
                judge1 = np.zeros((lon_len, lat_len, depth))

                for k in range(len(loc_old[0])):
                    judge1[int(loc_old[0][k]), int(loc_old[1][k]), int(loc_old[2][k])] = 1

                # loop for all mhws at t2
                for i3 in range(len(loc_now)):
                    judge2 = np.zeros((lon_len, lat_len, depth))
                    loc_now_x = loc_now[i3]['xloc'][0]
                    loc_now_y = loc_now[i3]['yloc'][0]
                    loc_now_z = loc_now[i3]['zloc'][0]
                    for k2 in range(len(loc_now_x)):
                        judge2[int(loc_now_x[k2]), int(loc_now_y[k2]), int(loc_now_z[k2])] = 1

                    # ------------------- Important Part ! ----------------
                    # compute overlap with a fixed t1 when loop for t2
                    overlap[i3] = np.sum(np.logical_and(judge1, judge2)) / min((np.sum(judge2)), (np.sum(judge1)))
                    # -----------------------------------------------------------------------------------------------

                idx = np.where(overlap >= alpha)[0]
                # =========================================================
                # maybe find more than one mhws at t2 overlap>0.5
                # for splitting
                if len(idx) > 0:
                    if len(idx) > 1:
                        print('---------------------------------------------------------')
                        print('get in spliting')
                        search[i2]['day'].append(day)
                        search[i2]['xloc'].append(mhw_xloc[idx[0]])
                        search[i2]['yloc'].append(mhw_yloc[idx[0]])
                        search[i2]['zloc'].append(mhw_zloc[idx[0]])

                        for i4 in range(1, idx.shape[0]):
                            search[i2]['xloc'][-1] = search[i2]['xloc'][-1] + mhw_xloc[idx[i4]]
                            search[i2]['yloc'][-1] = search[i2]['yloc'][-1] + mhw_yloc[idx[i4]]
                            search[i2]['zloc'][-1] = search[i2]['zloc'][-1] + mhw_zloc[idx[i4]]

                    else:
                        search[i2]['day'].append(day)
                        search[i2]['xloc'].append(mhw_xloc[idx[0]])
                        search[i2]['yloc'].append(mhw_yloc[idx[0]])
                        search[i2]['zloc'].append(mhw_zloc[idx[0]])

                #  ===================== End splitting ==================

                # =========================================================
                # for merging
                count[idx] = count[idx] + 1

        #  ======================== Begin merging ===========================
        # find previous tracks
        count_det_mhw = np.where(count > 1)[0]
        if len(count_det_mhw) > 0:
            print('---------------------------------------------------------')
            print('get in merging')
            idx_now = count_det_mhw
            old = {}
            for m in range(idx_now.shape[0]):
                # find t2 MHW in which t1 tracks before
                c = 0
                loc_now_xloc = loc_now[idx_now[m]]['xloc'][0]
                loc_now_yloc = loc_now[idx_now[m]]['yloc'][0]
                loc_now_zloc = loc_now[idx_now[m]]['zloc'][0]
                for n in range(len(search)):
                    # 取出已经关联上的t2时的热浪
                    search_xloc = search[n]['xloc'][-1]
                    search_yloc = search[n]['yloc'][-1]
                    search_zloc = search[n]['zloc'][-1]

                    l = np.zeros((lon_len, lat_len, depth))
                    for j5 in range(len(loc_now_xloc)):
                        l[int(loc_now_xloc[j5]), int(loc_now_yloc[j5]), int(loc_now_zloc[j5])] = 1

                    ll = np.zeros((lon_len, lat_len, depth))
                    for j5 in range(len(search_xloc)):
                        ll[int(search_xloc[j5]), int(search_yloc[j5]), int(search_zloc[j5])] = 1

                    if day == search[n]['day'][-1] and ((len(loc_now_xloc) == len(search_xloc) and loc_now_xloc == search_xloc) or
                                                        sum(np.where(np.logical_and(l==1, ll==1))[0]) == sum(np.where(l==1)[0])):
                        old[(str(m), str(c))] = n
                        c = c + 1

            if m > 1:
                print('Please take care:', i)

            # for number conservation
            for i5 in range(idx_now.shape[0]):
                new_fusion = []
                new_fusion_dict = {'day': None, 'xloc': [], 'yloc': [], 'zloc': []}
                new_fusion_loc = []
                include_x = loc_now[idx_now[i5]]['xloc'][0]
                include_y = loc_now[idx_now[i5]]['yloc'][0]
                include_z = loc_now[idx_now[i5]]['zloc'][0]

                mean_x, mean_y, mean_z = [], [], []
                for j in range(len([key for key in old.keys() if key[0] == str(i5)])):
                    # part I for those overlapped area
                    past_2d = np.zeros((lon_len,lat_len,depth))
                    for j3 in range(len(search[old[(str(i5), str(j))]]['xloc'][-2])):
                        past_2d[int(search[old[(str(i5), str(j))]]['xloc'][-2][j3]),
                                int(search[old[(str(i5), str(j))]]['yloc'][-2][j3]),
                                int(search[old[(str(i5), str(j))]]['zloc'][-2][j3])] = 1

                    now_2d = np.zeros((lon_len,lat_len,depth))
                    for j3 in range(len(search[old[(str(i5), str(j))]]['xloc'][-1])):
                        now_2d[int(search[old[(str(i5), str(j))]]['xloc'][-1][j3]),
                               int(search[old[(str(i5), str(j))]]['yloc'][-1][j3]),
                               int(search[old[(str(i5), str(j))]]['zloc'][-1][j3])] = 1

                    past_now_idx = np.where(np.logical_and(past_2d, now_2d))
                    for j3 in range(past_now_idx[0].shape[0]):
                        test = np.logical_and(np.logical_and(include_x == past_now_idx[0][j3], include_y == past_now_idx[1][j3]), include_z == past_now_idx[2][j3])
                        include_x = [-(j + 1) if flag else value for value, flag in zip(include_x, test)]
                        include_y = [-(j + 1) if flag else value for value, flag in zip(include_y, test)]
                        include_z = [-(j + 1) if flag else value for value, flag in zip(include_z, test)]

                    mean_x.append(np.mean(search[old[(str(i5), str(j))]]['xloc'][-2]))
                    mean_y.append(np.mean(search[old[(str(i5), str(j))]]['yloc'][-2]))
                    mean_z.append(np.mean(search[old[(str(i5), str(j))]]['zloc'][-2]))

                #  part II for others
                include_x_0 = [x for x in include_x if x > -1]
                include_y_0 = [y for y in include_y if y > -1]
                include_z_0 = [z for z in include_z if z > -1]
                others = np.where(np.array(include_x) >= 0)[0]
                if len(include_x_0) > 0:
                    a = []
                    for i6 in range(len(mean_x)):
                        a.append(np.sqrt((include_x_0 - mean_x[i6]) ** 2 + (include_y_0 - mean_y[i6]) ** 2 + (include_z_0 - mean_z[i6]) ** 2))

                    t = np.argmin(a, axis=0)
                    for w in range(len(t)):
                        include_x[others[w]] = -(t[w] + 1)
                        include_y[others[w]] = -(t[w] + 1)
                        include_z[others[w]] = -(t[w] + 1)



                for j in range(len([key for key in old.keys() if key[0] == str(i5)])):
                    if len(search[old[(str(i5), str(j))]]['xloc'][-1]) == len(include_x):
                        index_include_x = [index for index, value in enumerate(include_x) if value == -(j + 1)]
                        search[old[(str(i5), str(j))]]['xloc'][-1] = [search[old[(str(i5), str(j))]]['xloc'][-1][index]
                                                                      for index in index_include_x]
                        index_include_y = [index for index, value in enumerate(include_y) if value == -(j + 1)]
                        search[old[(str(i5), str(j))]]['yloc'][-1] = [search[old[(str(i5), str(j))]]['yloc'][-1][index]
                                                                      for index in index_include_y]
                        index_include_z = [index for index, value in enumerate(include_z) if value == -(j + 1)]
                        search[old[(str(i5), str(j))]]['zloc'][-1] = [search[old[(str(i5), str(j))]]['zloc'][-1][index]
                                                                      for index in index_include_z]

                    else:
                        #  The next is for both split and merge e.g. day 285
                        extract1 = np.zeros((lon_len, lat_len, depth))
                        for j1 in range(len(search[old[(str(i5), str(j))]]['xloc'][-1])):
                            extract1[int(search[old[(str(i5), str(j))]]['xloc'][-1][j1]),
                                     int(search[old[(str(i5), str(j))]]['yloc'][-1][j1]),
                                     int(search[old[(str(i5), str(j))]]['zloc'][-1][j1])] = 1

                        extract2 = np.zeros((lon_len, lat_len, depth))
                        for j1 in range(len(loc_now[idx_now[i5]]['xloc'][0])):
                            extract2[int(loc_now[idx_now[i5]]['xloc'][0][j1]),
                                     int(loc_now[idx_now[i5]]['yloc'][0][j1]),
                                     int(loc_now[idx_now[i5]]['zloc'][0][j1])] = 1

                        indice1 = np.where((extract1 == 1) & (extract2 == 0))
                        indice2 = np.where((extract1 == 1) & (extract2 == 1))

                        index_include_x = [index for index, value in enumerate(include_x) if value == -(j + 1)]
                        search[old[(str(i5), str(j))]]['xloc'][-1] = np.concatenate((indice1[0], indice2[0][index_include_x])).tolist()
                        index_include_y = [index for index, value in enumerate(include_y) if value == -(j + 1)]
                        search[old[(str(i5), str(j))]]['yloc'][-1] = np.concatenate((indice1[1], indice2[1][index_include_y])).tolist()
                        index_include_z = [index for index, value in enumerate(include_z) if value == -(j + 1)]
                        search[old[(str(i5), str(j))]]['zloc'][-1] = np.concatenate((indice1[2], indice2[2][index_include_z])).tolist()


                        new_fusion_loc.append(j)
                        new_fusion_dict['xloc'].append(indice2[0][index_include_x].tolist())
                        new_fusion_dict['yloc'].append(indice2[1][index_include_y].tolist())
                        new_fusion_dict['zloc'].append(indice2[2][index_include_z].tolist())
                        new_fusion.append(new_fusion_dict)
                        new_fusion_dict = {'day': None, 'xloc': [], 'yloc': [], 'zloc': []}

                # ------------------------------------------------------------------------------------------------------------
                # The following is for those unconnected according to distance
                if len(new_fusion) == 0:
                    # if there is no new split, loop for all xloc and yloc
                    # to find connectivity
                    # ------------------------------------------------------------------------
                    for j in range(len([key for key in old.keys() if key[0] == str(i5)])):
                        split_fusion = np.zeros((lon_len, lat_len, depth))
                        for k3 in range(len(search[old[(str(i5), str(j))]]['xloc'][-1])):
                            split_fusion[int(search[old[(str(i5), str(j))]]['xloc'][-1][k3]),
                                         int(search[old[(str(i5), str(j))]]['yloc'][-1][k3]),
                                         int(search[old[(str(i5), str(j))]]['zloc'][-1][k3])] = 1

                        D_labels, D_num_features = ndimage.label(split_fusion, structure=structure)
                        D_pixel_idx_list = ndimage.find_objects(D_labels)

                        # if find unconnected contours
                        if D_num_features > 1:
                            b = []
                            for k7 in range(D_num_features):
                                b.append(len(np.where(split_fusion[D_pixel_idx_list[k7]] > 0)[0]))
                            p = max(b)
                            q = b.index(p)
                            for k4 in range(D_num_features):
                                if k4 != q:
                                    conn_indice = np.where(D_labels == (k4 + 1))
                                    for k5 in range(len([key for key in old.keys() if key[0] == str(i5)])):
                                        if k5 != j:
                                            append_x = search[old[(str(i5), str(k5))]]['xloc'][-1] + conn_indice[0].tolist()
                                            append_y = search[old[(str(i5), str(k5))]]['yloc'][-1] + conn_indice[1].tolist()
                                            append_z = search[old[(str(i5), str(k5))]]['zloc'][-1] + conn_indice[2].tolist()

                                            raw_search = np.zeros((lon_len, lat_len, depth))
                                            for k6 in range(len(search[old[(str(i5), str(k5))]]['xloc'][-1])):
                                                raw_search[int(search[old[(str(i5), str(k5))]]['xloc'][-1][k6]),
                                                           int(search[old[(str(i5), str(k5))]]['yloc'][-1][k6]),
                                                           int(search[old[(str(i5), str(k5))]]['zloc'][-1][k6])] = 1

                                            D_raw_labels, D_raw_num_features = ndimage.label(raw_search, structure=structure)

                                            test_fusion = np.zeros((lon_len, lat_len, depth))
                                            for k6 in range(len(append_x)):
                                                test_fusion[int(append_x[k6]), int(append_y[k6]), int(append_z[k6])] = 1
                                            D_1_labels, D_1_num_features = ndimage.label(test_fusion, structure=structure)


                                            if D_1_num_features <= D_raw_num_features:
                                                # put the merged to the new one
                                                search[old[(str(i5), str(k5))]]['xloc'][-1] = append_x
                                                search[old[(str(i5), str(k5))]]['yloc'][-1] = append_y
                                                search[old[(str(i5), str(k5))]]['zloc'][-1] = append_z

                                                # extract the x,y from the old
                                                extract3 = np.zeros((lon_len, lat_len, depth))
                                                for j2 in range(conn_indice[0].shape[0]):
                                                    extract3[conn_indice[0][j2], conn_indice[1][j2], conn_indice[2][j2]] = 1
                                                extract4 = np.zeros((lon_len, lat_len, depth))
                                                for j2 in range(len(search[old[(str(i5), str(j))]]['xloc'][-1])):
                                                    extract4[int(search[old[(str(i5), str(j))]]['xloc'][-1][j2]),
                                                             int(search[old[(str(i5), str(j))]]['yloc'][-1][j2]),
                                                             int(search[old[(str(i5), str(j))]]['zloc'][-1][j2])] = 1
                                                ex_indice = np.where(np.logical_and(extract3 == 0, extract4 == 1))
                                                search[old[(str(i5), str(j))]]['xloc'][-1] = ex_indice[0].tolist()
                                                search[old[(str(i5), str(j))]]['yloc'][-1] = ex_indice[1].tolist()
                                                search[old[(str(i5), str(j))]]['zloc'][-1] = ex_indice[2].tolist()
                                                break # if find a track that the rest is connected, then append into it and stop loop
                                                # i.e. append to the first one

                    b = []
                else:
                    # if there are new splits, loop for those old xloc and yloc
                    # to find connectivity
                    # ------------------------------------------------------------------------
                    for j in range(len([key for key in old.keys() if key[0] == str(i5)])):
                        if len(np.where(np.array(new_fusion_loc) == j)[0]) == 0:
                            split_fusion = np.zeros((lon_len, lat_len, depth))
                            for k3 in range(len(search[old[(str(i5), str(j))]]['xloc'][-1])):
                                split_fusion[int(search[old[(str(i5), str(j))]]['xloc'][-1][k3]),
                                             int(search[old[(str(i5), str(j))]]['yloc'][-1][k3]),
                                             int(search[old[(str(i5), str(j))]]['zloc'][-1][k3])] = 1
                            D_labels, D_num_features = ndimage.label(split_fusion, structure=structure)
                            D_pixel_idx_list = ndimage.find_objects(D_labels)

                        else:
                            split_fusion = np.zeros((lon_len, lat_len, depth))
                            for k3 in range(len(new_fusion[np.where(np.array(new_fusion_loc) == j)[0][0]]['xloc'][0])):
                                split_fusion[new_fusion[np.where(np.array(new_fusion_loc) == j)[0][0]]['xloc'][0][k3],
                                             new_fusion[np.where(np.array(new_fusion_loc) == j)[0][0]]['yloc'][0][k3],
                                             new_fusion[np.where(np.array(new_fusion_loc) == j)[0][0]]['zloc'][0][k3]] = 1

                            D_labels, D_num_features = ndimage.label(split_fusion, structure=structure)
                            D_pixel_idx_list = ndimage.find_objects(D_labels)

                        if len(D_pixel_idx_list) > 1:
                            b = []
                            for k7 in range(len(D_pixel_idx_list)):
                                b.append(len(np.where(split_fusion[D_pixel_idx_list[k7]] > 0)[0]))
                            p = max(b)
                            q = b.index(p)
                            for k4 in range(len(D_pixel_idx_list)):
                                if k4 != q:
                                    conn_indice =  np.where(D_labels == (k4 + 1))
                                    for k5 in range(len([key for key in old.keys() if key[0] == str(i5)])):
                                        if k5 != j:
                                            append_x = search[old[(str(i5), str(k5))]]['xloc'][-1] + conn_indice[0].tolist()
                                            append_y = search[old[(str(i5), str(k5))]]['yloc'][-1] + conn_indice[1].tolist()
                                            append_z = search[old[(str(i5), str(k5))]]['zloc'][-1] + conn_indice[2].tolist()

                                            raw_search = np.zeros((lon_len, lat_len, depth))
                                            for k6 in range(len(search[old[(str(i5), str(k5))]]['xloc'][-1])):
                                                raw_search[int(search[old[(str(i5), str(k5))]]['xloc'][-1][k6]),
                                                           int(search[old[(str(i5), str(k5))]]['yloc'][-1][k6]),
                                                           int(search[old[(str(i5), str(k5))]]['zloc'][-1][k6])] = 1

                                            D_raw_labels, D_raw_num_features = ndimage.label(raw_search, structure=structure)

                                            test_fusion = np.zeros((lon_len, lat_len, depth))
                                            for k6 in range(len(append_x)):
                                                test_fusion[int(append_x[k6]), int(append_y[k6]), int(append_z[k6])] = 1
                                            D_1_labels, D_1_num_features = ndimage.label(test_fusion, structure=structure)


                                            if D_1_num_features <= D_raw_num_features:
                                                # put the merged to the new one
                                                search[old[(str(i5), str(k5))]]['xloc'][-1] = append_x
                                                search[old[(str(i5), str(k5))]]['yloc'][-1] = append_y
                                                search[old[(str(i5), str(k5))]]['zloc'][-1] = append_z

                                                # extract the x,y from the old
                                                extract3 = np.zeros((lon_len, lat_len, depth))
                                                for j2 in range(conn_indice[0].shape[0]):
                                                    extract3[conn_indice[0][j2], conn_indice[1][j2], conn_indice[2][j2]] = 1
                                                extract4 = np.zeros((lon_len, lat_len, depth))
                                                for j2 in range(len(search[old[(str(i5), str(j))]]['xloc'][-1])):
                                                    extract4[int(search[old[(str(i5), str(j))]]['xloc'][-1][j2]),
                                                             int(search[old[(str(i5), str(j))]]['yloc'][-1][j2]),
                                                             int(search[old[(str(i5), str(j))]]['zloc'][-1][j2])] = 1
                                                ex_indice = np.where(np.logical_and(extract3 == 0, extract4 == 1))
                                                search[old[(str(i5), str(j))]]['xloc'][-1] = ex_indice[0].tolist()
                                                search[old[(str(i5), str(j))]]['yloc'][-1] = ex_indice[1].tolist()
                                                search[old[(str(i5), str(j))]]['zloc'][-1] = ex_indice[2].tolist()
                                                break # if find a track that the rest is connected, then append into it and stop loop
                                                # i.e. append to the first one

                    b = []
        # ================================ End Merging ==============================

        # =================================================================
        # for new tracks
        mhw_xloc = [value for index, value in enumerate(mhw_xloc) if count[index] == 0]
        mhw_yloc = [value for index, value in enumerate(mhw_yloc) if count[index] == 0]
        mhw_zloc = [value for index, value in enumerate(mhw_zloc) if count[index] == 0]
        if len(mhw_xloc) > 0:
            for i7 in range(len(mhw_xloc)):
                search_dict = {'day': [], 'xloc': [], 'yloc': [], 'zloc': []}
                search_dict['day'].append(day)
                search_dict['xloc'].append(mhw_xloc[i7])
                search_dict['yloc'].append(mhw_yloc[i7])
                search_dict['zloc'].append(mhw_zloc[i7])
                search.append(search_dict)

        # == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == =

        # == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == == =
        # for moved tracks
        moved = []
        for i2 in range(len(search)):
            moved.append(search[i2]['day'][-1] <= day - 1)
            if moved[i2]:
                search[i2]['xloc'] = [list(map(lambda x: x + 1, row)) for row in search[i2]['xloc']]
                search[i2]['yloc'] = [list(map(lambda x: x + 1, row)) for row in search[i2]['yloc']]
                search[i2]['zloc'] = [list(map(lambda x: x + 1, row)) for row in search[i2]['zloc']]
                tracks.append(search[i2])

        # remove tracks from open track array
        search = [value for value, flag in zip(search, moved) if not flag]

        # =================================================================

# Add tracks that are still in the search array at the end of the
# time-series
for i in range(len(search)):
    search[i]['xloc'] = [list(map(lambda x: x + 1, row)) for row in search[i]['xloc']]
    search[i]['yloc'] = [list(map(lambda x: x + 1, row)) for row in search[i]['yloc']]
    search[i]['zloc'] = [list(map(lambda x: x + 1, row)) for row in search[i]['zloc']]
    tracks.append(search[i])

# %%%%%%%%%%%%%%%%%%%% remove tracks shorter than cut_off days %%%%%%%%%%%%%%
short = []
for i in range(len(tracks)):
    short.append(len(tracks[i]['day']) < cut_off)
tracks = [value for value, flag in zip(tracks, short) if not flag]

np.save('./mhw_tracks/MHW_tracks_3d_200m_1x1_60_125_coef_0.6.npy', tracks)
