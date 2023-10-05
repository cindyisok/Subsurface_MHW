# Subsurface_MHW
This is the package for subsurface MHW tracking method.

We first detect pointwise MHWs following Hobday et al., 2016 using main_detect_hobday.m and detect.m.
Then, correlation coeffcients are calculated within fixed boxes for each point by cal_correlation.m.
The KNN (K-nearest neighbor) cluster is carried out by KNN_time_para_2020_example.py and a MPI parrallel is adopted in this algorithm.
Last, we track the MHWs in the time domain by mhw_overlap_tracks_4d_125.m.

Welcome to use and modify!
