# Subsurface_MHW
This is the package for subsurface MHW tracking method, which is now being improved and updating by open software Python. Users can look forward to subsequent versions.


The current version of this package is organized as follows.


①Detect pointwise MHWs following Hobday et al., 2016 (main_detect_hobday.m and detect.m).


②Correlation coeffcients are calculated within fixed boxes for each point (cal_correlation.m).


③The KNN (K-nearest neighbor) cluster is carried out by KNN_time_para_2020_example.py and a MPI parrallel is adopted in this algorithm.


④Track the MHWs in the time domain (mhw_overlap_tracks_4d_125.m).


Welcome to use and modify!


# Reference
SUN D., F. -R. Li*, Z. Jing, S. -J. Hu, B. -H. Zhang, 2023: Frequent Marine Heatwaves Hidden below the Surface of the Global Ocean. Nature Geoscience, https://www.nature.com/articles/s41561-023-01325-w.