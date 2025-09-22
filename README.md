# Subsurface\_MHW

This is the package for subsurface MHW tracking method, which has been writtern by matlab and python.

The current version of this package is organized as follows.



①Detect pointwise MHWs following Hobday et al., 2016 (main\_detect\_hobday.m and detect.m).



②Correlation coeffcients are calculated within fixed boxes for each point (cal\_correlation.m).



③The correlation-based KNN (K-nearest neighbor) cluster is carried out by KNN\_time\_para\_2020\_example.py and a MPI parrallel is adopted in this algorithm.



④Track the MHWs in the time domain (mhw\_overlap\_tracks\_4d\_125.m and mhw\_overlap\_tracks\_4d\_125.py).


Usage of mhw\_overlap\_tracks\_4d\_125.py:

python mhw\_overlap\_tracks\_4d\_125.py --start\_year 1993 --end\_year 2020 --data\_path 'your/mhw\_data/path' --mask 'map\_mask.mat' --alpha 0.5 --longitude\_length 360 --latitude\_length 120 --depth 20 --cut\_off 5 --judgment 'region' --output 'MHW\_tracks\_3d\_200m\_1x1\_60\_125\_coef\_0.6\_final'



Welcome to use and modify!



# Reference

Sun D., F. -R. Li, Z. Jing, S. -J. Hu, B. -H. Zhang, 2023: Frequent Marine Heatwaves Hidden below the Surface of the Global Ocean. Nature Geoscience, https://www.nature.com/articles/s41561-023-01325-w.



# Open Results

Due to the substantial computational demands of the third step (correlation-based KNN), and in response to reader requests, we are making publicly available the results derived from the 1993-2020 GLORYS dataset presented in our NG paper, which include four-dimensional marine heatwave (MHW) events in upper 200 m of the global ocean. We are also releasing results based on the 1982-2020 OISST dataset, which contain surface-only MHW events (three-dimensional MHW events) used in our PiO paper.

These datasets are accessible for download via the following Baidu Disk link:

🔗 Link: https://pan.baidu.com/s/1qKzg2XA4g7HwbYn4oFBGpQ?pwd=mhwt

🔑 Password: mhwt


Should you have any questions, please feel free to contact me at sundi@ouc.edu.cn.

