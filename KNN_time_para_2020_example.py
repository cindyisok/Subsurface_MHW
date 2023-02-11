# -*- coding: utf-8 -*-
"""
Created on Fri Nov 25 19:16:44 2022

@author: Cindy
"""

import numpy as np
import scipy.io as scio
import cc3d
import h5py
from mpi4py import MPI
from sys import stdout

#
# load data
#
data = scio.loadmat('map_ecco2.mat')
sst = data['sst']

len_yr = 366
depth = 20
beta = 0.6


# ###############################################################################
#
# main function of maximum correlation cluster
#
def maximum_correlation_cluster(rank,sst, mask):
    
    i1 = 6
    i2 = 360-5+1
    j1 = 6
    j2 = 120-5+1
    mask_3d = np.zeros([360, 120, depth])
    
    for i in range(i1, i2):
       # print('############## successful! ####################')
        print('my rank is' + str(rank))
        print(i)
        stdout.flush()
        
        for j in range(j1, j2):
            
            if sst[i-1, j-1] == 0:
                
                datapath='/home/jingzhao/sundi/data/r_correlation_ecco2_real_trend/'
                data = scio.loadmat(datapath+'r_'+str(i)+'_'+str(j)+'.mat')
                
                r = data['r']
                r[np.isnan(r)] = 0
            
                for n in range(depth):
                    
                    t = np.where(r[61 + n*121 - 1, :] > beta)
                    
                    if len(t[0]) >= 1:
                        nt = len(t[0])
                        judge_knn = np.zeros((nt, 3)).astype('int')
    
                        [judge_knn[:, 0], judge_knn[:, 1], judge_knn[:, 2]] = \
                                         np.unravel_index(t[0], (11, 11, depth), order='F')
    
                        mask_now = np.zeros((11, 11, depth))
    
                        for h in range(np.shape(judge_knn)[0]):
                            mask_now[judge_knn[h, 0], judge_knn[h, 1], judge_knn[h, 2]] = 1
                        D1, N = cc3d.connected_components(mask_now, connectivity=26, return_N=True)
    
                        for h1 in range(1, N+1):
                            if D1[5, 5, n] == h1:
                                t_now = np.where(D1 == h1)
                                break
                            
                        judge_knn1 = np.zeros(
                            (np.shape(t_now)[1], 4)).astype('int')
                        for m in range(0, np.shape(t_now)[1]):
                            judge_knn1[m, 0] = t_now[0][m]
                            judge_knn1[m, 1] = t_now[1][m]
                            judge_knn1[m, 2] = t_now[2][m]
                            
                        judge_knn1[:, 0] = j - 5 + judge_knn1[:, 0] - 1
                        judge_knn1[:, 1] = i - 5 + judge_knn1[:, 1] - 1
    
                        for m in range(0, np.shape(t_now)[1]):
                            judge_knn1[m, 3] = mask[judge_knn1[m, 1], judge_knn1[m, 0], 
                                                        judge_knn1[m, 2]]
                        
                        if np.sum(judge_knn1[:, 3]) >= np.shape(t_now)[1]/2:
                            mask_3d[i-1, j-1, n] = 1
                            
                        
                    else:
                        mask_3d[i-1, j-1, n] = mask[i-1, j-1, n]
            
            
        #
        # the upper and lower boundary
        #
        for j in range(1, j1):
            
            if sst[i-1, j-1] == 0:
                
                datapath='/home/jingzhao/sundi/data/r_correlation_ecco2_real_trend/'
                data = scio.loadmat(datapath+'r_'+str(i)+'_'+str(j)+'.mat')
                
                r = data['r']
                r[np.isnan(r)] = 0
                
                for n in range(depth):
                    
                    t = np.where(r[j + 11*5 + n*121 - 1,:] > beta);
                    
                    if len(t[0]) >= 1:
                        nt = len(t[0])
                        judge_knn = np.zeros((nt, 3)).astype('int')
    
                        [judge_knn[:, 0], judge_knn[:, 1], judge_knn[:, 2]] = \
                                         np.unravel_index(t[0], (11, 11, depth), order='F')
    
                        mask_now = np.zeros((11, 11, depth))
    
                        for h in range(np.shape(judge_knn)[0]):
                            mask_now[judge_knn[h, 0], judge_knn[h, 1], judge_knn[h, 2]] = 1
                        D1, N = cc3d.connected_components(mask_now, connectivity=26, return_N=True)
    
                        for h1 in range(1, N+1):
                            if D1[j-1, 5, n] == h1:
                                t_now = np.where(D1 == h1)
                                break
                            
                        judge_knn1 = np.zeros(
                            (np.shape(t_now)[1], 4)).astype('int')
                        for m in range(0, np.shape(t_now)[1]):
                            judge_knn1[m, 0] = t_now[0][m]
                            judge_knn1[m, 1] = t_now[1][m]
                            judge_knn1[m, 2] = t_now[2][m]
    
                        judge_knn1[:, 1] = i - 5 + judge_knn1[:, 1] - 1
    
                        for m in range(0, np.shape(t_now)[1]):
                            judge_knn1[m, 3] = mask[judge_knn1[m, 1], judge_knn1[m, 0], 
                                                        judge_knn1[m, 2]]
                        
                        if np.sum(judge_knn1[:, 3]) >= np.shape(t_now)[1]/2:
                            mask_3d[i-1, j-1, n] = 1
    
                    else:
                        mask_3d[i-1, j-1, n] = mask[i-1, j-1, n]
                        
                        

        for j in range(j2, 120+1):
            
            if sst[i-1, j-1] == 0:
                
                datapath='/home/jingzhao/sundi/data/r_correlation_ecco2_real_trend/'
                data = scio.loadmat(datapath+'r_'+str(i)+'_'+str(j)+'.mat')
                
                r = data['r']
                r[np.isnan(r)] = 0
    
                for n in range(depth):
                    
                    t = np.where(r[j-109 + 11*5 + n*121 - 1,:] > beta);
                    
                    if len(t[0]) >= 1:
                        nt = len(t[0])
                        judge_knn = np.zeros((nt, 3)).astype('int')
    
                        [judge_knn[:, 0], judge_knn[:, 1], judge_knn[:, 2]] = \
                                         np.unravel_index(t[0], (11, 11, depth), order='F')
    
                        mask_now = np.zeros((11, 11, depth))
    
                        for h in range(np.shape(judge_knn)[0]):
                            mask_now[judge_knn[h, 0], judge_knn[h, 1], judge_knn[h, 2]] = 1
                        D1, N = cc3d.connected_components(mask_now, connectivity=26, return_N=True)
    
                        for h1 in range(1, N+1):
                            if D1[j-j2+6, 5, n] == h1:
                                t_now = np.where(D1 == h1)
                                break
                            
                        judge_knn1 = np.zeros(
                            (np.shape(t_now)[1], 4)).astype('int')
                        for m in range(0, np.shape(t_now)[1]):
                            judge_knn1[m, 0] = t_now[0][m]
                            judge_knn1[m, 1] = t_now[1][m]
                            judge_knn1[m, 2] = t_now[2][m]
    
                        judge_knn1[:, 0] = 120 - 10 + judge_knn1[:, 0] - 1
                        judge_knn1[:, 1] = i - 5 + judge_knn1[:, 1] - 1
    
                        for m in range(0, np.shape(t_now)[1]):
                            judge_knn1[m, 3] = mask[judge_knn1[m, 1], judge_knn1[m, 0], 
                                                        judge_knn1[m, 2]]
                        
                        if np.sum(judge_knn1[:, 3]) >= np.shape(t_now)[1]/2:
                            mask_3d[i-1, j-1, n] = 1
    
                    else:
                        mask_3d[i-1, j-1, n] = mask[i-1, j-1, n]
                    
                                       

    for i in range(1, i1):
        
        for j in range(j1, j2):           
            
            if sst[i-1, j-1] == 0:
                
                datapath='/home/jingzhao/sundi/data/r_correlation_ecco2_real_trend/'
                data = scio.loadmat(datapath+'r_'+str(i)+'_'+str(j)+'.mat')
                
                r = data['r']
                r[np.isnan(r)] = 0
    
                for n in range(depth):
                    
                    t = np.where(r[6 + 11*(i-1) + n*121 - 1, :] > beta)
                    
                    if len(t[0]) >= 1:
                        nt = len(t[0])
                        judge_knn = np.zeros((nt, 3)).astype('int')
    
                        [judge_knn[:, 0], judge_knn[:, 1], judge_knn[:, 2]] = \
                                         np.unravel_index(t[0], (11, 11, depth), order='F')
    
                        mask_now = np.zeros((11, 11, depth))
    
                        for h in range(np.shape(judge_knn)[0]):
                            mask_now[judge_knn[h, 0], judge_knn[h, 1], judge_knn[h, 2]] = 1
                        D1, N = cc3d.connected_components(mask_now, connectivity=26, return_N=True)
    
                        for h1 in range(1, N+1):
                            if D1[5, i-1, n] == h1:
                                t_now = np.where(D1 == h1)
                                break
                            
                        judge_knn1 = np.zeros((np.shape(t_now)[1], 4)).astype('int')

                        for m in range(0, np.shape(t_now)[1]):
                            judge_knn1[m, 0] = t_now[0][m]
                            judge_knn1[m, 1] = t_now[1][m]
                            judge_knn1[m, 2] = t_now[2][m]
    
                        judge_knn1[:, 0] = j - 5 + judge_knn1[:, 0] - 1
    
                        for m in range(0, np.shape(t_now)[1]):
                            judge_knn1[m, 3] = mask[judge_knn1[m, 1], judge_knn1[m, 0], 
                                                        judge_knn1[m, 2]]
                        
                        if np.sum(judge_knn1[:, 3]) >= np.shape(t_now)[1]/2:
                            mask_3d[i-1, j-1, n] = 1
    
                    else:
                        mask_3d[i-1, j-1, n] = mask[i-1, j-1, n]
                        
                        

        for j in range(1, j1):
            
            if sst[i-1, j-1] == 0:
                
                datapath='/home/jingzhao/sundi/data/r_correlation_ecco2_real_trend/'
                data = scio.loadmat(datapath+'r_'+str(i)+'_'+str(j)+'.mat')
                
                r = data['r']
                r[np.isnan(r)] = 0
    
                for n in range(depth):
                    
                    t = np.where(r[j + 11*(i-1) + n*121 - 1, :] > beta)
                    
                    if len(t[0]) >= 1:
                        nt = len(t[0])
                        judge_knn = np.zeros((nt, 3)).astype('int')
    
                        [judge_knn[:, 0], judge_knn[:, 1], judge_knn[:, 2]] = \
                                         np.unravel_index(t[0], (11, 11, depth), order='F')
    
                        mask_now = np.zeros((11, 11, depth))
    
                        for h in range(np.shape(judge_knn)[0]):
                            mask_now[judge_knn[h, 0], judge_knn[h, 1], judge_knn[h, 2]] = 1
                        D1, N = cc3d.connected_components(mask_now, connectivity=26, return_N=True)
    
                        for h1 in range(1, N+1):
                            if D1[j-1, i-1, n] == h1:
                                t_now = np.where(D1 == h1)
                                break
                            
                        judge_knn1 = np.zeros(
                            (np.shape(t_now)[1], 4)).astype('int')
                        for m in range(0, np.shape(t_now)[1]):
                            judge_knn1[m, 0] = t_now[0][m]
                            judge_knn1[m, 1] = t_now[1][m]
                            judge_knn1[m, 2] = t_now[2][m]
    
    
                        for m in range(0, np.shape(t_now)[1]):
                            judge_knn1[m, 3] = mask[judge_knn1[m, 1], judge_knn1[m, 0], 
                                                        judge_knn1[m, 2]]
                        
                        if np.sum(judge_knn1[:, 3]) >= np.shape(t_now)[1]/2:
                            mask_3d[i-1, j-1, n] = 1
    
                    else:
                        mask_3d[i-1, j-1, n] = mask[i-1, j-1, n]
                        
        
        
        for j in range(j2, 120+1):
            
            if sst[i-1, j-1] == 0:
                
                datapath='/home/jingzhao/sundi/data/r_correlation_ecco2_real_trend/'
                data = scio.loadmat(datapath+'r_'+str(i)+'_'+str(j)+'.mat')
                
                r = data['r']
                r[np.isnan(r)] = 0
    
                for n in range(depth):
                    
                    t = np.where(r[j-109 + 11*(i-1) + n*121 - 1, :] > beta)
                    
                    if len(t[0]) >= 1:
                        nt = len(t[0])
                        judge_knn = np.zeros((nt, 3)).astype('int')
    
                        [judge_knn[:, 0], judge_knn[:, 1], judge_knn[:, 2]] = \
                                         np.unravel_index(t[0], (11, 11, depth), order='F')
    
                        mask_now = np.zeros((11, 11, depth))
    
                        for h in range(np.shape(judge_knn)[0]):
                            mask_now[judge_knn[h, 0], judge_knn[h, 1], judge_knn[h, 2]] = 1
                        D1, N = cc3d.connected_components(mask_now, connectivity=26, return_N=True)
    
                        for h1 in range(1, N+1):
                            if D1[j-j2+6, i-1, n] == h1:
                                t_now = np.where(D1 == h1)
                                break
                            
                        judge_knn1 = np.zeros((np.shape(t_now)[1], 4)).astype('int')
                        for m in range(0, np.shape(t_now)[1]):
                            judge_knn1[m, 0] = t_now[0][m]
                            judge_knn1[m, 1] = t_now[1][m]
                            judge_knn1[m, 2] = t_now[2][m]
    
                        judge_knn1[:, 0] = 120 - 10 + judge_knn1[:, 0] - 1
    
                        for m in range(0, np.shape(t_now)[1]):
                            judge_knn1[m, 3] = mask[judge_knn1[m, 1], judge_knn1[m, 0], 
                                                        judge_knn1[m, 2]]
                        
                        if np.sum(judge_knn1[:, 3]) >= np.shape(t_now)[1]/2:
                            mask_3d[i-1, j-1, n] = 1
    
                    else:
                        mask_3d[i-1, j-1, n] = mask[i-1, j-1, n]



    #
    # the right boundary
    #
    for i in range(i2, 360+1):
        for j in range(j1, j2):
            
            if sst[i-1, j-1] == 0:
                
                datapath='/home/jingzhao/sundi/data/r_correlation_ecco2_real_trend/'
                data = scio.loadmat(datapath+'r_'+str(i)+'_'+str(j)+'.mat')
                
                r = data['r']
                r[np.isnan(r)] = 0
    
                for n in range(depth):
                    
                    t = np.where(r[6 + 11*(i-350) + n*121 - 1, :] > beta)
                    
                    if len(t[0]) >= 1:
                        nt = len(t[0])
                        judge_knn = np.zeros((nt, 3)).astype('int')
    
                        [judge_knn[:, 0], judge_knn[:, 1], judge_knn[:, 2]] = \
                                         np.unravel_index(t[0], (11, 11, depth), order='F')
    
                        mask_now = np.zeros((11, 11, depth))
    
                        for h in range(np.shape(judge_knn)[0]):
                            mask_now[judge_knn[h, 0], judge_knn[h, 1], judge_knn[h, 2]] = 1
                        D1, N = cc3d.connected_components(mask_now, connectivity=26, return_N=True)
    
                        for h1 in range(1, N+1):
                            if D1[5, i-355+5, n] == h1:
                                t_now = np.where(D1 == h1)
                                break
                            
                        judge_knn1 = np.zeros((np.shape(t_now)[1], 4)).astype('int')
                        for m in range(0, np.shape(t_now)[1]):
                            judge_knn1[m, 0] = t_now[0][m]
                            judge_knn1[m, 1] = t_now[1][m]
                            judge_knn1[m, 2] = t_now[2][m]
    
                        judge_knn1[:, 0] = j - 5 + judge_knn1[:, 0] - 1
                        judge_knn1[:, 1] = 360 - 10 + judge_knn1[:, 1] - 1
    
                        for m in range(0, np.shape(t_now)[1]):
                            judge_knn1[m, 3] = mask[judge_knn1[m, 1], judge_knn1[m, 0], 
                                                        judge_knn1[m, 2]]
                        
                        if np.sum(judge_knn1[:, 3]) >= np.shape(t_now)[1]/2:
                            mask_3d[i-1, j-1, n] = 1
    
                    else:
                        mask_3d[i-1, j-1, n] = mask[i-1, j-1, n]
        # print(d)                       
                        
        for j in range(1, j1):
            
            if sst[i-1, j-1] == 0:
                
                datapath='/home/jingzhao/sundi/data/r_correlation_ecco2_real_trend/'
                data = scio.loadmat(datapath+'r_'+str(i)+'_'+str(j)+'.mat')
                
                r = data['r']
                r[np.isnan(r)] = 0
    
                for n in range(depth):
                    
                    t = np.where(r[j + 11*(i-350) + (n-1)*121, :] > beta)
                    
                    if len(t[0]) >= 1:
                        nt = len(t[0])
                        judge_knn = np.zeros((nt, 3)).astype('int')
    
                        [judge_knn[:, 0], judge_knn[:, 1], judge_knn[:, 2]] = \
                                         np.unravel_index(t[0], (11, 11, depth), order='F')
    
                        mask_now = np.zeros((11, 11, depth))
    
                        for h in range(np.shape(judge_knn)[0]):
                            mask_now[judge_knn[h, 0], judge_knn[h, 1], judge_knn[h, 2]] = 1
                        D1, N = cc3d.connected_components(mask_now, connectivity=26, return_N=True)
    
                        for h1 in range(1, N+1):
                            if D1[j-1, i-355+5, n] == h1:
                                t_now = np.where(D1 == h1)
                                break
                            
                        judge_knn1 = np.zeros(
                            (np.shape(t_now)[1], 4)).astype('int')
                        for m in range(0, np.shape(t_now)[1]):
                            judge_knn1[m, 0] = t_now[0][m]
                            judge_knn1[m, 1] = t_now[1][m]
                            judge_knn1[m, 2] = t_now[2][m]
    
                        judge_knn1[:, 1] = 360 - 10 + judge_knn1[:, 1] - 1
    
                        for m in range(0, np.shape(t_now)[1]):
                            judge_knn1[m, 3] = mask[judge_knn1[m, 1], judge_knn1[m, 0], 
                                                        judge_knn1[m, 2]]
                        
                        if np.sum(judge_knn1[:, 3]) >= np.shape(t_now)[1]/2:
                            mask_3d[i-1, j-1, n] = 1
    
                    else:
                        mask_3d[i-1, j-1, n] = mask[i-1, j-1, n]
                        
        
        
        for j in range(j2, 120+1):
            
            if sst[i-1, j-1] == 0:
                
                datapath='/home/jingzhao/sundi/data/r_correlation_ecco2_real_trend/'
                data = scio.loadmat(datapath+'r_'+str(i)+'_'+str(j)+'.mat')
                
                r = data['r']
                r[np.isnan(r)] = 0
    
                for n in range(depth):
                    
                    t = np.where(r[j-109 + 11*(i-350) + n*121 - 1, :] > beta)
                    
                    if len(t[0]) >= 1:
                        nt = len(t[0])
                        judge_knn = np.zeros((nt, 3)).astype('int')
    
                        [judge_knn[:, 0], judge_knn[:, 1], judge_knn[:, 2]] = \
                                         np.unravel_index(t[0], (11, 11, depth), order='F')
    
                        mask_now = np.zeros((11, 11, depth))
    
                        for h in range(np.shape(judge_knn)[0]):
                            mask_now[judge_knn[h, 0], judge_knn[h, 1], judge_knn[h, 2]] = 1
                        D1, N = cc3d.connected_components(mask_now, connectivity=26, return_N=True)
    
                        for h1 in range(1, N+1):
                            if D1[j-j2+6, i-355+5, n] == h1:
                                t_now = np.where(D1 == h1)
                                break
                            
                        judge_knn1 = np.zeros((np.shape(t_now)[1], 4)).astype('int')
                        for m in range(0, np.shape(t_now)[1]):
                            judge_knn1[m, 0] = t_now[0][m]
                            judge_knn1[m, 1] = t_now[1][m]
                            judge_knn1[m, 2] = t_now[2][m]
    
                        judge_knn1[:, 0] = 120 - 10 + judge_knn1[:, 0] - 1
                        judge_knn1[:, 1] = 360 - 10 + judge_knn1[:, 1] - 1
                        
                        for m in range(0, np.shape(t_now)[1]):
                            judge_knn1[m, 3] = mask[judge_knn1[m, 1], judge_knn1[m, 0], 
                                                        judge_knn1[m, 2]]
                        
                        if np.sum(judge_knn1[:, 3]) >= np.shape(t_now)[1]/2:
                            mask_3d[i-1, j-1, n] = 1
    
                    else:
                        mask_3d[i-1, j-1, n] = mask[i-1, j-1, n]
                    
    return mask_3d


#
#
# ##################################################################################


#
# MPI parallel
#


#
# call
# 

comm = MPI.COMM_WORLD

rank = comm.rank
size = comm.size

ave, res = divmod(len_yr, size)

count = [ave+1 if p < res else ave for p in range(size)]
displ = [sum(count[:p]) for p in range(size)]

var='Pres'

if rank == 0:
    data = h5py.File('/home/jingzhao/sundi/result/masks/mask_day_2020_1x1_60_200m_real_trend_verti_CFSR.mat', mode='r')
    mask_day = data['mask_day'][:]
    mask_day = mask_day.transpose(1, 3, 2, 0)
    s, p, pp, l = mask_day.shape
else:
    s, p, pp, l = 0, 0, 0, 0

# p for lon, pp for lat, l for depth and s for days

l = comm.bcast(l, root=0)
p = comm.bcast(p, root=0)
pp = comm.bcast(pp, root=0)
s = comm.bcast(s, root=0)

if rank != 0:
    mask_day = np.zeros((s,p*pp*l))

mask_ = np.zeros((count[rank], p*pp*l))
                                       
comm.Scatterv([mask_day.flatten(), np.array(count)*p*pp*l, np.array(displ)*p*pp*l, MPI.DOUBLE], mask_, root=0)
            
for i in range(count[rank]):

    mask = mask_[i].reshape(p,pp,l)
    tempdata = maximum_correlation_cluster(rank, sst, mask)

    scio.savemat('./ecco2_real_trend/2020/test%03d.mat'%(i+1+displ[rank]),{'mask':tempdata})
            
            
            
            
            
            
            
