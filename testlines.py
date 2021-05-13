#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 10 13:00:55 2020

@author: akh
"""
import numpy as np 
from plotmod import plotlines
import pandas as pd
from htMOD import AKH_HT as HT 
import matplotlib.pyplot as plt
from sklearn.cluster import DBSCAN
import matplotlib.lines as mlines
from sklearn.preprocessing import scale

true=pd.read_csv('/home/akh/myprojects/Linking-and-Clustering-Dikes/test_rand1.csv')
theta, r, xc, yc= HT(true, xc=0, yc=0)


fig,ax=plt.subplots(2)

plotlines(true,'r', ax[0], linewidth=3)

ax[1].scatter(theta,r,c='red', edgecolor="black", s=100)
ax[1].set_ylabel('Rho (m)')
ax[1].set_xlabel('Theta (deg)')

true=pd.read_csv('/home/akh/myprojects/Linking-and-Clustering-Dikes/test_rand2.csv')
theta, r, xc, yc= HT(true, xc=0, yc=0)

plotlines(true,'b', ax[0], linewidth=3)

ax[1].scatter(theta,r, c="blue", edgecolor="black", s=100)
ax[1].set_ylabel('Rho (m)')
ax[1].set_xlabel('Theta (deg)')


true=pd.read_csv('/home/akh/myprojects/Linking-and-Clustering-Dikes/test_rand3.csv')
theta, r, xc, yc= HT(true, xc=0, yc=0)

plotlines(true,'g', ax[0], linewidth=3)

ax[1].scatter(theta,r, c="green", edgecolor="black", s=100)
ax[1].set_ylabel('Rho (m)')
ax[1].set_xlabel('Theta (deg)')
ax[1].set_xlim([-90,90])
ax[0].plot(0,0, 'w*', mec='black')


# fig,ax=plt.subplots()
# lines=pd.DataFrame(columns=['Xstart', 'Ystart', 'Xend', 'Yend', 'p', 'theta'])
# for i in range(0,len(nclusters)): 
#     #print(i)
#     mask=[clustering.labels_==i][0]
   
#     rrange=max(r[mask])-min(r[mask])
#     trange=max(theta[mask])-min(theta[mask])
#     avgrho=np.average(r[mask])
#     avgtheta=np.average(theta[mask])
#     xstart_all=np.array(true['Xstart'].loc[mask])
#     xend_all=np.array(true['Xend'].loc[mask])
#     ystart_all=np.array(true['Ystart'].loc[mask])
#     yend_all=np.array(true['Yend'].loc[mask])
    
    
    
#     xs=np.concatenate((xstart_all, xend_all), axis=0)
#     ys=np.concatenate((ystart_all,yend_all), axis=0)
#     ystart=max(ys)
#     xstart=max(xs)
    
#     xend=min(xs)
#     yend=min(ys)
#     xmid=(xstart+xend)/2
#     ymid=(ystart+yend)/2
#     print(xmid, ymid)
#     l=np.sqrt( (xstart-xend)**2+(ystart-yend)**2)
    
#     a = np.cos(np.deg2rad(avgtheta))
#     b = np.sin(np.deg2rad(avgtheta))
    
#     x0 = (a * avgrho) + xc
#     y0 = (b * avgrho) + yc
    
#     print(x0,y0)
#     x1 = xmid - l/2 * (-b)
#     y1 = ymid - l/2 * (a)
#     x2 = xmid + l/2 * (-b)
#     y2 = ymid + l/2 * (a)
#     ax.add_line(mlines.Line2D([x1, x2], [y1, y2], color='k'))
    
#     print([x1, x2], [y1, y2])
    
# plotlines(true, 'r', ax)