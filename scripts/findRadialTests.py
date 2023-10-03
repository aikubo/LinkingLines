#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov  1 16:25:42 2021

@author: akh
"""
import matplotlib.pyplot as plt
from plotmod import plotlines, labelcolors, plotbyAngle, BA_HT
from examineMod import examineClusters, plotlabel
from PrePostProcess import *
from matplotlib import cm
from findRadialCenters import sweepCenters, detectLocalPeaks
from matplotlib.lines import Line2D  

from matplotlib import cm
from skimage.morphology import reconstruction

from scipy.ndimage import gaussian_filter

from skimage import img_as_float
from skimage.filters import threshold_otsu, threshold_local
import numpy as np
import matplotlib.pyplot as plt
from scipy import ndimage as ndi
from plotmod import plotlines
from skimage.segmentation import watershed
from skimage.feature import peak_local_max

from skimage import feature
from htMOD import HT_center, AKH_HT
from findRadialCenters import sweepCenters, detectLocalPeaks

lines=pd.read_csv('/home/akh/myprojects/Linking-and-Clustering-Dikes/syntheticRadial_testfindRadial_doubledFalse.csv')
theta, rho, xc, yc= AKH_HT(lines)
lines['AvgRho']=rho
lines['AvgTheta']=theta

grids=[1000,2000, 3000, 5000, 7000]
thres=[500,1000,2000, 3000, 5000, 7000, 10000]

gridsArr,thresArr=np.meshgrid(grids,thres)

fig,ax=plt.subplots(1,2)
plotlines(lines, 'k', ax[0], center=True)
ax[1].scatter(theta,rho, c=lines['label'], cmap=cm.turbo)
foundArr=np.empty_like(gridsArr)
xc,yc=HT_center(lines)
trueCenters=np.array( [0, 40000,-100000,200000, 10000, 0,40000,100000 ,-300000,500000]).reshape(2,5)

for i in range(len(grids)):
    for j in range(len(thres)):
        
        err, dikes, thetaRange, thetaSTD, xs, ys=sweepCenters(lines, grids[i], thres[j], xc,yc)
        lines,labels,counts, centerData=detectLocalPeaks(err,dikes, lines,grids[i])
        print("Found", len(np.unique(labels))-1, "radial centers")
    
    
        found=0
        for k in np.unique(labels): 
            if k == 0:
                continue 
            c=np.array([centerData['CenterX'].iloc[k], centerData['CenterY'].iloc[k]]).reshape(2,1)
            centerError=np.sqrt(np.einsum('ij,ij->j', trueCenters-c, trueCenters-c))
            print(k, np.argmin(centerError), "error", np.min(centerError) )
            
            if np.min(centerError) < 1000:
                print( "Found True Radial Center")
                found=found+1 
                
        print("Found", found, "out of 5")
        foundArr[j,i]=found
fig,ax=plt.subplots()
c=ax.pcolor(gridsArr, thresArr, foundArr, shading='auto')
plt.colorbar(c, label="Radial Centers out of 5")
ax.set_xlabel("Grid Spacing (m)")
ax.set_ylabel("Threshold (m)")