#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 27 14:50:14 2021

@author: akh
"""
import pandas as pd
from htMOD import AKH_HT as HT
from sklearn.preprocessing import scale
import numpy as np 
from clusterMod import *
import matplotlib.pyplot as plt
from plotmod import plotlines, labelcolors, plotbyAngle, BA_HT
from examineMod import examineClusters, plotlabel,extendLines
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
from plotmod import plotlines, HT3D
from skimage.segmentation import watershed
from skimage.feature import peak_local_max

from skimage import feature
from htMOD import HT_center
from findRadialCenters import sweepCenters, detectLocalPeaks

lines=pd.read_csv("/home/akh/myprojects/Linking-and-Clustering-Dikes/AllCRBLinked.csv")
xc,yc=HT_center(lines)
gridM=1000

lines=MidtoPerpDistance(lines, xc, yc)

#X=(np.vstack(( np.deg2rad(theta), rho, lines['PerpOffsetDist'])).T)
midist=np.log(lines['PerpOffsetDist'].values)

HT3D(lines['AvgTheta'],lines['AvgRho'],midist, midist)

# err, dikes, thetaRange, thetaSTD, xs, ys=sweepCenters(lines, gridM, 5000, xc,yc, plot=True)
# lines,labels,counts, centerData=detectLocalPeaks(err,dikes, lines,gridM, plot=True)


# xc=475201.3737670694
# yc=4976341.597868606

dikeset=pd.read_csv("/home/akh/myprojects/Linking-and-Clustering-Dikes/dikedata/allCRB_dikes.csv")
dikeset=WKTtoArray(dikeset)

swarms=['CJDS', 'Ice_Harbor', 'Monument_Dikes',
        'Steens_dikes']

allCRB={}
fig,ax=plt.subplots(1,2)
xc,yc=HT_center(dikeset)
col=['#30123BFF', '#1AE4B6FF', '#FABA39FF', '#7A0403FF']
allLinked=pd.DataFrame()

b=0
for i in swarms: 
    
    df=dikeset[dikeset['layer']==i]

    theta, rho, xc, yc= HT(df, xc=xc, yc=yc)
    df['rho']=rho
    df['theta']=theta
    
    trange=2
    rrange=5000 
    
    stdT=np.std(theta)
    stdRho=np.std(rho)
    d2=trange/stdT
    dfromCJDS=0.07184050241647064
    d=d2
    clustering=HT_AGG(df,d)
    df['Labels']=clustering.labels_
    lines,IC=examineClusters(df)
    lines['Layer']=[i]*len(lines)
    allCRB[i]=df
    allCRB[i+"Linked"]=lines
    allCRB[i+"Center"]=[xc,yc]
    allLinked=allLinked.append(lines)
    plotlines(lines, col[b], ax[0])
    ax[1].scatter(theta, rho, color=col[b], alpha=0.6, edgecolor='k')
    #longlines=extendLines(lines, save=True, name=str(i)+'_extendedlines.csv')
    b=b+1

ax[1].set_ylabel("Rho (m)")
ax[1].set_xlabel("Theta ($^\circ$)")
custom_lines = [Line2D([0], [0], color=col[0], lw=4),
        Line2D([0], [0], color=col[1], lw=4),
        Line2D([0], [0], color=col[2], lw=4),
        Line2D([0], [0], color=col[3], lw=4)]

ax[0].plot(xc,yc, 'r*')
ax[1].legend(custom_lines, swarms)


