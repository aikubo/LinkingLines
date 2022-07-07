#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 12 16:38:08 2022

@author: akh
"""

import pandas as pd 
import matplotlib.pyplot as plt 
import numpy as np 
from htMOD import AKH_HT as HT
from htMOD import MidtoPerpDistance
from clusterMod import HT_AGG_custom as AggHT
from clusterMod import HT_AGG_custom2
from examineMod import examineClusters
from plotmod import *
import scipy.cluster.hierarchy as sch
from PrePostProcess import * 

dikeset=pd.read_csv('/home/akh/myprojects/Linking-and-Clustering-Dikes/dikedata/deccandata/deccan_central_cutouff_2331775y_processed.csv')
xc,yc=HT_center(dikeset)
# theta,rho,xc,yc=HT(dikeset)
# dikeset['theta']=theta
# dikeset['rho']=rho
dikeset=MidtoPerpDistance(dikeset, xc, yc)

dtheta=3 
drho=500
dikeset, Z=HT_AGG_custom2(dikeset, dtheta, drho)

lines,evaluation=examineClusters(dikeset)
lines.to_csv('/home/akh/myprojects/Linking-and-Clustering-Dikes/dikedata/deccandata/lines_deccan_central_cutoff.csv')
#np.save("/home/akh/myprojects/Linking-and-Clustering-Dikes/dikedata/deccandata/deccan_central_distanceMatrix", Z)

Z=np.load("/home/akh/myprojects/Linking-and-Clustering-Dikes/dikedata/deccandata/deccan_central_distanceMatrix.npy")
lines=pd.read_csv('/home/akh/myprojects/Linking-and-Clustering-Dikes/dikedata/deccandata/lines_deccan_central_cutoff.csv')


fig, ax, p, b=persistancePlot(Z)

fig,ax=plt.subplots(2,2)

DotsHT(fig, ax[0,0], lines, ColorBy='PerpOffsetDist')
ax[0,1].hist(lines['PerpOffsetDist'], bins=40)
ax[0,1].set_xlabel("PerpOffsetDist (m)")
ax[0,1].set_ylabel("Counts")
psplit=2175000#np.mean(lines['PerpOffsetDist'].values)

lines1=lines[ lines['PerpOffsetDist'] > psplit]
lines2=lines[ lines['PerpOffsetDist'] < psplit]

DotsHT(fig, ax[1,0], lines1, ColorBy='PerpOffsetDist')
DotsHT(fig, ax[1,1], lines2, ColorBy='PerpOffsetDist')

#fig, ax=HT3D(lines)

fig,ax=plt.subplots()
plotlines(lines, 'k', ax, ColorBy='PerpOffsetDist', alpha=0.6)


