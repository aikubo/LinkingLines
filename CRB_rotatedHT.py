#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 19 12:46:26 2021

@author: akh
"""
import pandas as pd
from htMOD import AKH_HT
from htMOD import rotateData2
from sklearn.preprocessing import scale
import numpy as np 
from clusterMod import *
import matplotlib.pyplot as plt
from plotmod import plotlines, labelcolors, plotbyAngle, BA_HT
from examineMod import examineClusters
from jitteringHTcenter import moveHTcenter

regular=pd.read_csv('/home/akh/myprojects/Linking-and-Clustering-Dikes/dikedata/dikeset_ptheta.csv')

linked=pd.read_csv('/home/akh/myprojects/Linking-and-Clustering-Dikes/dikedata/CRBLinked0621.csv')
dikeset=rotateData2(regular, 30)

theta, rho, xc, yc= AKH_HT(dikeset)
dikeset['rho']=rho
dikeset['theta']=theta

#fig,ax=plt.subplots()
#plotlines(dikeset1,'k', ax, center=True)
#plotlines(dikeset,'r', ax, center=True)
trange=2 
rrange=5000 

stdT=np.std(theta)
stdRho=np.std(rho)
d2=trange/stdT

clustering=HT_AGG(dikeset,d2)
dikeset['Labels']=clustering.labels_
dikeset['p']=rho
dikeset['theta']=theta
lines,IC=examineClusters(dikeset)

#generate a kmz
colorsSegments=labelcolors(dikeset['Labels'])
colorsDikes=labelcolors(lines['Label'])
errorAnalysis(lines)
errorAnalysis(linked)

label=22
plotlabel(dikeset,label)
label=22
plotlabel(regular,label)

#fig,ax=plotbyAngle(dikeset, lines, 20)

#)# Length and Width plots
# f,a=plt.subplots(2)
# a[0].hist(Deccanlines['R_Length'], bins=np.arange(0,100000,5000))
# a[0].set_xlabel('Dike Cluster Length')

# a[1].hist(Deccanlines['R_Width'],bins=np.arange(0,2000,100))
# a[1].set_xlabel('Dike Cluster Width')
# plt.tight_layout()

# f2,a2=plt.subplots(3)
# a2[0].scatter(Deccanlines['R_Length'], Deccanlines['R_Width'])
# a2[0].set_xlabel('Length')
# a2[0].set_ylabel('Width')

# a2[2].scatter(Deccanlines['Size'], Deccanlines['R_Width'])
# a2[2].set_xlabel('Size')
# a2[2].set_ylabel('Width')

# a2[1].scatter(Deccanlines['Size'], Deccanlines['R_Length'])
# a2[1].set_xlabel('Size')
# a2[1].set_ylabel('Length')