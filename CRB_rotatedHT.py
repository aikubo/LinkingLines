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
from clusterMod import scalesimple
regular=pd.read_csv('/home/akh/myprojects/Linking-and-Clustering-Dikes/dikedata/dikeset_ptheta.csv')

CRBl=pd.read_csv('/home/akh/myprojects/Linking-and-Clustering-Dikes/dikedata/CRBLinked0621.csv')
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
fig,ax=plt.subplots(2,4)
d1=2/np.std(regular['theta'])
clustering=HT_AGG(dikeset,d1)
dikeset['Labels']=clustering.labels_
dikeset['p']=rho
dikeset['theta']=theta
lines,IC=examineClusters(dikeset)

# ax[0,0].hist(CRBl['AvgTheta'])
# ax[0,1].hist(CRBl['AvgRho'])
# ax[0,2].scatter(CRBl['AvgTheta'], CRBl['AvgRho'])

# ax[1,0].hist(lines['AvgTheta'])
# ax[1,1].hist(lines['AvgRho'])

# dist=np.sqrt( (CRBl['AvgTheta']- lines['AvgTheta'])**2 + (CRBl['AvgRho']- lines['AvgRho'])**2)
# ax[1,2].scatter(lines['AvgTheta'], lines['AvgRho'], c=dist)

# # #generate a kmz
# colorsSegments=labelcolors(dikeset['Labels'],cm.turbo)
# colorsDikes=labelcolors(lines['Label'],cm.turbo)
errorAnalysis(lines, dikeset)
errorAnalysis(CRBl, regular)

#errorAnalysis(linked)

#c, err=checkClusterChange(CRBl, lines)


f,a=plt.subplots(3)
X=scalesimple(regular['rho'], regular['theta'])
a[0].hist(X[:,0])


X2=scalesimple(dikeset['rho'], dikeset['theta'])
a[2].hist(X2[:,0])
print(np.sum(X==X2))


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