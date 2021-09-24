#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 17 11:23:39 2021

@author: akh
"""
import pandas as pd
from htMOD import AKH_HT as HT
from sklearn.preprocessing import scale
import numpy as np 
from clusterMod import *
import matplotlib.pyplot as plt
from plotmod import plotlines, labelcolors, plotbyAngle, BA_HT
from examineMod import examineClusters, plotlabel
from PrePostProcess import *

from findRadialCenters import sweepCenters, detectLocalPeaks
'''
dikeset=pd.read_csv('/home/akh/myprojects/Linking-and-Clustering-Dikes/dikedata/dikeset_ptheta.csv')
dikeset=giveID(dikeset)
theta, rho, xc, yc= HT(dikeset)
dikeset['rho']=rho
dikeset['theta']=theta

trange=2 
rrange=5000 

stdT=np.std(theta)
stdRho=np.std(rho)
d2=trange/stdT

clustering=HT_AGG(dikeset,d2)
dikeset['Labels']=clustering.labels_
lines,IC=examineClusters(dikeset)
dikeset.to_csv('/home/akh/myprojects/Linking-and-Clustering-Dikes/dikedata/dikeset_ptheta.csv',index=False)

#generate a kmz
#colorsSegments=labelcolors(dikeset['Labels'])
#colorsDikes=labelcolors(lines['Label'])

#errorAnalysis(lines)

#lines2= writeToQGIS(lines, "CRBLinkedQGIS.csv")

# label=22
# plotlabel(dikeset,label)
#fig, ax= DotsHT(dikeset,lines)
'''
xc=475201.3737670694
yc=4976341.597868606
lines=pd.read_csv('/home/akh/myprojects/Linking-and-Clustering-Dikes/dikedata/CRBLinked0621.csv')
#lines.to_csv('/home/akh/myprojects/Linking-and-Clustering-Dikes/dikedata/CRBLinked0621.csv')

err, dikes,thetaRange, xs, ys=sweepCenters(lines, 5000, 5000, xc,yc)
lines,labels,counts=detectLocalPeaks(err,dikes, lines, plot=True)


col, colshort=labelcolors(lines["MainCatagory"], cm.turbo)

fig,ax=plt.subplots()
plotlines(lines, col, ax[0])

# custom_lines = [Line2D([0], [0], color=colshort[0], lw=4),
#         Line2D([0], [0], color=colshort[1], lw=4),
#         Line2D([0], [0], color=colshort[2], lw=4)]


# fig,ax=plt.subplots(1,4)
# #plotlines(Splines, col, ax[0])
# ax[0].legend(custom_lines, ["Linear", "SpanishPeaks", "DikeMountain"])

# for i in np.unique(labels):
#     angles=Splines[Splines['MainCatagory']==i]['AvgTheta']
#     astd=np.std(angles)
#     sortedAngles=np.sort(angles)
#     AverageAngularSpacing=np.mean(np.abs(sortedAngles[0:-1]-sortedAngles[1:]))
#     avgangle=np.mean(angles)
#     ax[1].scatter(avgangle, astd, c=colshort[i])
#     ax[2].scatter(astd, AverageAngularSpacing,c=colshort[i])
#     print("Label:", i)
#     print("Angular Spacing:", AverageAngularSpacing)
    

#TopHTSection(lines, 5000, 3)
#fig,ax=plotbyAngle(dikeset, lines, 20)

## Length and Width plots
# f,a=plt.subplots(2)
# a[0].hist(lines['R_Length'], bins=np.arange(0,100000,5000))
# a[0].set_xlabel('Dike Cluster Length')

# a[1].hist(lines['R_Width'],bins=np.arange(0,2000,100))
# a[1].set_xlabel('Dike Cluster Width')
# plt.tight_layout()

# f2,a2=plt.subplots(3)
# a2[0].scatter(lines['R_Length'], lines['R_Width'])
# a2[0].set_xlabel('Length')
# a2[0].set_ylabel('Width')

# a2[2].scatter(lines['Size'], lines['R_Width'])
# a2[2].set_xlabel('Size')
# a2[2].set_ylabel('Width')

# a2[1].scatter(lines['Size'], lines['R_Length'])
# a2[1].set_xlabel('Size')
# a2[1].set_ylabel('Length')
# fig,ax, h1, h2=BA_HT(dikeset, lines, rstep=2000)

plt.tight_layout()


