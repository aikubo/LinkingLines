#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 10 13:45:32 2021

@author: akh
"""
import pandas as pd
from htMOD import AKH_HT as HT
from htMOD import transformXstart
from sklearn.preprocessing import scale
import numpy as np 
from clusterMod import *
import matplotlib.pyplot as plt
from plotmod import plotlines, labelcolors, plotbyAngle, BA_HT, DotsHT
from examineMod import examineClusters, plotlabel
from PrePostProcess import *
from fitRectangle import endpoints2
from matplotlib import cm
import matplotlib.colors as colors
from matplotlib.lines import Line2D      
from plotmod import plotlines
import seaborn as sns

sns.set_context("talk")
# dikeset=pd.read_csv('/home/akh/myprojects/Linking-and-Clustering-Dikes/dikedata/DikeMountain_Peaks_3857WKT.csv')
# #dikeset=pd.read_csv('/home/akh/myprojects/Linking-and-Clustering-Dikes/dikedata/Peaks_3857.csv')
# dikeset=WKTtoArray(dikeset)
# dikeset=giveID(dikeset)
# 
# dikeset['rho']=rho
# dikeset['theta']=theta

dikeset=pd.read_csv('/home/akh/myprojects/Linking-and-Clustering-Dikes/dikedata/DikeMountain_SpanishPeaks_3857.csv')
theta, rho, xc, yc= HT(dikeset)
trange=2 
rrange=5000 

stdT=np.std(theta)
stdRho=np.std(rho)
d2=trange/stdT

clustering=HT_AGG(dikeset,d2)
dikeset['Labels']=clustering.labels_
dikeset['p']=rho
dikeset['theta']=theta
Splines,IC=examineClusters(dikeset)

Splines=transformXstart(Splines)

#fig, ax= DotsHT(Splines, ColorBy="Xstart")

from findRadialCenters import sweepCenters, detectLocalPeaks

err, dikes, thetaRange, thetaSTD, xs, ys=sweepCenters(Splines, 1000, 1000, xc,yc)
Splines,labels,counts, centerData=detectLocalPeaks(err,dikes, Splines,1000)


col, colshort=labelcolors(Splines["MainCatagory"], cm.viridis)
# custom_lines = [Line2D([0], [0], color=colshort[0], lw=4),
#         Line2D([0], [0], color=colshort[1], lw=4),
#         Line2D([0], [0], color=colshort[2], lw=4)]
'''

#fffc2b
#ce11b8
#3c72ce

'''
col=[ '#fffc2b', '#ce11b8', '#3c72ce']

fig,ax=plt.subplots()
#plotlines(Splines, col, ax[0])
#ax[0].legend(custom_lines, ["Linear", "SpanishPeaks", "DikeMountain"])
ax.set_ylabel('Rho (m)')
ax.set_xlabel('Theta ($^\circ$)')


fig2,ax1=plt.subplots(3,2)
#plotlines(Splines, col, ax[0])
#ax[0].legend(custom_lines, ["Linear", "SpanishPeaks", "DikeMountain"])
ax.set_ylabel('Rho (m)')
ax.set_xlabel('Theta ($^\circ$)')

ax1[0][1].set_title("Rho (m)")
ax1[0][0].set_title("Theta ($^\circ$)")
for i in np.unique(labels):
    mask=Splines['MainCatagory']==i
    angles=Splines[Splines['MainCatagory']==i]['AvgTheta']
    a=ax1[i][0]
    a2=ax1[i][1]
    a.hist(angles, bins=np.arange(-90,90,4), color=col[i], alpha=0.6)
    a2.hist( Splines['AvgRho'].loc[mask], color=col[i], alpha=0.6, bins=np.arange(min(rho), max(rho), 1000))
    
for i in np.unique(labels):
    a=ax1[i][0]
    a2=ax1[i][1]
    mask=Splines['MainCatagory']==i
    angles=Splines[Splines['MainCatagory']==i]['AvgTheta']
    astd=np.std(angles)
    sortedAngles=np.sort(angles)
    AverageAngularSpacing=np.mean(np.abs(sortedAngles[0:-1]-sortedAngles[1:]))
    avgangle=np.mean(angles)
    ax.scatter(Splines['AvgTheta'].loc[mask], Splines['AvgRho'].loc[mask], c=col[i])
    a.hist(angles, bins=np.arange(-90,90,4), color=col[i], histtype="step")
    a2.hist( Splines['AvgRho'].loc[mask], color=col[i],histtype="step", bins=np.arange(min(rho), max(rho), 1000) )
    #ax[2].scatter(astd, AverageAngularSpacing,c=colshort[i])
    print("Label:", i)
    print("Dikes:", np.sum(Splines['MainCatagory']==i))
    print( "Average Angle:" , avgangle)
    print("Angular Spacing:", AverageAngularSpacing)
    
plt.tight_layout()

# import imagepers

    
    

#writeToQGIS(Splines, "CatagorizedSpanishPeaks_10032021.csv")


# g0=imagepers.persistence(err)
# print(g0)
# fig = plt.figure()
# ax = fig.add_subplot(111)
# ax.set_title("Peristence diagram")
# ax.plot([0,100], [0,100], '-', c='grey')
# for i, homclass in enumerate(g0):
#     p_birth, bl, pers, p_death = homclass
#     if pers <= 1.0:
#         continue
    
#     x, y = bl, bl-pers
#     ax.plot([x], [y], '.', c='b')
#     ax.text(x, y+2, str(i+1), color='b')
# ax.set_xlabel("Birth level")
# ax.set_ylabel("Death level")

# im=err

# fig,ax=plt.subplots(1,2)
# ax[0].set_title("Loci of births")
# for i, homclass in enumerate(g0):
#     p_birth, bl, pers, p_death = homclass
#     if pers <= 20.0:
#         continue
#     y, x = p_birth
#     ax[0].plot(xs[x], ys[y], '.', c='b')
#     ax[0].text(xs[x], ys[y]+5000, str(i+1), color='b')
    
# ax[0].pcolor(xs, ys, err, cmap=cm.Reds, shading='auto')
# ax[1].set_xlim((0,im.shape[1]))
# ax[1].set_ylim((0,im.shape[0]))
# plt.gca().invert_yaxis()

Splines.to_csv('/home/akh/myprojects/Linking-and-Clustering-Dikes/dikedata/SpanishSilverLinked0621.csv')