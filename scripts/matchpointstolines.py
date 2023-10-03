#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 18 10:57:19 2021

@author: akh
"""
import pandas as pd
from htMOD import AKH_HT as HT
from sklearn.preprocessing import scale
import numpy as np 
from clusterMod import *
import matplotlib.pyplot as plt
from plotmod import plotlines, labelcolors, plotbyAngle, BA_HT, HThist
from examineMod import examineClusters
import seaborn as sns
from fitRectangle import allpoints
from pyproj import Proj

dfSegments=pd.read_csv('/home/akh/myprojects/Linking-and-Clustering-Dikes/dikedata/dikeset_ptheta.csv')

theta, rho, xc, yc= HT(dfSegments)
dfSegments['rho']=rho
dfSegments['theta']=theta

xc,yc=HT_center(dfSegments) 
dfLinked=pd.read_csv('/home/akh/myprojects/Linking-and-Clustering-Dikes/CRB_linked.csv')

points=pd.read_csv('/home/akh/myprojects/Linking-and-Clustering-Dikes/Scott Hughes Dike Locations.csv')
#points=pd.read_csv('/home/akh/myprojects/Linking-and-Clustering-Dikes/Hughes Data Classified 1.csv')
points.sort_values(by='FormationID')
npoints=len(points)

myProj = Proj("+proj=utm +zone=11, +ellps=WGS84 +datum=WGS84 +units=m")
#myProj = Proj("+proj=utm +zone=11, +ellps=WGS84 +datum=NAD83 +units=m")
#lon1, lat1 = myProj(df['Xstart'].values, df['Ystart'].values, inverse = True)
#lon2, lat2 = myProj(df['Xend'].values, df['Yend'].values, inverse = True)

linelocSegments=np.empty([npoints,4])
fig1,ax1=plt.subplots(2)
#ax[0].scatter(points['longitude'], points['latitude'], c='r')
#plotlines(df, 'k', ax[0], myProj=myProj)

UTMx,UTMy=myProj(points['Long.'],points['Lat.'])
ax1[0].scatter(UTMx,UTMy, c='r')
plotlines(dfSegments, 'k', ax1[0])

x1=dfSegments['Xstart'].values
x2=dfSegments['Xend'].values
y1=dfSegments['Ystart'].values
y2=dfSegments['Yend'].values
b=np.linspace(0,50,10)
    
for i in range(npoints):
    #print(i)
    dist=abs((x2-x1)*(y1-UTMy[i])-(x1-UTMx[i])*(y2-y1))/np.sqrt((x2-x1)**2+(y2-y1)**2)
    loc=np.argmin(dist)
    
    linelocSegments[i,0]=loc
    linelocSegments[i,1]=np.min(dist)
    linelocSegments[i,2]=dfSegments['theta'].iloc[loc]
    linelocSegments[i,3]=dfSegments['rho'].iloc[loc]
    
h1=ax1[1].hist(linelocSegments[:,1],bins=b)
ax1[1].set_xlabel("Distance to Closest Dike (m)")
ax1[1].set_ylabel('Points')
print(np.sum(linelocSegments[:,1]<10))

formations=(np.unique(points['FormationID']))

names=['SaddleMtns', 'Wanapum', 'Grande Ronde', 'Imnaha', 'PG']


f1,a2=plt.subplots(1,len(formations), sharey=True)
print(formations, "formations")
c,c2=labelcolors(points['FormationID'])
for i in range(len(formations)):
    print(i)
    mask=(points['FormationID']==formations[i])
    dikes=linelocSegments[mask,2]
    c=a2[i].hist(dikes, color=c2[i])
    a2[i].set_title(names[i])
    a2[i].set_xlim([-90,90])
    a2[i].set_xlabel("Angle (degrees)")
a2[0].set_ylabel("Counts")


f1,a2=plt.subplots(1,len(formations), sharey=True)
print(formations, "formations")
for i in range(len(formations)):
    print(i)
    mask=(points['FormationID']==formations[i])
    dikes=linelocSegments[mask,3]
    c=a2[i].hist(dikes, color=c2[i])
    a2[i].set_title(names[i])
    a2[i].set_xlim([min(linelocSegments[:,3]),max(linelocSegments[:,3])])
    a2[i].set_xlabel("Rho (m)")
a2[0].set_ylabel("Counts")

# fig2,ax2=plt.subplots(2)

# plotlines(dfLinked, 'k', ax2[0])
# ax2[0].scatter(UTMx,UTMy, c='r')
# x1=dfLinked['Xstart'].values
# x2=dfLinked['Xend'].values
# y1=dfLinked['Ystart'].values
# y2=dfLinked['Yend'].values
# linelocDikes=np.empty([npoints,2])
    
# for i in range(npoints):
#     #print(i)
#     dist=abs((x2-x1)*(y1-UTMy[i])-(x1-UTMx[i])*(y2-y1))/np.sqrt((x2-x1)**2+(y2-y1)**2)
#     loc=np.argmin(dist)
#     linelocDikes[i,0]=loc
#     linelocDikes[i,1]=np.min(dist)
    
# h2=ax2[1].hist(linelocDikes[:,1],bins=b)
# ax2[1].set_xlabel("Distance to Closest Dike (m)")
# ax2[1].set_ylabel('Points')
# print(np.sum(linelocDikes[:,1]<10))

# LinkedBetter=((linelocDikes[:,1]<10) & ~(linelocSegments[:,1]<10))
# fig3,ax3=plt.subplots()

# ax3.scatter(UTMx[LinkedBetter],UTMy[LinkedBetter], c='r')
# plotlines(dfSegments.iloc[linelocSegments[(linelocSegments[:,1]<10),0]], 'k', ax3)
# plotlines(dfLinked.iloc[linelocDikes[LinkedBetter,0]], 'g', ax3)