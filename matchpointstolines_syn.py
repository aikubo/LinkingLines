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

def fragmentDikes(df, maxL=20000, ndikesMax=None, distortion=0):
    np.random.seed(5)
    dfFrag=pd.DataFrame(columns=df.columns)
    ndikes=0
    if ndikesMax is None: 
        ndikesMax=2000
        
       
    for i in range(len(df)):
        nSegments=np.random.randint(3)
        if ndikes < ndikesMax: 
            for j in range(nSegments):
                high=max(0, df['Xend'].iloc[i])
                low=min(0, df['Xend'].iloc[i])
                xrange=np.random.randint(low,high=high, size=2)
                m=df['Slope'].iloc[i]*(1+np.random.rand()*distortion)
                yrange=m*xrange
                L=np.sqrt((xrange[0]-xrange[1])**2+(yrange[0]-yrange[1])**2)
                if L > maxL: 
                    continue
                
                if max(abs(yrange)) > df['Yend'].max() :
                    continue 
                dfFrag=dfFrag.append(pd.DataFrame({'Xstart':xrange[0], 'Xend': xrange[1], 
                                                   'Ystart': yrange[0], 'Yend':yrange[1], 
                                                   'Slope':m, 'Length':L
                                                       }, index=[0]), ignore_index=True)
                ndikes=ndikes+1
                
    # a=np.full(len(dfFrag), False)
    # a[:int(len(dfFrag)/mask)]=True 
    # np.random.shuffle(a)
    # dfFrag=dfFrag.iloc[a]
    
    
    return dfFrag 

dikelength=[50000, -50000]

center=np.array([0,0])
ndikes=50
angles=np.linspace(-90,90,ndikes)
m=np.tan(angles)
Xstart=np.ones(ndikes*2)*0
Ystart=np.ones(ndikes*2)*0
Xend1=dikelength[0]/np.sqrt(1+m**2)
Xend=np.append(Xend1, dikelength[1]/np.sqrt(1+m**2))
m2=np.append(m,m)
Yend=m2*Xend
labels=np.arange(0,ndikes*2)

df=pd.DataFrame({"label":labels,'Xstart':Xstart, 'Xend': Xend, 'Ystart': Ystart, 'Yend':Yend, 'Slope':m2})

npoints=70
points=np.random.rand(npoints,2)*dikelength[0]/2
lineloc=np.empty([npoints])
fig,ax=plt.subplots()
ax.scatter(points[:,0], points[:,1], c='r')

plotlines(df, 'k', ax)
x1=df['Xstart'].values
x2=df['Xend'].values
y1=df['Ystart'].values
y2=df['Yend'].values

    
for i in range(npoints):
    dist=abs((x2-x1)*(y1-points[i,1])-(x1-points[i,0])*(y2-y1))/np.sqrt((x2-x1)**2+(y2-y1)**2)
    loc=np.argmin(dist)
    print(labels[loc])