#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 22 11:52:19 2022

@author: akh
"""
from PrePostProcess import * 
import pandas as pd 
import numpy as np 
from htMOD import AKH_HT as HT 
from htMOD import HT_center, MidtoPerpDistance, moveHTcenter
from plotmod import *
import matplotlib.pyplot as plt 
from clusterMod import HT_AGG_custom
import scipy.cluster.hierarchy as sch
from examineMod import examineClusters, checkoutCluster, TopHTSection
import os
import labellines
from matplotlib.patches import Rectangle

import seaborn as sns

d=['dikedata/deccandata/All_deccan_3857_preprocessed.csv',
   'dikedata/crb/CJDS_FebStraightened.csv',
   'dikedata/spanish peaks/SpanishPeaks_3857_preprocessed.csv'
   ]

l=['dikedata/deccandata/AllDeccan_3_2000_LINKED.csv',
   'dikedata/crb/CJDS_Lines_3_500_March11.csv',
   'dikedata/spanish peaks/SpanishPeaks_3857_LINKED_2_2000.csv']

c=['r', 'g', 'b']
moveX=[0, -7.0e5, -1.0e6]
moveY=[0, 3.5e5, 5.24e5]

fig,ax=plt.subplots()
fig.set_size_inches( 12,6)

fig2,ax2=plt.subplots()
ax2.set_ylabel('Swarm Area (km $^{2}$)')
ax2.set_xlabel("Number of Dikes")
ax2.set_yscale('log')
averagewidth=[1,8,8]
name=['Deccan', 'CJDS', 'Spanish Peaks']
dfAll=pd.DataFrame()

for dikeset,lines, color, xnew, ynew, w, n in zip(d,l,c, moveX, moveY, averagewidth, name): 
    
    df=pd.read_csv(dikeset)
    
    

    xc,yc=HT_center(df)
    df2=moveHTcenter(df)
    df2=moveHTcenter(df2, xc=xnew, yc=ynew)
    df3=pd.read_csv(lines)
    df3['SwarmID']=n
        
    df3['R_Length'] = pd.to_numeric(df3['R_Length'], errors='coerce')
    df3['R_Width'] = pd.to_numeric(df3['R_Width'], errors='coerce')
    df3['Size'] = pd.to_numeric(df3['Size'], errors='coerce')
    
    xs=[ min(df2.min().iloc[ [0,2]]), max(df2.max().iloc[ [0,2]])]
    ys=[ min(df2.min().iloc[ [1,3]]), max(df2.max().iloc[ [1,3]])]
    width=np.ptp(xs)
    height=np.ptp(ys)
    Xedges=np.array([xs[0], xs[0], xs[1], xs[1], xs[0]])
    Yedges=np.array([ys[1], ys[0], ys[0], ys[1], ys[1]])
    ax.plot(Xedges, Yedges, color, linewidth=4)
    rect=Rectangle( (xs[0], ys[0]), width, height,facecolor=color, alpha=0.3, label=n)
    ax.add_patch(rect)
    ax.legend(loc='upper right')
    FixCartesianLabels(ax)
  
    
    
    n=len(df)
    plotlines(df2, color, ax, linewidth=0.5, alpha=0.4)
    ax2.plot(n, width*height/1000**2, "*", color=color, label=name)
    ax2.plot(n, df['seg_length'].sum()*w/1000**2, "p", color=color )
    a=df3['R_Length'].values*df3['R_Width'].values
    a[np.isnan(a)]=0
    a=np.sum(a)/1000**2
    ax2.plot(len(df3), a, 'D', color=color)
    print(width*height/1000**2,  df['seg_length'].sum()*w/1000**2 ,a)
    
    df3=df3.loc[df3['Size']>4]
    
    a=df3['R_Length'].values*df3['R_Width'].values
    a[np.isnan(a)]=0
    a=np.sum(a)/1000**2
    ax2.plot(len(df3), a, 'v', color=color)
    
    dfAll=pd.concat([dfAll, df3]).reset_index(drop=True)

dfAll['R_Length']=dfAll['R_Length']/1000

f, ax = plt.subplots(figsize=(7, 6))
sns.boxplot(x='SwarmID', y='R_Length', data=dfAll)
sns.despine(offset=5, trim=True)