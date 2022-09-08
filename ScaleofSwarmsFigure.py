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

# d=['dikedata/deccandata/AllDeccan_PreProcessed.csv',
#    'dikedata/crb/allCRB_dikes_PreProcessed.csv',
#    'dikedata/spanish peaks/SpanishPeaks_3857_preprocessed.csv'
#    ]

# l=['dikedata/deccandata/AllDeccanLinked_24_08_2022.csv',
#    'dikedata/crb/AllCRBLinked_24_08_2022.csv',
#    'dikedata/spanish peaks/SpanishPeaks_Complete_3_2013.csv']

# c=['r', 'g', 'b']
# moveX=[0, -7.0e5, -1.0e6]
# moveY=[0, 3.5e5, 5.24e5]

# fig,ax=plt.subplots()
# fig.set_size_inches( 12,6)

# fig2,ax2=plt.subplots(3, sharex=True)
# ax2[0].set_ylabel('Dike Length (km)')
# #ax2[0].set_xlabel("Number of Dikes")
# ax2[0].set_yscale('log')
# ax2[0].set_xscale('log')

# ax2[1].set_ylabel('Swarm Area (km $^{2}$)')
# #ax2[1].set_xlabel("Number of Dikes")
# ax2[1].set_yscale('log')
# ax2[1].set_xscale('log')

# ax2[2].set_ylabel('Erupted Volume (km $^{3}$)')
# ax2[2].set_xlabel("Number of Dikes")
# ax2[2].set_yscale('log')
# ax2[2].set_xscale('log')




# averagewidth=[2,8,8]
# name=['Deccan', 'CJDS', 'Spanish Peaks']
# dfAll=pd.DataFrame() 
# EruptedVolume=[ 1.3e6, 210000, 50] #50,210000, 1.3e6] #km
# # CRBG estimate: Kasbohm 2018
# # Deccan estimate: Jay and Widdowson 2008
# # Spanish peaks erupted volume:
    
# TrustOn=False  # only use the clustered dikes that are trusted
# for dikeset,lines, color, xnew, ynew, w, label, ev in zip(d,l,c, moveX, moveY, averagewidth, name, EruptedVolume): 
    
#     df=pd.read_csv(dikeset)
#     #df=DikesetReProcess(df)
    

#     xc,yc=HT_center(df)
#     df2=moveHTcenter(df)
#     df2=moveHTcenter(df2, xc=xnew, yc=ynew)
#     df3=pd.read_csv(lines)
#     df3['SwarmID']=label
        
#     df3['R_Length'] = pd.to_numeric(df3['R_Length'], errors='coerce')
#     df3['R_Width'] = pd.to_numeric(df3['R_Width'], errors='coerce')
#     df3['Size'] = pd.to_numeric(df3['Size'], errors='coerce')
#     #m=df3['TrustFilter']==1
#     if TrustOn:
#         df3=df3[m]
    
#     xs=[ min(df2.min().iloc[ [0,2]]), max(df2.max().iloc[ [0,2]])]
#     ys=[ min(df2.min().iloc[ [1,3]]), max(df2.max().iloc[ [1,3]])]
#     width=np.ptp(xs)
#     height=np.ptp(ys)
#     Xedges=np.array([xs[0], xs[0], xs[1], xs[1], xs[0]])
#     Yedges=np.array([ys[1], ys[0], ys[0], ys[1], ys[1]])
#     ax.plot(Xedges, Yedges, color, linewidth=4)
#     rect=Rectangle( (xs[0], ys[0]), width, height,facecolor=color, alpha=0.3, label=label)
#     ax.add_patch(rect)
#     ax.legend(loc='upper right')
#     FixCartesianLabels(ax)
  

#     n=len(df)
#     #plotlines(df2, color, ax, linewidth=0.5, alpha=0.4)
    
#     #comparison plots
    
#     # segment Length
#     ax2[0].plot(n,  df['seg_length'].mean(), "s", color=color, label=label)
#     ax2[0].plot(n,  df3['R_Length'].mean(), "P", color=color, label=label)
    
#     # width*length of swarm area
#     ax2[1].plot(n,  width*height/1000**2, "*", color=color, label=label)
#     # linked dike length * width (not packet width) (only 4+segments)
    
#     l=df3['R_Length'].mean()
#     ax2[1].plot(n,  l*w/1000**2, "p", color=color, label=label)
    
#     #erupted volume 
#     ax2[2].plot(n,  ev, "*", color=color, label=label)
#     ax2[2].legend()
#     plt.tight_layout()
    
#     # ax2.plot(n, width*height/1000**2, "*", color=color, label=name)
#     # ax2.plot(n, df['seg_length'].sum()*w/1000**2, "p", color=color )
#     # a=df3['R_Length'].values*df3['R_Width'].values
#     # a[np.isnan(a)]=0
#     # a=np.sum(a)/1000**2
#     # ax2.plot(len(df3), a, 'D', color=color)
#     # print(width*height/1000**2,  df['seg_length'].sum()*w/1000**2 ,a)
    
#     # df3=df3.loc[df3['Size']>4]
    
#     # a=df3['R_Length'].values*df3['R_Width'].values
#     # a[np.isnan(a)]=0
#     # a=np.sum(a)/1000**2
#     # ax2.plot(len(df3), a, 'v', color=color)
    
#     dfAll=pd.concat([dfAll, df3]).reset_index(drop=True)
from matplotlib.gridspec import GridSpec

dikeset=pd.read_csv('/home/akh/myprojects/Linking-and-Clustering-Dikes/dikedata/AllDatasets25_08_2022.csv')
dfAll=pd.read_csv('/home/akh/myprojects/Linking-and-Clustering-Dikes/dikedata/AllDatasetsLinked25_08_2022.csv')
dfAll['R_Length']=dfAll['R_Length']/1000

dfAllTrusted=dfAll[ dfAll['TrustFilter']==1]
# f, ax = plt.subplots(2,1, figsize=(12, 6))
# ax[0].set_yscale('symlog')
# ax[1].set_yscale('symlog')
# sns.boxplot(x='Region', y='R_Length', data=dfAllTrusted, ax=ax[0])
# sns.boxplot(x='SwarmID', y='R_Length', data=dfAllTrusted, ax=ax[1])
# sns.despine(offset=5, trim=True)
sns.set_theme(style='ticks')
fig = plt.figure(figsize=(12,6), constrained_layout=True)

gs = GridSpec(1, 4, figure=fig)
ax1 = fig.add_subplot(gs[0, 0])
ax2 = fig.add_subplot(gs[0, 1:])

ax1.set_yscale('log')
ax2.set_xscale('log')
sns.boxplot(x='Region', y='R_Length', data=dfAll, ax=ax1, palette="vlag")

sns.boxplot(x='R_Length', y='SwarmID', data=dfAll, ax=ax2, palette="vlag")
# sns.stripplot(x='R_Length', y='SwarmID',data=dfAll, ax=ax[1],
#               size=.8, color=".05", linewidth=0)
ax2.xaxis.grid(True)
ax2.set(ylabel="")
sns.despine(offset=5)
labelSubplots([ax1, ax2], fontsize=26)
fig.savefig('/home/akh/myprojects/Linking-and-Clustering-Dikes/dikedata/All_Lengths.png')
fig.savefig('/home/akh/myprojects/Linking-and-Clustering-Dikes/dikedata/All_Lengths.svg', dpi=600)


sns.set_theme(style='ticks')
fig = plt.figure(figsize=(12,6), constrained_layout=True)

gs = GridSpec(1, 4, figure=fig)
ax1 = fig.add_subplot(gs[0, 0])
ax2 = fig.add_subplot(gs[0, 1:])

ax1.set_yscale('log')
ax2.set_xscale('log')
sns.boxplot(x='Region', y='R_Length', data=dfAllTrusted, ax=ax1, palette="vlag")

sns.boxplot(x='R_Length', y='SwarmID', data=dfAllTrusted, ax=ax2, palette="vlag")
# sns.stripplot(x='R_Length', y='SwarmID',data=dfAll, ax=ax[1],
#               size=.8, color=".05", linewidth=0)
ax2.xaxis.grid(True)
ax2.set(ylabel="")
sns.despine(offset=5)
labelSubplots([ax1, ax2], fontsize=26)
fig.savefig('/home/akh/myprojects/Linking-and-Clustering-Dikes/dikedata/All_LengthsTrusted.png')
fig.savefig('/home/akh/myprojects/Linking-and-Clustering-Dikes/dikedata/All_LengthsTrusted.svg', dpi=600)