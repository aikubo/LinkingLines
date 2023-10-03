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
from PublicationFigures import *
import seaborn as sns
from matplotlib.gridspec import GridSpec

d=['dikedata/deccandata/AllDeccan_PreProcessed.csv',
    'dikedata/crb/allCRB_dikes_PreProcessed.csv',
    'dikedata/spanish peaks/SpanishPeaks_3857_preprocessed.csv'
    ]

l=['dikedata/deccandata/AllDeccanLinked_euc_18_10_2022.csv',
    'dikedata/crb/AllCRBLinked_euc_18_10_2022.csv',
    'dikedata/spanish peaks/SpanishPeaks_Complete_euc_3_2013.csv']

c=['r', 'g', 'b']
moveX=[0, -7.0e5, -1.0e6]
moveY=[0, 3.5e5, 5.24e5]

fig,ax=plt.subplots()
fig.set_size_inches( 12,6)

fig2,ax2=plt.subplots(3, sharex=True)
ax2[0].set_ylabel('Dike Length (km)')
#ax2[0].set_xlabel("Number of Dikes")
ax2[0].set_yscale('log')
ax2[0].set_xscale('log')

ax2[1].set_ylabel('Swarm Area (km $^{2}$)')
#ax2[1].set_xlabel("Number of Dikes")
ax2[1].set_yscale('log')
ax2[1].set_xscale('log')

ax2[2].set_ylabel('Erupted Volume (km $^{3}$)')
ax2[2].set_xlabel("Number of Dikes")
ax2[2].set_yscale('log')
ax2[2].set_xscale('log')




averagewidth=[2,8,8]
name=['Deccan', 'CJDS', 'Spanish Peaks']
dfAll=pd.DataFrame() 
EruptedVolume=[ 1.3e6, 210000, 14] #50,210000, 1.3e6] #km
# CRBG estimate: Kasbohm 2018
# Deccan estimate: Jay and Widdowson 2008
# Spanish peaks erupted volume: 14 km^3
    
figAngle,axAngle=plt.subplots(3, sharex=True)
figAngle.set_size_inches(4,9)


TrustOn=False  # only use the clustered dikes that are trusted
for dikeset,lines, color, xnew, ynew, w, label, ev, histax in zip(d,l,c, moveX, moveY, averagewidth, name, EruptedVolume, axAngle): 
    
    df=pd.read_csv(dikeset)
    #df=DikesetReProcess(df)
    

    xc,yc=HT_center(df)
    df2=moveHTcenter(df)
    df2=moveHTcenter(df2, xc=xnew, yc=ynew)
    df3=pd.read_csv(lines)
    df3['Region']=label
        
    df3['R_Length'] = pd.to_numeric(df3['R_Length'], errors='coerce')
    df3['R_Width'] = pd.to_numeric(df3['R_Width'], errors='coerce')
    df3['R_Width']=df3['R_Width']+w
    df3['Size'] = pd.to_numeric(df3['Size'], errors='coerce')
    #m=df3['TrustFilter']==1
    df3=df3.assign(NormPerpOffsetDist=df3['PerpOffsetDist'].values/np.mean(df3['PerpOffsetDist'].values))
    if TrustOn:
        df3=df3[m]
    
    xs=[ min(df2.min().iloc[ [0,2]]), max(df2.max().iloc[ [0,2]])]
    ys=[ min(df2.min().iloc[ [1,3]]), max(df2.max().iloc[ [1,3]])]
    width=np.ptp(xs)
    height=np.ptp(ys)
    Xedges=np.array([xs[0], xs[0], xs[1], xs[1], xs[0]])
    Yedges=np.array([ys[1], ys[0], ys[0], ys[1], ys[1]])
    ax.plot(Xedges, Yedges, color, linewidth=4)
    rect=Rectangle( (xs[0], ys[0]), width, height,facecolor=color, alpha=0.3, label=label)
    ax.add_patch(rect)
    ax.legend(loc='upper right')
    FixCartesianLabels(ax)

    AngleHistograms(df, df3, ax=histax)
    histax.set_title(label)
    


    n=len(df)
    #plotlines(df2, color, ax, linewidth=0.5, alpha=0.4)
    
    #comparison plots
    
    # segment Length
    ax2[0].plot(n,  df['seg_length'].mean(), "s", color=color, label=label)
    ax2[0].plot(n,  df3['R_Length'].mean(), "P", color=color, label=label)
    
    # width*length of swarm area
    ax2[1].plot(n,  width*height/1000**2, "*", color=color, label=label)
    # linked dike length * width (not packet width) (only 4+segments)
    
    l=df3['R_Length'].mean()
    ax2[1].plot(n,  l*w/1000**2, "p", color=color, label=label)
    
    #erupted volume 
    ax2[2].plot(n,  ev, "*", color=color, label=label)
    ax2[2].legend()
    plt.tight_layout()
    
    
    # ax2.plot(n, width*height/1000**2, "*", color=color, label=name)
    # ax2.plot(n, df['seg_length'].sum()*w/1000**2, "p", color=color )
    # a=df3['R_Length'].values*df3['R_Width'].values
    # a[np.isnan(a)]=0
    # a=np.sum(a)/1000**2
    # ax2.plot(len(df3), a, 'D', color=color)
    # print(width*height/1000**2,  df['seg_length'].sum()*w/1000**2 ,a)
    
    # df3=df3.loc[df3['Size']>4]
    
    # a=df3['R_Length'].values*df3['R_Width'].values
    # a[np.isnan(a)]=0
    # a=np.sum(a)/1000**2
    # ax2.plot(len(df3), a, 'v', color=color)
    
    dfAll=pd.concat([dfAll, df3]).reset_index(drop=True)

axAngle[2].set_xlabel('Theta ($^\circ$)')
for i in axAngle:
    i.set_ylabel('Normalized Counts')
c=axAngle[2]
c.legend()
labelSubplots(axAngle, fontsize=18)
plt.tight_layout()
figAngle.savefig('/home/akh/myprojects/Linking-and-Clustering-Dikes/dikedata/AngleHistograms.png', dpi=600)
figAngle.savefig('/home/akh/myprojects/Linking-and-Clustering-Dikes/dikedata/AngleHistograms.pdf', dpi=600)


dikeset=pd.read_csv('/home/akh/myprojects/Linking-and-Clustering-Dikes/dikedata/AllDatasets_18_10_2022.csv')
dfAll=pd.read_csv('/home/akh/myprojects/Linking-and-Clustering-Dikes/dikedata/AllDatasetsLinked_18_10_2022.csv')
dfAll=dfAll.assign(R_Width=dfAll['R_Width']+8)
dfAll['Dike Cluster Length (km)']=dfAll['R_Length']/1000
dfAll['Dike Cluster Width (m)']=dfAll['R_Width']+2
dikeset['Segment Length (km)']=dikeset['seg_length']/1000
dfAll=dfAll.assign(Status=[ 'Filtered Cluster' if i==1 else 'Cluster' for i in dfAll['TrustFilter']])
dfAllTrusted=dfAll[ dfAll['TrustFilter']==1]
dfAll=dfAll[dfAll['Linked']==1]

for i in dfAll['SwarmID'].unique():
    
    
    print(i)
    print('PerpOffsetDist')
    print( 'Mean', dikeset[dikeset['SwarmID']==i]['PerpOffsetDist'].mean())
    print( 'Median', dikeset[dikeset['SwarmID']==i]['PerpOffsetDist'].median())
    print( 'Max', dikeset[dikeset['SwarmID']==i]['PerpOffsetDist'].max())
    print( 'Std', dikeset[dikeset['SwarmID']==i]['PerpOffsetDist'].std())
    
    
    print('segments')
    print(len(dikeset[dikeset['SwarmID']==i]))
    print( 'Mean', dikeset[dikeset['SwarmID']==i]['seg_length'].mean())
    print( 'Median',dikeset[dikeset['SwarmID']==i]['seg_length'].median())
    print( 'Max', dikeset[dikeset['SwarmID']==i]['seg_length'].max())


    # print('Linked')
    # print(len(dfAll[dfAll['SwarmID']==i]))
    # print( 'Mean', dfAll[dfAll['SwarmID']==i]['R_Length'].mean())
    # print( 'Median', dfAll[dfAll['SwarmID']==i]['R_Length'].median())
    # print( 'Max', dfAll[dfAll['SwarmID']==i]['R_Length'].max())

    
    # print('Linked Trusted')
    # print(len(dfAll[dfAll['SwarmID']==i]))
    # print( 'Mean', dfAllTrusted[dfAllTrusted['SwarmID']==i]['R_Length'].mean())
    # print( 'Median', dfAllTrusted[dfAllTrusted['SwarmID']==i]['R_Length'].median())
    # print( 'Max', dfAllTrusted[dfAllTrusted['SwarmID']==i]['R_Length'].max())


    
          

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

sns.boxplot(x='Region', y='Segment Length (km)', data=dikeset, ax=ax1, palette="vlag")

sns.boxplot(x='Segment Length (km)', y='SwarmID', data=dikeset, ax=ax2, palette="vlag")
# sns.stripplot(x='R_Length', y='SwarmID',data=dfAll, ax=ax[1],
#               size=.8, color=".05", linewidth=0)
ax2.xaxis.grid(True)
ax2.set(ylabel="")
sns.despine(offset=5)
labelSubplots([ax1, ax2], fontsize=26)
fig.savefig('/home/akh/myprojects/Linking-and-Clustering-Dikes/dikedata/All_seg_Lengths.png')
fig.savefig('/home/akh/myprojects/Linking-and-Clustering-Dikes/dikedata/All_seg_Lengths.svg', dpi=600)


fig,ax=plt.subplots(1,3)
fig.set_size_inches(12,5)
#temp=pd.read_csv('/home/akh/myprojects/Linking-and-Clustering-Dikes/dikedata/deccandata/AllDeccan_PreProcessed.csv')
#DotsHT(fig, ax[0],temp, ColorBy='SwarmID', palette="hls", rhoScale=True)

j=1
palettes=['flare', 'RdBu']
for i, p in enumerate([ 'Deccan','CRBG', 'SpanishPeaks']):
    temp=dikeset[dikeset['Region']==p]
    if p=='Deccan':
        temp=pd.read_csv('/home/akh/myprojects/Linking-and-Clustering-Dikes/dikedata/deccandata/AllDeccan_PreProcessed.csv')
    theta,rho, x, y=HT(temp)
    temp=temp.assign(theta=theta, rho=rho)
    ax[i].set_xlabel('Theta ($^\circ$)')
    ax[i].set_ylabel('Rho (km)')
    sns.scatterplot(temp, x='theta', y=rho/1000, hue='theta', style='SwarmID', palette='viridis', ax=ax[i], alpha=0.4, edgecolor="black")
    #DotsHT(fig, ax[j],temp, ColorBy='SwarmID', palette="hls", rhoScale=True)
    j=j+1
plt.tight_layout()


fig.savefig('/home/akh/myprojects/Linking-and-Clustering-Dikes/dikedata/All_HT3_coloredbyTheta.png')


fig,ax=plt.subplots(1,3)
fig.set_size_inches(12,5)
#temp=pd.read_csv('/home/akh/myprojects/Linking-and-Clustering-Dikes/dikedata/deccandata/AllDeccan_PreProcessed.csv')
#DotsHT(fig, ax[0],temp, ColorBy='SwarmID', palette="hls", rhoScale=True)

j=1
palettes=['flare', 'RdBu']
for i, p in enumerate([ 'Deccan','CRBG', 'SpanishPeaks']):
    temp=dikeset[dikeset['Region']==p]
    if p=='Deccan':
        temp=pd.read_csv('/home/akh/myprojects/Linking-and-Clustering-Dikes/dikedata/deccandata/AllDeccan_PreProcessed.csv')
    theta,rho, x, y=HT(temp)
    temp=temp.assign(theta=theta, rho=rho)
    ax[i].set_xlabel('Theta ($^\circ$)')
    ax[i].set_ylabel('Rho (km)')
    sns.scatterplot(temp, x='theta', y=rho/1000, hue='SwarmID', style='SwarmID', legend='brief', palette='viridis', ax=ax[i], alpha=0.4, edgecolor="black")
    #DotsHT(fig, ax[j],temp, ColorBy='SwarmID', palette="hls", rhoScale=True)
    j=j+1
plt.tight_layout()


fig.savefig('/home/akh/myprojects/Linking-and-Clustering-Dikes/dikedata/All_HT3_swarmid.png')
fig.savefig('/home/akh/myprojects/Linking-and-Clustering-Dikes/dikedata/All_HT3_swarmid.svg', dpi=600)


sns.set_theme(style='ticks')
fig = plt.figure(figsize=(12,6), constrained_layout=True)

dfAll=dfAll.assign(TwistAngle= dfAll['EnEchelonAngleDiff'].values)

dikesetTemp=pd.DataFrame({'Length (km)':dikeset['Segment Length (km)'].values, 'Status': ['Segments']*len(dikeset),'Region': dikeset['Region'], 'Size':np.ones(len(dikeset)),'SwarmID':dikeset['SwarmID'], 'Dike Cluster Width (m)':[10]*len(dikeset)})
dfAllTemp=pd.DataFrame({'Length (km)':dfAll['Dike Cluster Length (km)'].values, 'Status': ['Cluster']*len(dfAll),'Region': dfAll['Region'], 'SwarmID':dfAll['SwarmID'], 'Size': dfAll['Size'], 'Dike Cluster Width (m)':dfAll['R_Width'].values+1})
dfAllTrustedTemp=pd.DataFrame({'Length (km)':dfAllTrusted['Dike Cluster Length (km)'].values, 'Status': ['Filtered Cluster']*len(dfAllTrusted),'Region': dfAllTrusted['Region'],  'Size': dfAllTrusted['Size'],'SwarmID':dfAllTrusted['SwarmID'], 'Dike Cluster Width (m)':dfAllTrusted['R_Width'].values+1})
All_AllDf=pd.concat((dfAllTemp, dfAllTrustedTemp), ignore_index=True)
All_AllDf=pd.concat( (dikesetTemp, All_AllDf), ignore_index=True)
All_lines=pd.concat((dfAllTemp, dfAllTrustedTemp), ignore_index=True)

SegsTrusted=pd.concat( (dikesetTemp, dfAllTrustedTemp))
dfLIP=dfAllTrusted.iloc[np.in1d(dfAllTrusted['SwarmID'].values, ['Deccan:Central','Deccan:Narmada-Tapi','Deccan:Saurashtra','CRBG:CJDS'])]


sns.set_theme(style='ticks')

fig,ax=PubBoxnWhisker(All_AllDf, hue='Status', l='Length (km)')

fig.savefig('/home/akh/myprojects/Linking-and-Clustering-Dikes/dikedata/PubWidthLengthsBNW.pdf', dpi=600)

fig,ax=PubBoxnWhisker(SegsTrusted, hue='Status')

fig.savefig('/home/akh/myprojects/Linking-and-Clustering-Dikes/dikedata/PubSegsTrusted_WidthLengthsBNW.pdf', dpi=600)




fig1,ax=PubWidthsVsLengths(dfAllTrusted , hue='SwarmID', style='Region')
fig1.savefig('/home/akh/myprojects/Linking-and-Clustering-Dikes/dikedata/WidthVsLength_trusted2.pdf', dpi=600)

fig,ax=PubWidthsVsLengths(dfLIP, hue='SwarmID', style='Region')
fig.savefig('/home/akh/myprojects/Linking-and-Clustering-Dikes/dikedata/WidthVsLength_LIPSonly2.pdf', dpi=600)

fig,ax=PubWidthsVsLengths(dfAllTrusted, hue='SwarmID', style='Region')
fig.savefig('/home/akh/myprojects/Linking-and-Clustering-Dikes/dikedata/WidthVsLength_byswarm2.pdf', dpi=600)


m=All_lines['Status']=='Filtered Cluster'

g=sns.jointplot(data=All_lines, x='Size', y='Length (km)', alpha=0.6, hue='Status', ylim=(-100,1250), xlim = (1,30))

slope, intercept, r_value, p_value, std_err = stats.linregress(All_lines['Size'].values, All_lines['Length (km)'].values)
x=np.linspace(2, All_lines['Size'].max())
y=slope*x+intercept 

g.ax_joint.plot(x,y, 'b-.', alpha=0.8)
g.ax_joint.annotate(f'$m = {slope:.1f}, R^2 = {r_value:.2f}$',
                xy=(0.6, 0.6), xycoords='axes fraction',
                ha='left', va='center',
                bbox={'boxstyle': 'round', 'fc': 'powderblue', 'ec': 'navy'})
   

slope, intercept, r_value, p_value, std_err = stats.linregress(All_lines['Size'].loc[m].values, All_lines['Length (km)'].loc[m].values)
x=np.linspace(3, All_lines['Size'].max())
y=slope*x+intercept 

g.ax_joint.plot(x,y,  '-.',color='orange', alpha=0.9)
g.ax_joint.annotate(f'$m = {slope:.1f}, R^2 = {r_value:.2f}$',
                xy=(0.6, 0.7), xycoords='axes fraction',
                ha='left', va='center',
                bbox={'boxstyle': 'round', 'fc': 'peachpuff', 'ec': 'coral'})

#g.ax_joint.set_yscale('log')
g.savefig("/home/akh/myprojects/Linking-and-Clustering-Dikes/dikedata/SizebyLength_All.pdf", transparent=True)



g=sns.scatterplot(x='Overlap', y='TwistAngle',  hue='Region', data=dfAll, palette=p)
#g.save("/home/akh/myprojects/Linking-and-Clustering-Dikes/dikedata/Twist_Overlap_All.pdf", transparent=True,dpi=600)


# fig,ax=plt.subplots(); ax.scatter(dfAll['NormPerpOffsetDist'], dfAll['R_Width'], alpha=0.4)
# #ax.set_yscale('log')
# #ax.set_xscale('log')
# ax.set_ylabel('Dike Cluster Width (m)')
# ax.set_xlabel('Normalized Perpendicular Offset')
# x=dfAll['NormPerpOffsetDist'].values
# y=dfAll['R_Width']
# slope, intercept, r_value, p_value, std_err = stats.linregress(x,y)
# x2=np.linspace(min(x), max(x))
# y2=slope*x2+intercept 
# ax.plot(x2,y2, 'r-')

# ax.annotate(f'$R^2 = {r_value:.2f}$',
#                 xy=(0.8, 0.85), xycoords='axes fraction',
#                 ha='left', va='center',
#                 bbox={'boxstyle': 'round', 'fc': 'white', 'ec': 'red'})
# fig.savefig("/home/akh/myprojects/Linking-and-Clustering-Dikes/dikedata/PerpOffset_Width.pdf", transparent=True,dpi=600)

plotScatterHist(dfAll, 'NormPerpOffsetDist', 'Dike Cluster Width (m)')

fig2,ax=PubTwistOverlap(dfAllTrusted)
fig2.savefig("/home/akh/myprojects/Linking-and-Clustering-Dikes/dikedata/TwistOverlap.pdf", dpi=600)

fig,ax=PubTwistOverlap(dfAllTrusted, partial=False)
fig.savefig("/home/akh/myprojects/Linking-and-Clustering-Dikes/dikedata/TwistOverlap_all.pdf", dpi=600)


calculated_overlap=dfAllTrusted['SegmentLSum'].values/dfAllTrusted['Size']-dfAllTrusted['StdRho'].values/2

fig, ax=plt.subplots()
sns.scatterplot(data=dfAllTrusted, x=dfAllTrusted['Overlap'].values*dfAllTrusted['SegmentLSum'].values, y=calculated_overlap, ax=ax)
ax.plot(calculated_overlap,calculated_overlap, 'r-.', label='1:1')
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_ylim((min(dfAll['Overlap'].values*dfAll['SegmentLSum'].values), max(dfAll['Overlap'].values*dfAll['SegmentLSum'].values)))
ax.set_xlim((min(dfAll['Overlap'].values*dfAll['SegmentLSum'].values), max(dfAll['Overlap'].values*dfAll['SegmentLSum'].values)))
ax.legend()
ax.set_xlabel('Observed Overlap')
ax.set_ylabel('Calculated Overlap (Pollard 1987)')