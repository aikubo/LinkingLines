#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 21 14:07:06 2021

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

sns.set_context("talk")
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

df=pd.DataFrame({'Xstart':Xstart, 'Xend': Xend, 'Ystart': Ystart, 'Yend':Yend, 'Slope':m2})

fig, ax =plt.subplots(2,2)


plotlines(df, 'k', ax[0][0])
theta1,rho1, xc, yc=AKH_HT(df.astype(float), xc=0, yc=0)
ax[1][0].scatter(theta1, rho1, edgecolor='black')
ax[0][0].plot(0,0, "*", mec='black', markersize=20)
ax[1][1].set_xlabel('Theta (deg)')
ax[1][0].set_xlabel('Theta (deg)')

ax[0][1].set_label('Perfect Radial Swarm')
ax[0][0].set_label('With Noise')

ax[1][0].set_ylabel('Rho (m)')

dfFrag1=fragmentDikes(df, distortion=1)
d=20000
plotlines(dfFrag1, 'r', ax[0][1], linewidth=2)
plotlines(df, 'k', ax[0][1], alpha=0.1)
theta2, rho2, xc2, yc2=AKH_HT(dfFrag1.astype(float), xc=0, yc=0)
ax[1][1].scatter(theta2, rho2, edgecolor='black')
ax[0][1].plot(0,0, "*", mec='black', markersize=20)

""" make gif of center change """
fig,ax=plt.subplots(1,2)
fig.set_size_inches(12,6)
rold=rho1
told=theta1
n=0

plotlines(df, 'grey', ax[0], center=True)
plotlines(dfFrag1, 'r', ax[0], linewidth=2)

ax[0].plot(0,0, "*", mec='black', markersize=20)
ax[1].scatter(theta1,rho1, edgecolor='black')
ax[1].set_ylim([-30000, 30000,])
ax[1].set_xlabel('Theta (deg)')
ax[1].set_ylabel('Rho (m)')

#title="Center at ["+str(0) +"m ,"+str(0)+" m]"
#ax[0].set_title(title)
n=0
plt.tight_layout()
name="radial"+str(n)+".png"
fig.savefig(name, dpi=600)
plt.tight_layout()
xcs=np.array([0,d,d,d,0,-d,-d,-d])
ycs=np.array([d,d,0,-d,-d,-d,0,d])


for ic in range(len(xcs)):

    i=xcs[ic]
    j=ycs[ic]
    fig,ax=plt.subplots(1,2)
    fig.set_size_inches(12,6)
    ax[1].scatter(told,rold, c='grey', s=20)
    
    plotlines(df, 'grey', ax[0], center=True)
    plotlines(dfFrag1, 'r', ax[0], linewidth=2)
    theta2, rho2, xc2, yc2=AKH_HT(dfFrag1.astype(float), xc=i, yc=j)
    ax[0].plot(i,j, "*", mec='black', markersize=20)
    ax[1].scatter(theta2,rho2, edgecolor='black')
    ax[1].set_ylim([-30000, 30000,])
    ax[1].set_xlabel('Theta (deg)')
    ax[1].set_ylabel('Rho (m)')
    
    #title="Center at ["+str(i) +"m ,"+str(j)+" m]|"
    #ax[0].set_title(title)
    rold=np.append(rold, rho2)
    told=np.append(told, theta2)
    n=n+1
    plt.tight_layout()
    name="radial"+str(n)+".png"
    fig.savefig(name, dpi=600)
""" """
# dfFrag2=fragmentDikes(df, ndikesMax=40, distortion=0)
# dfFrag3=fragmentDikes(df, ndikesMax=20, distortion=5)


# 
# theta2, rho2, xc2, yc2=AKH_HT(dfFrag1.astype(float), xc=20000, yc=20000)
# theta3, rho3, xc3, yc3=AKH_HT(dfFrag1.astype(float), xc=-20000, yc=20000)
# theta4, rho4, xc4, yc4=AKH_HT(dfFrag1.astype(float), xc=20000, yc=-20000)
# theta5, rho5, xc5, yc5=AKH_HT(dfFrag1.astype(float), xc=-20000, yc=-20000)

# plotlines(dfFrag1, 'r', ax[0][1], linewidth=3, center=True, xc=xc2, yc=yc2)

# plotlines(df, 'grey', ax[0][1], alpha=0.1)

# plotlines(dfFrag2, 'r', ax[0][2], linewidth=3)
# plotlines(df, 'grey', ax[0][2], alpha=0.1, center=True, xc=xc3, yc=yc3)
# plotlines(dfFrag3, 'r', ax[0][3], linewidth=3)
# plotlines(df, 'grey', ax[0][3], alpha=0.1, center=True, xc=xc4, yc=yc4)
                    
# 
# ax[1][1].scatter(theta2, rho2)
# ax[1][2].scatter(theta3, rho3)
# ax[1][3].scatter(theta4, rho4)
# ax[1][0].set_ylabel('Rho (m) ')

# for i in range(4):
#     ax[1][i].set_xlabel('Theta (deg)')

# ax[0][0].set_title('Centered')
# ax[0][1].set_title('[20000,20000]')
# ax[0][2].set_title('[-20000,20000]')
# ax[0][3].set_title('[20000,-20000]')
# ax[0][4].set_title('[-20000,-20000]')

# fig, ax =plt.subplots(2,4)
# # 1 complete info
# # 2 Short dikes only 
# # 3 Little Dikes only 
# # 4 Short few dikes 
# plotlines(df, 'k', ax[0][0])

# dfFrag1=fragmentDikes(df, maxL=5000)
# dfFrag2=fragmentDikes(df, maxL=10000, ndikesMax=20)
# dfFrag3=fragmentDikes(df, maxL=1000, ndikesMax=20)


# theta1,rho1, xc, yc=AKH_HT(df.astype(float))
# theta2, rho2, xc, yc=AKH_HT(dfFrag1.astype(float))
# theta3, rho3, xc, yc=AKH_HT(dfFrag2.astype(float))
# theta4, rho4, xc, yc=AKH_HT(dfFrag3.astype(float))


# plotlines(dfFrag1, 'r', ax[0][1], linewidth=3)
# plotlines(df, 'grey', ax[0][1], alpha=0.1)

# plotlines(dfFrag2, 'r', ax[0][2], linewidth=3)
# plotlines(df, 'grey', ax[0][2], alpha=0.1)
# plotlines(dfFrag3, 'r', ax[0][3], linewidth=3)
# plotlines(df, 'grey', ax[0][3], alpha=0.1)
                    
# ax[1][0].scatter(theta1, rho1)
# ax[1][1].scatter(theta2, rho2)
# ax[1][2].scatter(theta3, rho3)
# ax[1][3].scatter(theta4, rho4)
# ax[1][0].set_xlabel('Rho (m) ')
# for i in range(4):
#     ax[1][i].set_xlabel('Theta (deg)')

# ax[0][0].set_title('Complete Info')
# ax[0][1].set_title('Short Dikes only')
# ax[0][2].set_title('Few Dikes')
# ax[0][3].set_title('Few Short Dikes')


# dikelength=500000

# center=np.array([0,0])
# ndikes=100
# angles=np.random.normal(45, 2, ndikes)
# m=np.tan(angles)
# Xstart=np.random.randint(0,high=dikelength, size=ndikes)
# Ystart=np.random.randint(0,high=dikelength, size=ndikes)
# Xend=dikelength/np.sqrt(1+m**2)+Xstart
# Yend=m*Xend

# df=pd.DataFrame({'Xstart':Xstart, 'Xend': Xend, 'Ystart': Ystart, 'Yend':Yend, 'Slope':m})

# fig, ax =plt.subplots(2,2)
# plotlines(df, 'k', ax[0][0])

# dfFrag=fragmentDikes(df)

# theta1,rho1, xc, yc=AKH_HT(df.astype(float))
# theta2, rho2, xc, yc=AKH_HT(dfFrag.astype(float))

# plotlines(dfFrag, 'r', ax[0][1])
# ax[1][0].scatter(rho1, theta1)
# ax[1][1].scatter(rho2, theta2)
