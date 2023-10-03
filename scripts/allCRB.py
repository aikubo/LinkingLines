#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 27 14:50:14 2021

@author: akh
"""
from scipy.cluster import hierarchy
from scipy.spatial.distance import squareform
import pandas as pd
from htMOD import AKH_HT as HT
from sklearn.preprocessing import scale
import numpy as np 
from clusterMod import *
import matplotlib.pyplot as plt
from plotmod import plotlines, labelcolors, plotbyAngle, BA_HT, plotResults
from examineMod import examineClusters, plotlabel,extendLines
from PrePostProcess import *
from matplotlib import cm
from findRadialCenters import sweepCenters, detectLocalPeaks
from matplotlib.lines import Line2D  

from matplotlib import cm
from skimage.morphology import reconstruction

from scipy.ndimage import gaussian_filter

from skimage import img_as_float
from skimage.filters import threshold_otsu, threshold_local
import numpy as np
import matplotlib.pyplot as plt
from scipy import ndimage as ndi
from plotmod import plotlines, HT3D
from skimage.segmentation import watershed
from skimage.feature import peak_local_max

from skimage import feature
from htMOD import HT_center
from findRadialCenters import sweepCenters, detectLocalPeaks
from scipy.spatial.distance import pdist, squareform
import scipy.cluster.hierarchy as sch

from examineMod import persitance

# lines=pd.read_csv("/home/akh/myprojects/Linking-and-Clustering-Dikes/AllCRBLinked.csv")
# xc,yc=HT_center(lines)
# gridM=1000

# lines=MidtoPerpDistance(lines, xc, yc)

# #X=(np.vstack(( np.deg2rad(theta), rho, lines['PerpOffsetDist'])).T)
# midist=np.log(lines['PerpOffsetDist'].values)

# HT3D(lines['AvgTheta'],lines['AvgRho'],midist, midist)

# err, dikes, thetaRange, thetaSTD, xs, ys=sweepCenters(lines, gridM, 5000, xc,yc, plot=True)
# lines,labels,counts, centerData=detectLocalPeaks(err,dikes, lines,gridM, plot=True)
23

# xc=475201.3737670694
# yc=4976341.597868606

dikeset=pd.read_csv("/home/akh/myprojects/Linking-and-Clustering-Dikes/dikedata/allCRB_dikes.csv")
dikeset=completePreProcess(dikeset)

swarms=['CJDS', 'Ice_Harbor', 'Monument_Dikes',
        'Steens_dikes']

allCRB={}
#fig,ax=plt.subplots(1,2)
xc,yc=HT_center(dikeset)
dikeset= MidtoPerpDistance(dikeset, xc, yc)
midist=dikeset['PerpOffsetDist'].values

theta,rho, xc, yc=HT(dikeset)
dikeset['theta']=theta
dikeset['rho']=rho
col=['#30123BFF', '#1AE4B6FF', '#FABA39FF', '#7A0403FF']
allLinked=pd.DataFrame()

fig,ax=plt.subplots(1,2)
drho=dikeset['seg_length'].mean()

dtheta=1
drho=700
# fig, ax, Z1=persitance(dikeset)
#distancesPlot(dikeset)

# x,y=np.meshgrid( np.linspace(0,90,100), np.linspace(0, max(rho), 100))
# d=np.sqrt( (x/dtheta)**2 + (y/drho)**2)
# fig,ax=plt.subplots()
# c=ax.contourf(x,y,d, levels=[1, 10, 30, 100, 1000],  colors=('r','orange', 'y', 'g', 'b'), extend='max')
# fig.colorbar(c)


data,clusters, M=HT_AGG_custom(dikeset, dtheta, drho, linkage='complete')
lines,IC=examineClusters(data)



plotResults(lines)
dataOld=pd.read_csv('/home/akh/myprojects/Linking-and-Clustering-Dikes/AllCRBLinked.csv')

for i in [dataOld.loc[ dataOld['Linked']==1], lines.loc[ lines['Linked']==1]]:
    
    print("length mean:", i['R_Length'].mean())
    print("length max:", i['R_Length'].max())
    print("width mean:", i['R_Width'].mean())
    print("width max:", i['R_Width'].max())
    print("Cluster Size Mean:", i['Size'].mean())
    print("Cluster Size Max:", i['Size'].max())

hmm=lines[lines['ThetaRange']>10]

hmmSegments=pd.DataFrame()
for i in hmm['Label']:
    hmmSegments=hmmSegments.append(dikeset[dikeset['Labels']==i])

# for dtheta in [0.5, 0.75, 1.0, 1.25, 1.5] : #, drho in zip( [0.5, 1, 1.5, 2], [100,500,1000, 2000]):
#     data,clusters=HT_AGG_custom(dikeset, dtheta, drho)
#     lines,IC=examineClusters(data)
    
#     #n clusters
#     ax[0,0].scatter(dtheta, len(np.unique(data['Labels'])))
#     ax[1,0].scatter(drho, len(np.unique(data['Labels'])))
    
#     # avg length
#     ax[0,1].scatter(dtheta, lines['R_Length'].mean())
#     ax[1,1].scatter(drho, lines['R_Length'].mean())
    
#     #avg width
#     ax[0,2].scatter(dtheta, lines['R_Width'].mean())
#     ax[1,2].scatter(drho, lines['R_Width'].mean())
    
#     #avg size
#     ax[0,3].scatter(dtheta, lines['Size'].mean())
#     ax[1,3].scatter(drho, lines['Size'].mean())
    
    
# ax[0,0].set_title('Number of Clusters')
# ax[0,1].set_title('Average length')
# ax[0,2].set_title('Average width')
# ax[0,3].set_title('Size of Clusters')
# for i in range(4):
#     ax[0,i].set_xlabel('Theta Cutoff')
#     ax[1,i].set_xlabel('Rho Cutoff')

#b=0
# for i in swarms: 
    
#     df=dikeset[dikeset['layer']==i]

#     theta, rho, xc, yc= HT(df, xc=xc, yc=yc)
#     df['rho']=rho
#     df['theta']=theta
    
#     trange=2
#     rrange=5000 
    
#     stdT=np.std(theta)
#     stdRho=np.std(rho)
#     d2=trange/stdT
#     dfromCJDS=0.07184050241647064
#     d=d2
#     clustering=HT_AGG(df,d)
#     df['Labels']=clustering.labels_
#     lines,IC=examineClusters(df)
#     lines['Layer']=[i]*len(lines)
#     allCRB[i]=df
#     allCRB[i+"Linked"]=lines
#     allCRB[i+"Center"]=[xc,yc]
#     allLinked=allLinked.append(lines)
#     plotlines(lines, col[b], ax[0])
#     ax[1].scatter(theta, rho, color=col[b], alpha=0.6, edgecolor='k')
#     #longlines=extendLines(lines, save=True, name=str(i)+'_extendedlines.csv')
#     b=b+1

# ax[1].set_ylabel("Rho (m)")
# ax[1].set_xlabel("Theta ($^\circ$)")
# custom_lines = [Line2D([0], [0], color=col[0], lw=4),
#         Line2D([0], [0], color=col[1], lw=4),
#         Line2D([0], [0], color=col[2], lw=4),
#         Line2D([0], [0], color=col[3], lw=4)]

# ax[0].plot(xc,yc, 'r*')
# ax[1].legend(custom_lines, swarms)


