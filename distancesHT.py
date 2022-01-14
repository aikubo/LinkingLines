#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 20 11:33:24 2021

@author: akh

runfile('/home/akh/myprojects/Linking-and-Clustering-Dikes/distancesHT.py', wdir='/home/akh/myprojects/Linking-and-Clustering-Dikes')
"""
import numpy as np 
from plotmod import plotlines, DotsHT
import pandas as pd
from htMOD import AKH_HT as HT 
import matplotlib.pyplot as plt
from sklearn.cluster import DBSCAN
import matplotlib.lines as mlines
from sklearn.preprocessing import scale
from scipy.spatial.distance import pdist, squareform
from htMOD import MidtoPerpDistance
from examineMod import *
from matplotlib import cm
import scipy
import pylab
import scipy.cluster.hierarchy as sch
from scipy.spatial.distance import squareform
from sklearn.cluster import AgglomerativeClustering as AGG
from scipy.cluster.hierarchy import dendrogram

from PrePostProcess import segLength, giveHashID
def CylDistanceAngleModified(u, v):
    '''
    Calculates cylindrical distance between two points in 3D with angle modification

    Input:
        u: [theta, rho, z]
        v: [theta, rho, z]
    Output:
        dist: distance between u and v
    '''
    
    dist1 = u[0]-v[0]
    dist2 = dist1%180
    deltaT = min(dist1, dist2)
    return np.sqrt(u[1]**2 + v[1]**2 + 2*u[1]*v[1]*np.cos(np.deg2rad(deltaT)) + (u[2]-v[2])**2)

def PolarDistanceAngleModified(u, v):
    '''
    Calculates polar distance between two points in 2D with angle modification

    Input:  
        u: [theta, rho]
        v: [theta, rho]
    Output:
        dist: distance between u and v

    '''
    
    dist1 = abs(u[0]-v[0])
    dist2 = dist1%180
    deltaT = min(dist1, dist2)
    return np.sqrt(u[1]**2 + v[1]**2 + 2*u[1]*v[1]*np.cos(np.deg2rad(deltaT)))

def DistortionDistance(u,v):
    '''
    Calculates distance based on distortion calculation in 2D

    Input:
        u: [theta, rho]
        v: [theta, rho]
    Output:
        dist: distance between u and v

        For u==v, dist=0
        for u[0]==v[0], dist=u[1]-v[1]
        for u[1]==v[1] and dtheta-> 0, dist ->0
        )
    '''
    tol=0.0000000001
    alpha=2
    dtheta=abs(u[0]-v[0])%180
    d1=abs(v[1]-u[1]/np.cos(np.deg2rad(dtheta)))
    d2=abs(v[1]/np.cos(np.deg2rad(dtheta))-u[1])
    
    return abs((1+alpha*dtheta)*max(d1,d2))

def DistortionDistance2(u,v):
    '''
    Calculates distance based on distortion calculation in 2D

    Input:
        u: [theta, rho]
        v: [theta, rho]
    Output:
        dist: distance between u and v

        For u==v, dist=0
        for u[0]==v[0], dist=u[1]-v[1]
        for u[1]==v[1] and dtheta-> 0, dist ->0
        )
    '''
    dtheta=abs(u[0]-v[0])%180
    d1=v[1]-u[1]/np.cos(np.deg2rad(dtheta))
    d2=v[1]/np.cos(np.deg2rad(dtheta))-u[1]
    
    return abs(max(d1,d2))


def plotDendro(dist1, labels, title):

    #https://stackoverflow.com/questions/2982929/plotting-results-of-hierarchical-clustering-ontop-of-a-matrix-of-data-in-python    
    D=dist1
    condensedD = squareform(D)
    
    # Compute and plot first dendrogram.
    fig = pylab.figure(figsize=(8,8))
    
    ax1 = fig.add_axes([0.09,0.1,0.15,0.6])
    ax1.set_title(title)
    Y = sch.linkage(condensedD, method='complete')
    Z1 = sch.dendrogram(Y, labels=labels, orientation='left')
    #ax1.set_xticks()
    #ax1.set_yticks([])
    
    # Compute and plot second dendrogram.
    ax2 = fig.add_axes([0.3,0.76,0.6,0.2])
    Y = sch.linkage(condensedD, method='complete')
    Z2 = sch.dendrogram(Y,  labels=labels)
    #ax2.set_xticks(labels[Z1['leaves']])
    #ax2.set_yticks([])
    
    # Plot distance matrix.
    #add axis(left, bottom, width, height)
    axmatrix = fig.add_axes([0.3,0.1,0.6,0.6])
    idx1 = Z1['leaves']
    idx2 = Z2['leaves']
    D = D[idx1,:]
    D = D[:,idx2]
    im = axmatrix.matshow(D, aspect='auto', origin='lower', cmap=plt.cm.YlGnBu)
    axmatrix.set_xticks([])
    axmatrix.set_yticks([])
    
   # Plot colorbar.
    axcolor = fig.add_axes([0.91,0.1,0.02,0.6])
    pylab.colorbar(im, cax=axcolor)
    fig.show()
    #fig.savefig('dendrogram.png')
    return Z1

def fullTree(model, **kwargs): 
    # Create linkage matrix and then plot the dendrogram
    
    # create the counts of samples under each node
    counts = np.zeros(model.children_.shape[0])
    n_samples = len(model.labels_)
    for i, merge in enumerate(model.children_):
        current_count = 0
        for child_idx in merge:
            if child_idx < n_samples:
                current_count += 1  # leaf node
            else:
                current_count += counts[child_idx - n_samples]
        counts[i] = current_count
    
    linkage_matrix = np.column_stack([model.children_, model.distances_,
                                      counts]).astype(float)
    
    # Plot the corresponding dendrogram
    dendrogram(linkage_matrix, **kwargs)

def CyclicEuclidean(u,v):
    u[0]=u[0]+90
    v[0]=v[0]+90
    dtheta=min( (u[0]-v[0])%180, (v[0]-u[0])%180)
    return np.sqrt( (u[1]-v[1])**2 + dtheta**2)

dikes = pd.read_csv('/home/akh/myprojects/Linking-and-Clustering-Dikes/dikedata/test.csv')
dikes=segLength(dikes)
dikes=giveHashID(dikes)
theta, rho, xc, yc= HT(dikes, xc=0, yc=0)
rho=rho*10000
dikes['theta']=theta
dikes['rho']=rho
dikes= MidtoPerpDistance(dikes, xc, yc)
midist=dikes['PerpOffsetDist'].values

X3D = (np.vstack( (theta, rho, midist) )).T
X2D = (np.vstack( (theta, rho-np.mean(rho))) ).T

fig,ax=plt.subplots(2); plotlines(dikes, 'k', ax[0], ColorBy='label')
DotsHT(fig, ax[1], dikes, ColorBy="label")
ax[0].axis('equal')


# # cluster=AGG(n_clusters=7, affinity=DistortionDistance, linkage="complete").fit(X2D)
# # fullTree(cluster)
# dist1=squareform(pdist(X2D, metric='euclidean'))

dist2=squareform(pdist(X2D, CyclicEuclidean))

# dist25=squareform(pdist(X2D, DistortionDistance2))
#dist2[dist2>1000]=1000
# dist25[dist25>17]=17
# dist3=squareform(pdist(X2D, PolarDistanceAngleModified))


# #dist4=squareform(pdist(X2D, DistortionDistance))

# #Z1=plotDendro(dist1, dikes['label'].values, 'euclidean')

Z1=plotDendro(dist2, dikes['label'].values, 'CyclicEuclidean')

# #Z1=plotDendro(dist3, dikes['label'].values, 'Polar Distance')

# #Z1=plotDendro(dist25, dikes['label'].values, 'Distortion Distance 2')

from clusterMod import HT_AGG_custom

dikes2,clusters=HT_AGG_custom(dikes,10, CyclicEuclidean)



lines, IC=examineClusters(dikes2)
dikes['Labels']=dikes['label']
linesTrue, ic2=examineClusters(dikes)
fig,ax=plt.subplots(2)
plotlines(lines, 'k', ax[0])
DotsHT(fig, ax[1], lines, ColorBy="Linked", cmap=cm.viridis, marker="+")
ax[0].axis('equal')


checkClusterChange(linesTrue, lines)

""" makes graphs of rho vs dtheta and distance contours """
# x,y=np.meshgrid( np.linspace(0,300000), np.linspace(0,300000))
# rho=np.sqrt(x**2+y**2)
# f,ax=plt.subplots(1,3, sharey=True)
# i=0

# for t in [0.5, 1, 2]:   

#     d=abs((1+t)*rho*(1-1/np.cos(np.deg2rad(t))))
#     c=ax[i].imshow(d)
#     ax[i].set_title(r"$\delta \theta $ :" +str(t))
#     CS = ax[i].contour(x,y,d)
#     ax[i].set_xticks([0,100000,200000])
#     ax[i].clabel(CS, inline=True, fontsize=10)
#     ax[i].set_xlabel('X')
    
#     plt.colorbar(c, ax=ax[i], label='Distance')
#     i=i+1
    
# ax[0].set_ylabel('Y')
# ax[0].set_yticks([0,100000,200000])