#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov  5 10:20:50 2021

@author: akh
"""
from sklearn.cluster import OPTICS, cluster_optics_dbscan
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from plotmod import plotlines, labelcolors, plotbyAngle, BA_HT
from examineMod import examineClusters, plotlabel
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
from plotmod import plotlines, labelcolors
from skimage.segmentation import watershed
from skimage.feature import peak_local_max

from skimage import feature
from htMOD import HT_center, AKH_HT, MidtoPerpDistance
from findRadialCenters import sweepCenters, detectLocalPeaks
from htMOD import AKH_HT as HT
from clusterMod import scalesimple, scaleRobust, HT_AGG
from fitRadialCenters import *
from sklearn.preprocessing import scale, RobustScaler
from clusterMod import analyzeScaling
from sklearn.preprocessing import scale, RobustScaler, PowerTransformer

from sklearn.cluster import DBSCAN
from sklearn.neighbors import NearestNeighbors

lines = pd.read_csv(
    '/home/akh/myprojects/Linking-and-Clustering-Dikes/syntheticRadial_testfindRadial_doubledFalse.csv')
# lines=pd.read_csv("/home/akh/myprojects/Linking-and-Clustering-Dikes/AllCRBLinked.csv")
theta, rho, xc, yc = AKH_HT(lines)
lines['AvgRho'] = rho
lines['AvgTheta'] = theta

lines = MidtoPerpDistance(lines, xc, yc)
midist=lines['PerpOffsetDist'].values

X = (np.vstack( (theta, rho-np.mean(rho), midist-np.mean(midist) )) ).T


#analyzeScaling(theta, rho, midist)


# def cosineRhoDist(x,y):
#     return np.cos(x[0]-x[1])*abs(y[0]-y[1])

# def sineRhoDist(x,y):
#     return np.sin(x[0]-x[1])*abs(y[0]-y[1])

# def SqsineRhoDist(x,y):
#     return np.sqrt(abs(np.sin(x[0]-x[1]))+(y[0]-y[1])**2)

# def RadDist(x,y):
#     xr=0
#     yr=40000
#     xc=18084.226429950962
#     yc=34674.08558893058

#     rhoRadial=(xr-xc)*np.cos(np.deg2rad(x))+(yr-yc)*np.sin(np.deg2rad(x))
#     d=rhoRadial-y
#     return d[1]-d[0]

def CylDistanceAngleModified(u, v):
    # u =[ theta, rho, perpdist]
    
    dist1 = u[0]-v[0]
    dist2 = dist1%180
    deltaT = min(dist1, dist2)
    return np.sqrt(u[1]**2 + v[1]**2 + 2*u[1]*v[1]*np.cos(np.deg2rad(deltaT)) + (u[2]-v[2])**2)

def PolarDistanceAngleModified(u, v):
    # u =[ theta, rho, perpdist]
    
    dist1 = u[0]-v[0]
    dist2 = dist1%180
    deltaT = min(dist1, dist2)
    return np.sqrt(u[1]**2 + v[1]**2 + 2*u[1]*v[1]*np.cos(np.deg2rad(deltaT)))

nbrs = NearestNeighbors(n_neighbors=500, radius=0.4, metric=CylDistanceAngleModified, n_jobs=7).fit(X)
distances, indices = nbrs.kneighbors(X)


# clust=DBSCAN(eps=30000, min_samples=2, metric=CylDistanceAngleModified).fit(X)

# fig,ax=plt.subplots()
# ax.scatter(theta,rho, c=clust.labels_, cmap='turbo')

# labels=clust.labels_

# def dist(u,v):

#     return np.sqrt( (u[1]-v[1])**2*(1+np.sin(u[0]-v[0]+180))


# dikeset=pd.read_csv('/home/akh/myprojects/Linking-and-Clustering-Dikes/dikedata/DikeMountain_SpanishPeaks_3857.csv')
# theta, rho, xc, yc= HT(dikeset)
# trange=2
# rrange=5000

# stdT=np.std(theta)
# stdRho=np.std(rho)
# d2=trange/stdT

# clustering=HT_AGG(dikeset,d2)
# dikeset['Labels']=clustering.labels_
# dikeset['p']=rho
# dikeset['theta']=theta
# Splines,IC=examineClusters(dikeset)

#X=scaleRobust(rho, theta)



# clust = OPTICS(min_samples=2, n_jobs=3,cluster_method='xi', xi=.995, metric=CylDistanceAngleModified).fit(X)

# # # Run the fit

# space = np.arange(len(X))
# reachability = clust.reachability_[clust.ordering_]
# labels = clust.labels_[clust.ordering_]

# plt.figure(figsize=(10, 7))
# G = gridspec.GridSpec(3, 2)
# ax0 = plt.subplot(G[0, :])
# ax1 = plt.subplot(G[1, :])

# ax2 = plt.subplot(G[2, 0])
# ax3 = plt.subplot(G[2, 1])

# ax3.scatter(theta, rho, c=lines['label'], cmap=cm.turbo)
# ax3.set_title("True")

# # Reachability plot
# xi = reachability[1:]/reachability[0:-1]
# xi[xi < 0] = 0
# xi = np.insert(xi, 0, 0, axis=0)

# col, colshort = labelcolors(lines["label"], cm.turbo)
# ax0.set_ylabel("Reachability (epsilon distance), True Labels")
# ax0.set_title("Reachability Plot")
# for klass, color in zip(np.unique(lines["label"]), colshort):
#     Xk = space[lines['label'] == klass]
#     Rk = clust.reachability_[lines['label'] == klass]
#     ax0.plot(Xk, Rk, marker='.', color=color, alpha=0.3)


# colors = ["g.", "r.", "b.", "y.", "c."]
# ax1.set_ylabel("Reachability")
# ax1.set_title("Reachability Plot, True Labels")
# slope=1-xi

# for klass, color in zip(np.unique(labels), colors):
#     Xk = space[labels == klass]
#     Rk = reachability[labels == klass]
#     ax1.scatter(Xk, Rk, c=slope[labels == klass], alpha=0.3, cmap='turbo')


# ax1.plot(space[labels == -1], reachability[labels == -1], "k.", alpha=0.3)
# ax1.set_ylabel("Reachability")
# ax1.set_title("Reachability Plot")

# # OPTICS
# colors = ["g.", "r.", "b.", "y.", "c."]
# for klass, color in zip(range(0, 5), colors):
#     Xk = X[clust.labels_ == klass]
#     ax2.plot(Xk[:, 0], Xk[:, 1], color, alpha=0.3)
# ax2.plot(X[clust.labels_ == -1, 0], X[clust.labels_ == -1, 1], "k+", alpha=0.1)
# ax2.set_title("Automatic Clustering\nOPTICS")


# plt.tight_layout()
# plt.show()

# Centers = RadialFitLabels(lines['AvgTheta'], lines['AvgRho'], clust.labels_, xc, yc)
# trueCenters = np.array([0, 40000, -100000, 200000, 10000,
#                        0, 40000, 100000, -300000, 500000]).reshape(2, 5)
# found = 0
# for k in np.unique(labels):
#     if k == 0:
#         continue
#     c = np.array(Centers['Center'].iloc[k]).reshape(2, 1)
#     centerError = np.sqrt(np.einsum('ij,ij->j', trueCenters-c, trueCenters-c))
#     print(k, np.argmin(centerError), "error", np.min(centerError))

#     if np.min(centerError) < 1000:
#         print("Found True Radial Center")
#         found = found+1

# print("Found", found, "out of 5")
# from sklearn.cluster import DBSCAN