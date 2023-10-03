#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 12 13:16:04 2021

@author: akh
"""

from sklearn.cluster import OPTICS, cluster_optics_dbscan, DBSCAN
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
from sklearn import metrics
from skimage import img_as_float
from skimage.filters import threshold_otsu, threshold_local
import numpy as np
import matplotlib.pyplot as plt
from scipy import ndimage as ndi
from plotmod import plotlines
from skimage.segmentation import watershed
from skimage.feature import peak_local_max

from skimage import feature
from htMOD import HT_center, AKH_HT
from findRadialCenters import sweepCenters, detectLocalPeaks
from htMOD import AKH_HT as HT
from clusterMod import scalesimple, scaleRobust, HT_AGG
from fitRadialCenters import *
from sklearn.preprocessing import scale, RobustScaler

lines=pd.read_csv('/home/akh/myprojects/Linking-and-Clustering-Dikes/syntheticRadial_testfindRadial_doubledFalse.csv')
# lines=pd.read_csv("/home/akh/myprojects/Linking-and-Clustering-Dikes/AllCRBLinked.csv")
theta, rho, xc, yc= AKH_HT(lines)
lines['AvgRho']=rho
lines['AvgTheta']=theta
df=lines

gridM=5000
xs=np.arange( min( min(df['Xstart']), min(df['Xend']))-2000, max( max(df['Xstart']), max(df['Xend']))+2000, gridM)
ys=np.arange( min( min(df['Ystart']), min(df['Yend']))-2000, max( max(df['Ystart']), max(df['Yend']))+2000,  gridM)
data=np.empty((  len(lines), 2+len(ys)*len(xs)))
X=(np.vstack((theta, rho)).T)
data[:,:2]=X
# xr,yr=np.meshgrid( xs, ys)
# err=np.empty_like(xr)
# thetaRange=np.empty_like(xr)
# dikes = np.empty_like(xr, dtype=np.ndarray)
# thetaSTD=np.empty_like(xr)
# print(err.shape)
# dikelabel=np.arange(0,len(df))
k=2
for i in range(len(xs)):
    for j in range(len(ys)): 

        rhoRadial=(xs[i]-xc)*np.cos(np.deg2rad(theta))+(ys[j]-yc)*np.sin(np.deg2rad(theta))
        data[:,k]=rhoRadial
        k=k+1
        

transformer=RobustScaler().fit(data)
X=transformer.transform(data) 

# #############################################################################
# Compute DBSCAN
db = DBSCAN(eps=.9, min_samples=20).fit(X[:,2:])
core_samples_mask = np.zeros_like(db.labels_, dtype=bool)
core_samples_mask[db.core_sample_indices_] = True
labels = db.labels_

# Number of clusters in labels, ignoring noise if present.
n_clusters_ = len(set(labels)) - (1 if -1 in labels else 0)
n_noise_ = list(labels).count(-1)

print("Estimated number of clusters: %d" % n_clusters_)
print("Estimated number of noise points: %d" % n_noise_)

# #############################################################################

# Plot result
import matplotlib.pyplot as plt

# Black removed and is used for noise instead.
unique_labels = set(labels)
colors = [plt.cm.Spectral(each) for each in np.linspace(0, 1, len(unique_labels))]
for k, col in zip(unique_labels, colors):
    if k == -1:
        # Black used for noise.
        col = [0, 0, 0, 1]

    class_member_mask = labels == k

    xy = X[class_member_mask & core_samples_mask]
    plt.plot(
        xy[:, 0],
        xy[:, 1],
        "o",
        markerfacecolor=tuple(col),
        markeredgecolor="k",
        markersize=14,
    )

    xy = X[class_member_mask & ~core_samples_mask]
    plt.plot(
        xy[:, 0],
        xy[:, 1],
        "o",
        markerfacecolor=tuple(col),
        markeredgecolor="k",
        markersize=6,
    )

plt.title("Estimated number of clusters: %d" % n_clusters_)
plt.show()
# clust = OPTICS(min_samples=30, xi=0.5, n_jobs=3)

# # Run the fit
# clust.fit(X)


# space = np.arange(len(X))
# reachability = clust.reachability_[clust.ordering_]
# labels = clust.labels_[clust.ordering_]

# plt.figure(figsize=(10, 7))
# G = gridspec.GridSpec(3, 2)
# ax0=plt.subplot(G[0,:])
# ax1 = plt.subplot(G[1, :])

# ax2 = plt.subplot(G[2, 0])
# ax3 = plt.subplot(G[2, 1])

# ax3.scatter(theta,rho, c=lines['label'], cmap=cm.turbo)
# ax3.set_title("True")

# # Reachability plot
# xi=reachability[1:]/reachability[0:-1]
# xi[xi<0]=0
# xi = np.insert(xi, 0, 0, axis=0)

# col, colshort=labelcolors(lines["label"], cm.turbo)
# ax0.set_ylabel("Reachability (epsilon distance), True Labels")
# ax0.set_title("Reachability Plot")
# # for klass, color in zip(np.unique(lines["label"]), colshort):
# #     Xk = space[lines['label'].iloc[clust.ordering_] == klass]
# #     Rk = clust.reachability_[lines['label'] == klass]
# #     ax0.plot(Xk, Rk, marker='.', color=color, alpha=0.3)
    
    
# colors = ["g.", "r.", "b.", "y.", "c.", "p."]
# ax1.set_ylabel("Reachability")
# ax1.set_title("Reachability Plot, True Labels")

# for klass, color in zip(np.unique(labels), colors):
#     Xk = space[labels == klass]
#     Rk = reachability[labels == klass]
#     ax1.plot(Xk, Rk, color, alpha=0.3)
    
    
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

# Centers=RadialFit(lines['AvgTheta'], lines['AvgRho'], clust.labels_, xc,yc)
# trueCenters=np.array( [0, 40000,-100000,200000, 10000, 0,40000,100000 ,-300000,500000]).reshape(2,5)
# found=0
# for k in np.unique(labels): 
#     if k == 0:
#         continue 
#     c=np.array(Centers['Center'].iloc[k]).reshape(2,1)
#     centerError=np.sqrt(np.einsum('ij,ij->j', trueCenters-c, trueCenters-c))
#     print(k, np.argmin(centerError), "error", np.min(centerError) )
    
#     if np.min(centerError) < 1000:
#         print( "Found True Radial Center")
#         found=found+1 
        
# print("Found", found, "out of 5")