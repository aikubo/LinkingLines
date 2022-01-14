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
from plotmod import plotlines
from skimage.segmentation import watershed
from skimage.feature import peak_local_max

from skimage import feature
from htMOD import HT_center, AKH_HT, MidtoPerpDistance
from findRadialCenters import sweepCenters, detectLocalPeaks
from htMOD import AKH_HT as HT
from clusterMod import scalesimple, scaleRobust, HT_AGG
from fitRadialCenters import *
from sklearn.cluster import SpectralClustering
from plotmod import DotsHT
from scipy.spatial.distance import mahalanobis
from scipy.spatial.distance import cdist

lines=pd.read_csv('/home/akh/myprojects/Linking-and-Clustering-Dikes/syntheticRadial_testfindRadial_doubledFalse.csv')
# lines=pd.read_csv("/home/akh/myprojects/Linking-and-Clustering-Dikes/AllCRBLinked.csv")
theta, rho, xc, yc= AKH_HT(lines)
lines['AvgRho']=rho
lines['AvgTheta']=theta

lines=MidtoPerpDistance(lines, xc, yc)

X=(np.vstack(( np.deg2rad(theta), rho, lines['PerpOffsetDist'])).T)
y=cdist(X,X,'mahalanobis')

clustering = SpectralClustering(n_clusters=6,
                                assign_labels='kmeans',
                                n_jobs=3,
                                random_state=0,
                                affinity='precomputed').fit(y)

lines['label']=clustering.labels_
fig, ax=plt.subplots()

fig,ax=DotsHT(fig,ax,lines, ColorBy='label')
