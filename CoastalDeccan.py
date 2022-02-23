#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 25 10:25:01 2022

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



dikeset=pd.read_csv("/home/akh/myprojects/Linking-and-Clustering-Dikes/dikedata/Deccan_Coastal.csv")
dikeset=completePreProcess(dikeset)

allCRB={}
fig,ax=plt.subplots(1,2)
xc,yc=HT_center(dikeset)
dikeset= MidtoPerpDistance(dikeset, xc, yc)
midist=dikeset['PerpOffsetDist'].values

theta,rho, xc, yc=HT(dikeset)
dikeset['theta']=theta
dikeset['rho']=rho
col=['#30123BFF', '#1AE4B6FF', '#FABA39FF', '#7A0403FF']
allLinked=pd.DataFrame()

fig,ax=plt.subplots(1,2)
drho=500
dtheta=1

data,clusters, M=HT_AGG_custom(dikeset, dtheta, drho)
lines,IC=examineClusters(data)

plotResults(lines)
