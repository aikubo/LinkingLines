#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov  4 14:58:10 2021

@author: akh
"""
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
from htMOD import HT_center, AKH_HT
from findRadialCenters import sweepCenters, detectLocalPeaks
from htMOD import AKH_HT as HT
from clusterMod import scalesimple
from PrePostProcess import *
from htMOD import *

from tsne import tsne

lines=pd.read_csv('/home/akh/myprojects/Linking-and-Clustering-Dikes/syntheticRadial_testfindRadial_doubledFalse.csv')
theta, rho, xc, yc= AKH_HT(lines)
lines['AvgRho']=rho
lines['AvgTheta']=theta

X=scalesimple(rho, theta)

lines=MidtoPerpDistance(lines, xc, yc)
midist=lines['PerpOffsetDist'].values
from plotmod import HT3D

HT3D(theta, rho, lines['PerpOffsetDist'].values, lines['label'].values)
X = (np.vstack( (theta, rho, midist )) ).T

Y=tsne(X, initial_dims=3)
# # And the power (sig_fft is of complex dtype)
# power = np.abs(sig_fft)**2
# N=len(theta)
# T=np.mean(theta[1:-1]-theta[2:])
# # The corresponding frequencies
# sample_freq = np.linspace(0, 1/(2*T), int(N/2))

# # Plot the FFT power
# fig, ax = plt.subplots()
# yf=2.0/N * np.abs(sig_fft[:N//2])
# ax.plot(sample_freq, yf)
# plt.show()