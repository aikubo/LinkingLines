#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec  7 14:44:27 2021

@author: akh
"""
import pandas as pd
from htMOD import AKH_HT as HT
from sklearn.preprocessing import scale
import numpy as np 
from clusterMod import *
import matplotlib.pyplot as plt
from plotmod import plotlines, labelcolors, plotbyAngle, BA_HT
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
from fitRadialCenters import RadialFit
lines=pd.read_csv('/home/akh/myprojects/Linking-and-Clustering-Dikes/dikedata/allCRB_dikes.csv')
xc,yc=HT_center(lines)
gridM=1000
lines['LayerNumber']=pd.util.hash_array(lines['Layer'])
import seaborn as sns 

sns.set_context('talk')
lines=MidtoPerpDistance(lines, xc, yc)

Centers=RadialFit(lines, plot=True, ColorBy='LayerNumber')

