#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 10 11:50:44 2021

@author: akh
"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from scipy.ndimage import gaussian_filter
from skimage import data
from skimage import img_as_float
from skimage.morphology import reconstruction
from skimage.filters import threshold_otsu
from plotmod import plotlines, labelcolors, plotbyAngle, BA_HT, DotsHT
from matplotlib.lines import Line2D

Splines=pd.read_csv('/home/akh/myprojects/Linking-and-Clustering-Dikes/dikedata/SpanishSilverLinked0621.csv')

#fig, ax= DotsHT(Splines, ColorBy="Xstart")

from findRadialCenters import sweepCenters, detectLocalPeaks
xc=-11684130.47751338
yc=4503174.613590027
err, dikes, xs, ys=sweepCenters(Splines, 1000, 1000, xc,yc)
#xlist,ylist = detectLocalPeaks(err, persistance=20, plot=False)

labels=ArithmeticError
detectLocalPeaks(err, plot=False)

# image=img_as_float(err)
# image=gaussian_filter(image,1)

# seed=np.copy(image)
# seed[1:-1, 1:-1] = image.min()
# mask = image

# dilated = reconstruction(seed, mask, method='dilation')

# fig, (ax0, ax1, ax2, ax3) = plt.subplots(nrows=1,
#                                     ncols=4,
#                                     figsize=(8, 2.5),
#                                     sharex=True,
#                                     sharey=True)
# """ 
# after https://scikit-image.org/docs/stable/auto_examples/color_exposure/plot_regional_maxima.html#sphx-glr-auto-examples-color-exposure-plot-regional-maxima-py
# and 
# https://scikit-image.org/docs/stable/auto_examples/applications/plot_thresholding.html#:~:text=Thresholding%20is%20used%20to%20create,Histogram%2Dbased.

# """
# ax0.imshow(image, cmap='gray')
# ax0.set_title('original image')
# ax0.axis('off')

# ax1.imshow(dilated, vmin=image.min(), vmax=image.max(), cmap='gray')
# ax1.set_title('dilated')
# ax1.axis('off')

# ax2.imshow(image - dilated, cmap='gray')
# ax2.set_title('image - dilated')
# ax2.axis('off')
# plt.gca().invert_yaxis()
# fig.tight_layout()
# #fig, ax = try_all_threshold(image-dilated, figsize=(10, 8), verbose=False)
# thresh=threshold_otsu(image-dilated)
# 
# ax3.set_title('threshold- otsu')
# ax3.imshow(spots,cmap=plt.cm.gray)
# plt.gca().invert_yaxis()
# indicesx, indicesy=np.meshgrid( np.arange(0,err.shape[1]), np.arange(0,err.shape[0]))


# import numpy as np
# import matplotlib.pyplot as plt
# from scipy import ndimage as ndi

# from skimage.segmentation import watershed
# from skimage.feature import peak_local_max


# # Now we want to separate the two objects in image
# # Generate the markers as local maxima of the distance to the background
# distance = ndi.distance_transform_edt(spots)
# coords = peak_local_max(distance, footprint=np.ones((3, 3)), labels=spots)
# mask = np.zeros(distance.shape, dtype=bool)
# mask[tuple(coords.T)] = True
# markers, _ = ndi.label(spots)
# labels = watershed(-distance, markers, mask=spots)

# fig, axes = plt.subplots(ncols=2, figsize=(9, 3), sharex=True, sharey=True)
# ax = axes.ravel()

# ax[0].imshow(spots, cmap=plt.cm.gray)
# ax[0].set_title('Overlapping objects')
# ax[1].imshow(labels, cmap=plt.cm.nipy_spectral)
# ax[1].set_title('Separated objects')

# for a in ax:
#     a.set_axis_off()

# fig.tight_layout()
# plt.gca().invert_yaxis()
# plt.show()
# Splines['Catagory']=''
# catagories=["Linear", "Radial1", "Radial2"]

    
# col, colshort=labelcolors(Splines["Catagory"], cm.turbo)
# custom_lines = [Line2D([0], [0], color=colshort[0], lw=4),
#                 Line2D([0], [0], color=colshort[1], lw=4),
#                 Line2D([0], [0], color=colshort[2], lw=4)]


# fig,ax=plt.subplots()
# plotlines(Splines, col, ax)
# ax.legend(custom_lines, ["Linear", "SpanishPeaks", "DikeMountain"])