#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 17 11:23:39 2021

@author: akh
"""
import pandas as pd
from htMOD import AKH_HT as HT
from sklearn.preprocessing import scale
import numpy as np 
from clusterMod import *
import matplotlib.pyplot as plt
from plotmod import plotlines, labelcolors, plotbyAngle, BA_HT
from examineMod import examineClusters, plotlabel
from PrePostProcess import *
from matplotlib import cm
from findRadialCenters import sweepCenters, detectLocalPeaks

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



'''
dikeset=pd.read_csv('/home/akh/myprojects/Linking-and-Clustering-Dikes/dikedata/dikeset_ptheta.csv')
dikeset=giveID(dikeset)
theta, rho, xc, yc= HT(dikeset)
dikeset['rho']=rho
dikeset['theta']=theta

trange=2 
rrange=5000 

stdT=np.std(theta)
stdRho=np.std(rho)
d2=trange/stdT

clustering=HT_AGG(dikeset,d2)
dikeset['Labels']=clustering.labels_
lines,IC=examineClusters(dikeset)
dikeset.to_csv('/home/akh/myprojects/Linking-and-Clustering-Dikes/dikedata/dikeset_ptheta.csv',index=False)

#generate a kmz
#colorsSegments=labelcolors(dikeset['Labels'])
#colorsDikes=labelcolors(lines['Label'])

#errorAnalysis(lines)

#lines2= writeToQGIS(lines, "CRBLinkedQGIS.csv")

# label=22
# plotlabel(dikeset,label)
#fig, ax= DotsHT(dikeset,lines)
'''
xc=475201.3737670694
yc=4976341.597868606
lines=pd.read_csv('/home/akh/myprojects/Linking-and-Clustering-Dikes/dikedata/CRBLinked0621.csv')
#lines.to_csv('/home/akh/myprojects/Linking-and-Clustering-Dikes/dikedata/CRBLinked0621.csv')

err, dikes,thetaRange, thetaSTD, xs, ys=sweepCenters(lines, 1000, 100, xc,yc)
fig,ax=plt.subplots(1,2)
ax[1].hist(err.flatten())

ax[0].imshow(err)


gridM=5000
maxArea=np.pi*(15000**2)
#lines,labels,counts, Centers=detectLocalPeaks(err,dikes, lines, gridM, maxArea=maxChamber, plot=True)

# df=lines

# plot=True
# xs=np.arange( min( min(df['Xstart']), min(df['Xend']))-2000, max( max(df['Xstart']), max(df['Xend']))+2000, gridM)
# ys=np.arange( min( min(df['Ystart']), min(df['Yend']))-2000, max( max(df['Ystart']), max(df['Yend']))+2000,  gridM)
# xr,yr=np.meshgrid( xs, ys)

# image=img_as_float(err**2)
# image=gaussian_filter(image,1) # filter image

# seed=np.copy(image)
# seed[1:-1, 1:-1] = image.min()
# mask = image

# dilated = reconstruction(seed, mask, method='dilation') #dilate
# #blocksize=101

# #if thresholdType =='otsu':
# thresh=threshold_otsu(image-dilated) #otsu (automatic threshold)
# #else: 
# #    thresh=threshold_local(image-dilated, blocksize)
# indicesx, indicesy=np.meshgrid( np.arange(0,err.shape[1]), np.arange(0,err.shape[0]))
# spots=image-dilated> thresh
# # Now we want to separate the two objects in image
# # Generate the markers as local maxima of the distance to the background
# distance = ndi.distance_transform_edt(spots)
# coords = peak_local_max(distance, footprint=np.ones((3, 3)), labels=spots)
# mask = np.zeros(distance.shape, dtype=bool)
# mask[tuple(coords.T)] = True
# markers, _ = ndi.label(spots)
# labels = watershed(-distance, markers, mask=spots)
# df['MainCatagory']=''

# edges=feature.canny(spots)


# fig,ax=plt.subplots()
# ax.imshow(edges)

# #catagories=["Linear", "Radial1", "Radial2"]
# totalcounts=np.empty( (len(df), len(np.unique(labels))))
# areas=np.empty( len(np.unique(labels)))
# Centers=np.empty( (len(np.unique(labels)), 2))

# #rint(np.unique(labels), np.unique(labels,return_counts=True))
# for i in np.unique(labels):
#     name="Catagory"+str(i)+"Percent"
#     df[name]=''
#     mask=(labels==i)
#     areas[i]=np.sum(mask)*gridM**2 
#     Centers[i,:]=[np.mean(xr[mask]), np.mean(yr[mask])]
#     rdikes,counts=np.unique(np.concatenate(dikes[mask]).ravel(), return_counts=True)
    
#     totalcounts[rdikes.astype(int),i]=counts
#     if areas[i] > maxArea and i>0: 
#         print("Warning: Label", i, "exceeded max area")


#     df['CenterX']=Centers[i,0]
#     df['CenterY']=Centers[i,1]
    
    


# df["CenterX"]=''
# df["CenterY"]=''
# centerData=pd.DataFrame({"CenterX": Centers[:,0], "CenterY": Centers[:,1], "Area": areas})
# for i in range(len(totalcounts)):
#     icounts=(totalcounts[i]/areas)/(np.sum(totalcounts[i]/areas))
#     j=np.argmax(icounts)
    
#     #print(i, icounts[j], j)
#     if icounts[j] > 0.60:
#         df["MainCatagory"].iloc[i]=j
#     else:
#         df["MainCatagory"].iloc[i]=0
#         j=0
        
#     #df["MainCatagory"].iloc[i]=j
#     # print(i, df["MainCatagory"].iloc[i], j)
#     #totalcounts[i]=totalcounts[i]/np.sum(totalcounts[i])
#     df["CenterX"].iloc[i]=Centers[j,0]
#     df["CenterX"].iloc[i]=Centers[j,0]
    
#     for j in np.unique(labels):
#         name="Catagory"+str(j)+"Percent"
#         df[name].iloc[i]=(totalcounts[i][j]/(areas[j]))/(np.sum(totalcounts[i]/areas))

# if plot:
#     fig, (ax0, ax1, ax2, ax3) = plt.subplots(nrows=1,
#                                 ncols=4,
#                                 figsize=(8, 2.5),
#                                 sharex=True,
#                                 sharey=True)
#     ax2.imshow(image - dilated, cmap='gray')
#     ax2.set_title('image - dilated')
#     ax2.axis('off')
#     plt.gca().invert_yaxis()
#     fig.tight_layout()
#     #fig, ax = try_all_threshold(image-dilated, figsize=(10, 8), verbose=False)

#     ax3.set_title('threshold- otsu')
#     ax3.imshow(spots,cmap=plt.cm.gray)
#     plt.gca().invert_yaxis()
#     fig, axes = plt.subplots(ncols=2, figsize=(9, 3), sharex=True, sharey=True)
#     ax = axes.ravel()
    
#     ax[0].imshow(spots, cmap=plt.cm.gray)
#     ax[0].set_title('Overlapping objects')
#     ax[1].imshow(labels, cmap=plt.cm.nipy_spectral)
#     ax[1].set_title('Separated objects')
    
#     for a in ax:
#         a.set_axis_off()
    
#     fig.tight_layout()
#     plt.gca().invert_yaxis()
#     plt.show()



# col, colshort=labelcolors(lines["MainCatagory"], cm.turbo)

# fig,ax=plt.subplots()
# plotlines(lines, col, ax)

# col, colshort=labelcolors(lines["MainCatagory"], cm.viridis)
# # custom_lines = [Line2D([0], [0], color=colshort[0], lw=4),
# #         Line2D([0], [0], color=colshort[1], lw=4),
# #         Line2D([0], [0], color=colshort[2], lw=4)]


# #plotlines(Splines, col, ax[0])
# #ax[0].legend(custom_lines, ["Linear", "SpanishPeaks", "DikeMountain"])

# fig,ax=plt.subplots()
# #plotlines(Splines, col, ax[0])
# #ax[0].legend(custom_lines, ["Linear", "SpanishPeaks", "DikeMountain"])

# sns.set_context('talk')

# for i in np.unique(lines['MainCatagory']):
#     angles=lines[lines['MainCatagory']==i]['AvgTheta']
#     astd=np.std(angles)
#     sortedAngles=np.sort(angles)
#     AverageAngularSpacing=np.mean(np.abs(sortedAngles[0:-1]-sortedAngles[1:]))
#     avgangle=np.mean(angles)
#     ax.scatter(angles, lines[lines['MainCatagory']==i]['AvgRho'], c=colshort[i])
#     ax.set_xlabel('Theta ($^\circ$)')
#     ax.set_ylabel('Rho (m)')
#     #ax[2].scatter(astd, AverageAngularSpacing,c=colshort[i])
#     print("Label:", i)
#     print("Dikes:", np.sum(labels==i))
#     print( "Average Angle:" , avgangle)
#     print("Angular Spacing:", AverageAngularSpacing)
    
# custom_lines = [Line2D([0], [0], color=colshort[0], lw=4),
#         Line2D([0], [0], color=colshort[1], lw=4),
#         Line2D([0], [0], color=colshort[2], lw=4)]


# fig,ax=plt.subplots(1,4)
# #plotlines(Splines, col, ax[0])
# ax[0].legend(custom_lines, ["Linear", "SpanishPeaks", "DikeMountain"])

# for i in np.unique(labels):
#     angles=Splines[Splines['MainCatagory']==i]['AvgTheta']
#     astd=np.std(angles)
#     sortedAngles=np.sort(angles)
#     AverageAngularSpacing=np.mean(np.abs(sortedAngles[0:-1]-sortedAngles[1:]))
#     avgangle=np.mean(angles)
#     ax[1].scatter(avgangle, astd, c=colshort[i])
#     ax[2].scatter(astd, AverageAngularSpacing,c=colshort[i])
#     print("Label:", i)
#     print("Angular Spacing:", AverageAngularSpacing)
    

#TopHTSection(lines, 5000, 3)
#fig,ax=plotbyAngle(dikeset, lines, 20)

## Length and Width plots
# f,a=plt.subplots(2)
# a[0].hist(lines['R_Length'], bins=np.arange(0,100000,5000))
# a[0].set_xlabel('Dike Cluster Length')

# a[1].hist(lines['R_Width'],bins=np.arange(0,2000,100))
# a[1].set_xlabel('Dike Cluster Width')
# plt.tight_layout()

# f2,a2=plt.subplots(3)
# a2[0].scatter(lines['R_Length'], lines['R_Width'])
# a2[0].set_xlabel('Length')
# a2[0].set_ylabel('Width')

# a2[2].scatter(lines['Size'], lines['R_Width'])
# a2[2].set_xlabel('Size')
# a2[2].set_ylabel('Width')

# a2[1].scatter(lines['Size'], lines['R_Length'])
# a2[1].set_xlabel('Size')
# a2[1].set_ylabel('Length')
# fig,ax, h1, h2=BA_HT(dikeset, lines, rstep=2000)

plt.tight_layout()


