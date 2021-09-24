#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug  3 15:36:30 2021

@author: akh
"""

import pandas as pd

import numpy as np 

import matplotlib.pyplot as plt

from jitteringHTcenter import moveHTcenter, rotateHT
from matplotlib import cm
from skimage.morphology import reconstruction

from scipy.ndimage import gaussian_filter

from skimage import img_as_float
from skimage.filters import threshold_otsu
import numpy as np
import matplotlib.pyplot as plt
from scipy import ndimage as ndi
from plotmod import plotlines
from skimage.segmentation import watershed
from skimage.feature import peak_local_max

def findLocations(err,n=2): 
    loc=np.empty((n,2))
    for i in range(n):
        loc[i,:]=np.unravel_index(np.argmax(err, axis=None), err.shape)
        err[int(loc[i,0]), int(loc[i,1])]=0
        
    return loc

def sweepCenters(df, gridM, threshold, xc, yc):
     
    """ 
    Generates heatmap of radial center in dike data set 
    
     Parameters 
     ----------
     
     df : pandas.DataFrame
         dike data set must include columns ['Xstart', 'Ystart', 'Xend', 'Yend', 'AvgTheta', 'AvgRho']
     gridM: float 
         grid spacing for x and y
     threshold :  float
         error threshold between rhoRadial and AvgRho (m)
     xc: float
        x coordinate of center of Hough Transform 
     yc: float 
         y coordinate of center of Hough Transform
         
     
     Returns 
     -------
     err : numpy.ndarray 
         array of the sum of dikes that are above the threshold for each x,y potential center 
     dikes: numpy.ndarray
        
     xs : numpy.ndarray 
         array of x coordinate of radial centers
     ys : numpy.ndarray 
         array of y coodinate of radial centers
         
     """
    
    fig, ax=plt.subplots(1,2)

    theta=df['AvgTheta']

    xs=np.arange( min( min(df['Xstart']), min(df['Xend']))-2000, max( max(df['Xstart']), max(df['Xend']))+2000, gridM)
    ys=np.arange( min( min(df['Ystart']), min(df['Yend']))-2000, max( max(df['Ystart']), max(df['Yend']))+2000,  gridM)
    xr,yr=np.meshgrid( xs, ys)
    err=np.empty_like(xr)
    thetaRange=np.empty_like(xr)
    dikes = np.empty_like(xr, dtype=np.ndarray)
    print(err.shape)
    dikelabel=np.arange(0,len(df))
    for i in range(len(xs)):
        for j in range(len(ys)): 
            #print(i,j)
            rhoRadial=(xs[i]-xc)*np.cos(np.deg2rad(theta))+(ys[j]-yc)*np.sin(np.deg2rad(theta))
            thresholdArray=(threshold*(np.cos(np.deg2rad(theta)) + np.sin(np.deg2rad(theta))))**2
            mask= ((df['AvgRho']-rhoRadial)**2< thresholdArray)
            trange=np.max(theta[mask])-np.min(theta[mask])-2
            err[j,i]=np.sum(mask) #*trange**2
            thetaRange[j,i]=trange
            dikes[j,i]=dikelabel[mask]
    #plotlines(df, 'k', ax[0], center=True, alpha=0.4)
    
    c=ax[1].pcolor(xr, yr, err, cmap=cm.Reds, shading='auto') #,  norm=colors.PowerNorm(gamma=3))
    c=ax[0].pcolor(xr, yr, err, cmap=cm.Reds, shading='auto')
    plotlines(df, 'b', ax[1], center=True, alpha=0.1)
    cbar=fig.colorbar(c, ax=ax[1])

    # N,M=len(xs), len(ys)
    # Z = np.zeros((N, M))
    # for i, (x,y) in enumerate(product(xs,ys)):
    #     Z[np.unravel_index(i, (N,M))] = err([x,y])
        
        
    # fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
    # surf = ax.plot_surface(xs, ys, Z.T, cmap=cm.coolwarm,
    #                        linewidth=0, antialiased=False)
    # fig.colorbar(surf, shrink=0.5, aspect=5)
        
    return err, dikes, thetaRange, xs, ys

# def checkOneCenter(df, threshold, xc, yc, xr, yr):
#     fig, ax=plt.subplots(1,2)

#     theta=df['AvgTheta']

#     err=np.empty_like(theta)
    

#     rhoRadial=(xs[i]-xc)*np.cos(np.deg2rad(theta))+(ys[j]-yc)*np.sin(np.deg2rad(theta))
#     thresholdArray=threshold*(np.cos(np.deg2rad(theta)) + np.sin(np.deg2rad(theta)))
#     err= (df['AvgRho']-rhoRadial) #**2< threshold
    
#     plotlines(df, 'k', ax[0], center=True)
#     ax[0].plot(xr, yr, "y*")
#     ax[1].scatter(df['AvgTheta'], df['AvgRho'], c=err)
#     #c=ax[1].pcolor(xr, yr, err, cmap=cm.Reds, shading='auto')
#     #cbar=fig.colorbar(c, ax=ax[1])
    
    
#     return err, xs, ys


def detectLocalPeaks(err,dikes,df, plot=False):
    """ 
    after https://scikit-image.org/docs/stable/auto_examples/color_exposure/plot_regional_maxima.html#sphx-glr-auto-examples-color-exposure-plot-regional-maxima-py
    and 
    https://scikit-image.org/docs/stable/auto_examples/applications/plot_thresholding.html#:~:text=Thresholding%20is%20used%20to%20create,Histogram%2Dbased.
    
    """
    image=img_as_float(err)
    image=gaussian_filter(image,1) # filter image
    
    seed=np.copy(image)
    seed[1:-1, 1:-1] = image.min()
    mask = image
    
    dilated = reconstruction(seed, mask, method='dilation') #dilate
    thresh=threshold_otsu(image-dilated) #otsu (automatic threshold)
    indicesx, indicesy=np.meshgrid( np.arange(0,err.shape[1]), np.arange(0,err.shape[0]))
    spots=image-dilated> thresh
    # Now we want to separate the two objects in image
    # Generate the markers as local maxima of the distance to the background
    distance = ndi.distance_transform_edt(spots)
    coords = peak_local_max(distance, footprint=np.ones((3, 3)), labels=spots)
    mask = np.zeros(distance.shape, dtype=bool)
    mask[tuple(coords.T)] = True
    markers, _ = ndi.label(spots)
    labels = watershed(-distance, markers, mask=spots)
    df['MainCatagory']=''
    
    #catagories=["Linear", "Radial1", "Radial2"]
    totalcounts=np.empty( (len(df), len(np.unique(labels))))
    areas=np.empty(np.unique(labels))
    for i in np.unique(labels):
        name="Catagory"+str(i)+"Percent"
        df[name]=''
        mask=(labels==i)
        rdikes,counts=np.unique(np.concatenate(dikes[mask]).ravel(), return_counts=True)
        totalcounts[rdikes,i]=counts
        areas[i]=np.sum(mask)

        
    for i in range(len(totalcounts)):
        j=np.argmax(totalcounts[i]/areas)
        df["MainCatagory"].iloc[i]=j
        #totalcounts[i]=totalcounts[i]/np.sum(totalcounts[i])
        
        for j in np.unique(labels):
            name="Catagory"+str(j)+"Percent"
            df[name].iloc[i]=(totalcounts[i][j]/(areas[j]))/(np.sum(totalcounts[i]/areas))

    if plot:
        fig, (ax0, ax1, ax2, ax3) = plt.subplots(nrows=1,
                                    ncols=4,
                                    figsize=(8, 2.5),
                                    sharex=True,
                                    sharey=True)
        ax2.imshow(image - dilated, cmap='gray')
        ax2.set_title('image - dilated')
        ax2.axis('off')
        plt.gca().invert_yaxis()
        fig.tight_layout()
        #fig, ax = try_all_threshold(image-dilated, figsize=(10, 8), verbose=False)

        ax3.set_title('threshold- otsu')
        ax3.imshow(spots,cmap=plt.cm.gray)
        plt.gca().invert_yaxis()
        fig, axes = plt.subplots(ncols=2, figsize=(9, 3), sharex=True, sharey=True)
        ax = axes.ravel()
        
        ax[0].imshow(spots, cmap=plt.cm.gray)
        ax[0].set_title('Overlapping objects')
        ax[1].imshow(labels, cmap=plt.cm.nipy_spectral)
        ax[1].set_title('Separated objects')
        
        for a in ax:
            a.set_axis_off()
        
        fig.tight_layout()
        plt.gca().invert_yaxis()
        plt.show()

        
        # col, colshort=labelcolors(Splines["Catagory"], cm.turbo)
        # custom_lines = [Line2D([0], [0], color=colshort[0], lw=4),
        #         Line2D([0], [0], color=colshort[1], lw=4),
        #         Line2D([0], [0], color=colshort[2], lw=4)]
        
        
        # fig,ax=plt.subplots()
        # plotlines(Splines, col, ax)
        # ax.legend(custom_lines, ["Linear", "SpanishPeaks", "DikeMountain"])
    
    

    
    return df,labels,totalcounts


    
    
    
    # """ Find local peaks in data (err) 
    # This method employs the imagepers library by Stefan Huber Salzburg Univesity 
    # of Applied Science 
    # https://www.sthu.org/code/codesnippets/imagepers.html
    
    
    # Parameters 
    # ----------
    
    # err : numpy 2d array
    #     height data for height filtration 
    # pesistance :  float, default=20 
    #     persistance threshold 
    # plot : logical, default=True 
    #     turns on or off plots of persistance diagram and the loci
        
    
    # Returns 
    # -------
    # xlist : list 
    #     list of x indices in err where local peaks occur
        
    # ylist : list 
    #     list of y indices in err where local peaks occur
        
    # """
    # g0=imagepers.persistence(err)
    # xlist=[]
    # ylist=[]
    # # print(g0)
    # if plot:
    #     fig = plt.figure()
    #     ax = fig.add_subplot(111)
    #     ax.set_title("Peristence diagram")
    #     ax.plot([0,100], [0,100], '-', c='grey')
    # for i, homclass in enumerate(g0):
    #     p_birth, bl, pers, p_death = homclass
    #     if pers <= 1.0:
    #         continue
        
    #     if plot:
    #         x, y = bl, bl-pers
    #         ax.plot([x], [y], '.', c='b')
    #         ax.text(x, y+2, str(i+1), color='b')
            
    # if plot:
    #     ax.set_xlabel("Birth level")
    #     ax.set_ylabel("Death level")
    #     fig,ax=plt.subplots(1,2)
    #     ax[0].set_title("Loci of births")
    
    
    # for i, homclass in enumerate(g0):
    #     p_birth, bl, pers, p_death = homclass
    #     if pers <= persistance:
    #         continue
    #     y, x = p_birth
    #     xlist.append(x)
    #     ylist.append(y)
    #     print("local peak at", x,y)
    #     if plot:
    #         ax[0].plot(xs[x], ys[y], '.', c='b')
    #         ax[0].text(xs[x], ys[y]+5000, str(i+1), color='b')
            
            
        
    # if plot: 
    #     ax[0].pcolor(xs, ys, err, cmap=cm.Reds, shading='auto')
    #     ax[1].set_xlim((0,err.shape[1]))
    #     ax[1].set_ylim((0,err.shape[0]))
    #     plt.gca().invert_yaxis()
        

    # return xlist, ylist



# def drawCircles(err, xlist,ylist, xs, ys, rStart=2 ):
#     xr,yr=np.meshgrid( xs, ys)
    
    
#     def costFunction(err, R, x,y, xs): 
#         ss=np.sum(err[x-R:x+R,y-R:y+R])
#         area= np.pi*(xs[0]-xs[R])**2
        
#         return area/ss
#         # (np.pi*R**2)/sumpoints
        
#     for x,y in zip(xlist,ylist):
#         xval=xs[x]
#         yval=ys[y]
        
#         tol=
#         alpha=0.5
#         for r in range(min(x-err.shape[0], y-err.shape[1]))
#         #while costFuction(err, R, x,y, xs) > tol
#         Rnew=Rnew-alpha*
        
        