
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

from scipy.ndimage import gaussian_filter

import numpy as np
import matplotlib.pyplot as plt
from scipy import ndimage as ndi
from plotmod import plotlines


def findLocations(err,n=2): 
    loc=np.empty((n,2))
    for i in range(n):
        loc[i,:]=np.unravel_index(np.argmax(err, axis=None), err.shape)
        err[int(loc[i,0]), int(loc[i,1])]=0
        
    return loc

def sweepCenters(df, gridM, threshold, xc, yc, ThetaRangeWeight=0, plot=False, radialSpreadMinimum=5):
     
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
    ThetaRangeWeight: float optional
        weighting in thresholding as a function of the range in theta of lines 
        intersecting with possible radial center
    plot: boolean optional 
        show plots of radial centers and theta range
     
     Returns 
     -------
     err : numpy.ndarray 
         array of the sum of dikes that are above the threshold for each x,y potential center 
     dikes: numpy.ndarray,  dtype=np.ndarray
     thetaRange: numpy.ndarray
         array of min(dikes['theta'])-max(dikes['theta'])
     thetaSTD: numpy.ndarray
         array of np.std(dikes['theta'])
     xs : numpy.ndarray 
         array of x coordinate of radial centers
     ys : numpy.ndarray 
         array of y coodinate of radial centers
         
     """
    
    

    theta=df['AvgTheta'].values

    xs=np.arange( min( min(df['Xstart']), min(df['Xend']))-2000, max( max(df['Xstart']), max(df['Xend']))+2000, gridM)
    ys=np.arange( min( min(df['Ystart']), min(df['Yend']))-2000, max( max(df['Ystart']), max(df['Yend']))+2000,  gridM)
    xr,yr=np.meshgrid( xs, ys)
    err=np.empty_like(xr)
    thetaRange=np.empty_like(xr)
    dikes = np.empty_like(xr, dtype=np.ndarray)
    thetaSTD=np.empty_like(xr)
    print(err.shape)
    dikelabel=np.arange(0,len(df))
    for i in range(len(xs)):
        for j in range(len(ys)): 
            #print(i,j)
            err[j,i]=0
            rhoRadial=(xs[i]-xc)*np.cos(np.deg2rad(theta))+(ys[j]-yc)*np.sin(np.deg2rad(theta))
            thresholdArray=(threshold*(np.cos(np.deg2rad(theta)) + np.sin(np.deg2rad(theta))))**2
            mask= ((df['AvgRho']-rhoRadial)**2< thresholdArray)

            
            if np.sum(mask) < 5: ##less than 5 
                thetaRange[j,i]=0
                dikes[j,i]=np.empty([1])
                thetaRange[j,i]=0
                thetaSTD[j,i]=0
                continue
            trange=(np.max(theta[mask])-np.min(theta[mask])-2)
            tstd=(np.std(theta[mask]))
            thetaRange[j,i]=trange
            thetaSTD[j,i]=tstd
            
            bins=np.arange(min(theta[mask]), max(theta[mask]), 10) 
            counts=len(np.unique(np.digitize(theta[mask],bins)))
            
            if radialSpreadMinimum <= counts: 
                err[j,i]=np.sum(mask)*((counts)**ThetaRangeWeight)
                dikes[j,i]=dikelabel[mask]
            else:
                err[j,i]=0
                dikes[j,i]=np.empty([1])
                
    
    if plot:
        fig, ax=plt.subplots(1,2)
        c=ax[1].pcolor(xr, yr, err, cmap=cm.Reds, shading='auto') #,  norm=colors.PowerNorm(gamma=3))
        c=ax[0].pcolor(xr, yr, err, cmap=cm.Reds, shading='auto')
        plotlines(df, 'b', ax[1], center=True, alpha=0.8)
        cbar=fig.colorbar(c, ax=ax[1])
        cbar.set_label('Intersecting Dikes')
        ax[1].set_xlabel("X (m)")
        ax[0].set_ylabel("Y (m)")
        ax[0].set_xlabel("X (m)")
        
    
        # fig,ax=plt.subplots(1,2)
        # c=ax[0].pcolor(xr,yr, thetaSTD, cmap=cm.Blues, shading='auto')
        # c=ax[1].pcolor(xr,yr, thetaSTD, cmap=cm.Blues, shading='auto')
        # plotlines(df, 'k', ax[1], center=True, alpha=0.6)
        # fig.colorbar(c,ax=ax[1])
        # ax[1].set_title("Theta STD")
        # ax[1].set_xlim([min(xs), max(xs)])
        # ax[1].set_ylim([min(ys), max(ys)])
    # N,M=len(xs), len(ys)
    # Z = np.zeros((N, M))
    # for i, (x,y) in enumerate(product(xs,ys)):
    #     Z[np.unravel_index(i, (N,M))] = err([x,y])
        
        
    # fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
    # surf = ax.plot_surface(xs, ys, Z.T, cmap=cm.coolwarm,
    #                        linewidth=0, antialiased=False)
    # fig.colorbar(surf, shrink=0.5, aspect=5)
        
    return err, dikes, thetaRange, thetaSTD, xs, ys




def detectLocalPeaks(err,dikes,df,gridM, plot=False,  maxArea=2827433388):
    """ 
    after https://scikit-image.org/docs/stable/auto_examples/color_exposure/plot_regional_maxima.html#sphx-glr-auto-examples-color-exposure-plot-regional-maxima-py
    and 
    https://scikit-image.org/docs/stable/auto_examples/applications/plot_thresholding.html#:~:text=Thresholding%20is%20used%20to%20create,Histogram%2Dbased.
    
    """
    
    xs=np.arange( min( min(df['Xstart']), min(df['Xend']))-2000, max( max(df['Xstart']), max(df['Xend']))+2000, gridM)
    ys=np.arange( min( min(df['Ystart']), min(df['Yend']))-2000, max( max(df['Ystart']), max(df['Yend']))+2000,  gridM)
    xr,yr=np.meshgrid( xs, ys)
    
    image=img_as_float(err)
    image=gaussian_filter(image,1) # filter image
    
    seed=np.copy(image)
    seed[1:-1, 1:-1] = image.min()
    mask = image
    
    dilated = reconstruction(seed, mask, method='dilation') #dilate
    #blocksize=101
    
    #if thresholdType =='otsu':
    thresh=threshold_otsu(image-dilated) #otsu (automatic threshold)
    #else: 
    #    thresh=threshold_local(image-dilated, blocksize)
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
    areas=np.empty( len(np.unique(labels)))
    Centers=np.empty( (len(np.unique(labels)), 2))
    
    #rint(np.unique(labels), np.unique(labels,return_counts=True))
    for i in np.unique(labels):
        name="Catagory"+str(i)+"Percent"
        df[name]=''
        mask=(labels==i)
        areas[i]=np.sum(mask)*gridM**2 
        Centers[i,:]=[np.mean(xr[mask]), np.mean(yr[mask])]
        rdikes,counts=np.unique(np.concatenate(dikes[mask]).ravel(), return_counts=True)
        
        totalcounts[rdikes.astype(int),i]=counts
        if areas[i] > maxArea & i>0: 
            print("Warning: Label", i, "exceeded max area")


        df['CenterX']=Centers[i,0]
        df['CenterY']=Centers[i,1]
        
        


    df["CenterX"]=''
    df["CenterY"]=''
    centerData=pd.DataFrame({"CenterX": Centers[:,0], "CenterY": Centers[:,1], "Area": areas})
    for i in range(len(totalcounts)):
        icounts=(totalcounts[i]/areas)/(np.sum(totalcounts[i]/areas))
        j=np.argmax(icounts)
        
        #print(i, icounts[j], j)
        if icounts[j] > 0.60:
            df["MainCatagory"].iloc[i]=j
        else:
            df["MainCatagory"].iloc[i]=0
            j=0
            
        #df["MainCatagory"].iloc[i]=j
        # print(i, df["MainCatagory"].iloc[i], j)
        #totalcounts[i]=totalcounts[i]/np.sum(totalcounts[i])
        df["CenterX"].iloc[i]=Centers[j,0]
        df["CenterX"].iloc[i]=Centers[j,0]
        
        for j in np.unique(labels):
            name="Catagory"+str(j)+"Percent"
            df[name].iloc[i]=(totalcounts[i][j]/(areas[j]))/(np.sum(totalcounts[i]/areas))

    if plot:
        fig, (ax0, ax1, ax2, ax3) = plt.subplots(nrows=1,
                                    ncols=4,
                                    figsize=(8, 2.5),
                                    sharex=True,
                                    sharey=True)
        ax0.imshow(image - dilated, cmap='gray', origin="lower")
        ax0.set_title('image - dilated')
        ax0.axis('off')
        plt.gca().invert_yaxis()
        fig.tight_layout()
        #fig, ax = try_all_threshold(image-dilated, figsize=(10, 8), verbose=False)

        ax1.set_title('threshold- otsu')
        ax1.imshow(spots,cmap=plt.cm.gray, origin="lower")
        plt.gca().invert_yaxis()

        ax2.imshow(spots, cmap=plt.cm.gray, origin="lower")
        ax2.set_title('Overlapping objects')
        ax3.imshow(labels, cmap=plt.cm.nipy_spectral, origin="lower")
        ax3.set_title('Separated objects')
        
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
    
    

    
    return df,labels,totalcounts,centerData


    