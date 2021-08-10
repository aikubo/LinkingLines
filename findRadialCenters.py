#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug  3 15:36:30 2021

@author: akh
"""

import pandas as pd
from htMOD import AKH_HT as HT
from sklearn.preprocessing import scale
import numpy as np 
from clusterMod import *
import matplotlib.pyplot as plt
from plotmod import plotlines, labelcolors, plotbyAngle, BA_HT, HThist
from examineMod import examineClusters
import seaborn as sns
from jitteringHTcenter import moveHTcenter, rotateHT
from matplotlib import cm
import matplotlib.colors as colors
from itertools import product
    
def findLocations(err,n=2): 
    loc=np.empty((n,2))
    for i in range(n):
        loc[i,:]=np.unravel_index(np.argmax(err, axis=None), err.shape)
        err[int(loc[i,0]), int(loc[i,1])]=0
        
    return loc

def sweepCenters(df, gridM, threshold, xc, yc):
   """ Generates heatmap of radial center in dike data set 
   
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
    print(err.shape)
    for i in range(len(xs)):
        for j in range(len(ys)): 
            #print(i,j)
            rhoRadial=(xs[i]-xc)*np.cos(np.deg2rad(theta))+(ys[j]-yc)*np.sin(np.deg2rad(theta))
            thresholdArray=(threshold*(np.cos(np.deg2rad(theta)) + np.sin(np.deg2rad(theta))))**2
            mask= ((df['AvgRho']-rhoRadial)**2< thresholdArray)
            trange=np.max(theta[mask])-np.min(theta[mask])-2
            err[j,i]=np.sum(mask) #*trange**2
                
    plotlines(df, 'k', ax[0], center=True)
    c=ax[1].pcolor(xr, yr, err, cmap=cm.Reds, shading='auto') #,  norm=colors.PowerNorm(gamma=3))
    cbar=fig.colorbar(c, ax=ax[1])

    # N,M=len(xs), len(ys)
    # Z = np.zeros((N, M))
    # for i, (x,y) in enumerate(product(xs,ys)):
    #     Z[np.unravel_index(i, (N,M))] = err([x,y])
        
        
    # fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
    # surf = ax.plot_surface(xs, ys, Z.T, cmap=cm.coolwarm,
    #                        linewidth=0, antialiased=False)
    # fig.colorbar(surf, shrink=0.5, aspect=5)
        
    return err, xs, ys

def checkOneCenter(df, threshold, xc, yc, xr, yr):
    fig, ax=plt.subplots(1,2)

    theta=df['AvgTheta']

    err=np.empty_like(theta)
    

    rhoRadial=(xs[i]-xc)*np.cos(np.deg2rad(theta))+(ys[j]-yc)*np.sin(np.deg2rad(theta))
    thresholdArray=threshold*(np.cos(np.deg2rad(theta)) + np.sin(np.deg2rad(theta)))
    err= (df['AvgRho']-rhoRadial) #**2< threshold
    
    plotlines(df, 'k', ax[0], center=True)
    ax[0].plot(xr, yr, "y*")
    ax[1].scatter(df['AvgTheta'], df['AvgRho'], c=err)
    #c=ax[1].pcolor(xr, yr, err, cmap=cm.Reds, shading='auto')
    #cbar=fig.colorbar(c, ax=ax[1])
    
    
    return err, xs, ys

import numpy as np
from scipy.ndimage.filters import maximum_filter
from scipy.ndimage.morphology import generate_binary_structure, binary_erosion
import matplotlib.pyplot as pp
import imagepers

# #for some reason I had to reshape. Numpy ignored the shape header.
# paws_data = np.loadtxt("paws.txt").reshape(4,11,14)

# #getting a list of images
# paws = [p.squeeze() for p in np.vsplit(paws_data,4)]


def detectLocalPeaks(err, persistance=20, plot=True):

    
    """ Find local peaks in data (err) 
    This method employs the imagepers library by Stefan Huber Salzburg Univesity 
    of Applied Science 
    https://www.sthu.org/code/codesnippets/imagepers.html
    
    
    Parameters 
    ----------
    
    err : numpy 2d array
        height data for height filtration 
    pesistance :  float, default=20 
        persistance threshold 
    plot : logical, default=True 
        turns on or off plots of persistance diagram and the loci
        
    
    Returns 
    -------
    xlist : list 
        list of x indices in err where local peaks occur
        
    ylist : list 
        list of y indices in err where local peaks occur
        
    """
    g0=imagepers.persistence(err)
    xlist=[]
    ylist=[]
    # print(g0)
    if plot:
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.set_title("Peristence diagram")
        ax.plot([0,100], [0,100], '-', c='grey')
    for i, homclass in enumerate(g0):
        p_birth, bl, pers, p_death = homclass
        if pers <= 1.0:
            continue
        
        if plot:
            x, y = bl, bl-pers
            ax.plot([x], [y], '.', c='b')
            ax.text(x, y+2, str(i+1), color='b')
            
    if plot:
        ax.set_xlabel("Birth level")
        ax.set_ylabel("Death level")
        fig,ax=plt.subplots(1,2)
        ax[0].set_title("Loci of births")
    
    
    for i, homclass in enumerate(g0):
        p_birth, bl, pers, p_death = homclass
        if pers <= persistance:
            continue
        y, x = p_birth
        xlist.append(x)
        ylist.append(y)
        print("local peak at", x,y)
        if plot:
            ax[0].plot(xs[x], ys[y], '.', c='b')
            ax[0].text(xs[x], ys[y]+5000, str(i+1), color='b')
            
            
        
    if plot: 
        ax[0].pcolor(xs, ys, err, cmap=cm.Reds, shading='auto')
        ax[1].set_xlim((0,err.shape[1]))
        ax[1].set_ylim((0,err.shape[0]))
        plt.gca().invert_yaxis()
        

    return xlist, ylist