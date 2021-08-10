#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  1 12:49:51 2021

@author: akh
"""

import numpy as np 
import math 
import pandas as pd
from scipy import * 

#from skimage.transform import probabilistic_hough_line as ProbHough
#from sklearn.cluster import DBSCAN
import matplotlib.pyplot as plt

from pyproj import Proj
from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
import seaborn as sns
import matplotlib.colors as mcolors
from htMOD import HT_center
from fitRectangle import *

sns.set()
np.random.seed(5)

def RGBtoHex(vals, rgbtype=1):
  """Converts RGB values in a variety of formats to Hex values.

     @param  vals     An RGB/RGBA tuple
     @param  rgbtype  Valid valus are:
                          1 - Inputs are in the range 0 to 1
                        256 - Inputs are in the range 0 to 255

     @return A hex string in the form '#RRGGBB' or '#RRGGBBAA'
"""

  if len(vals)!=3 and len(vals)!=4:
    raise Exception("RGB or RGBA inputs to RGBtoHex must have three or four elements!")
  if rgbtype!=1 and rgbtype!=256:
    raise Exception("rgbtype must be 1 or 256!")

  #Convert from 0-1 RGB/RGBA to 0-255 RGB/RGBA
  if rgbtype==1:
    vals = [255*x for x in vals]

  #Ensure values are rounded integers, convert to hex, and concatenate
  return '#' + ''.join(['{:02X}'.format(int(round(x))) for x in vals])


def pltRec(lines, xc, yc, a): 
    xi,yi=endpoints2(lines)
    x0=xc
    y0=xc
    ang=-np.deg2rad(np.mean(lines['theta']))
    xp, yp= rotateXYShift(ang, xi,yi, x0,y0)
    #plotlines(lines, 'k.-', a)
    
    width=np.ptp(xp)
    length=np.ptp(yp)
    
    # if width>length :
    #     length=width
    #     width=length
    xc=(max(xp)-min(xp))/2 + min(xp)
    yc=(max(yp)-min(yp))/2 + min(yp)
    
    xr=xc+width/2
    xl=xc-width/2
    yu=yc+length/2
    yd=yc-length/2
    xs=np.append(xr,xl)
    ys=np.append(yu,yd)
    
    
    xpi, ypi=unrotate(ang, xp, yp, x0, y0)

    Xedges=np.array([xs[0], xs[0], xs[1], xs[1], xs[0]])
    Yedges=np.array([ys[1], ys[0], ys[0], ys[1], ys[1]])
    # a.plot(Xedges, Yedges, 'r.-')
    
    xs,ys=unrotate(ang,Xedges,Yedges,x0,y0)

    a.plot(xs, ys, 'r.-')
    for i in range(0,len(lines)):
         a.plot( [xi[i], xi[i+len(lines)]],  [yi[i], yi[i+len(lines)]], 'k.-')

def labelcolors(labels, colormap):
    n=len(np.unique(labels))
    c=colormap(np.linspace(0, 1, n))
    colors=[]
    colorsShort=[ RGBtoHex(c[i]) for i in range(n)]
    
    for i in range(len(labels)):
        cloc=np.where(np.unique(labels)==labels.iloc[i])[0][0]
        cval=RGBtoHex(c[cloc])
        colors.append(cval)
        
    return colors, colorsShort

def plotlines(data, col, ax, alpha=1, myProj=None, maskar=None, linewidth=1, center=False, xc=None, yc=None):

    #plots the line segments contained in data[maskar]
    # in the color col 
    # converts UTM to lat/long
    
    # data - data frame with columns 
    # 'Xstart', 'Xend', 'Yend', 'Ystart'
    # which are UTM coordinates 
    
    # maskar - arraylike, type logical mask of len(data)
    #       masking of data you want to plot
    
    # col   - string or RGB tuple 
    #        color you want all lines to be plotted
    
    # ax    - object you want to plot on 
    
    if maskar is not None:
        temp=data.loc[maskar]
        if not isinstance(col,str):
            col=col[maskar]
    else :
        temp=data
    
    for i in range(0,len(temp)):
        x1=temp['Xstart'].iloc[i]
        y1=temp['Ystart'].iloc[i]
        y2=temp['Yend'].iloc[i]
        x2=temp['Xend'].iloc[i]
        
        
        
        
        if myProj is not None:
            lon1, lat1 = myProj(x1, y1, inverse = True)
            lon2, lat2 = myProj(x2, y2, inverse = True)
        else: 
            lon1=x1
            lat1=y1
            lon2=x2
            lat2=y2
        
        LAT = [lat1, lat2]
        LONG = [lon1, lon2]
        
        if isinstance(col,list):
            colo=col[i]
        elif isinstance(col, str):
            colo=col
        
        ax.plot(LONG,LAT, c=colo, alpha=alpha, linewidth=linewidth)
        
    if center: 
        if xc is None or yc is None:
            xc,yc=HT_center(data)
        ax.plot(xc,yc, "*r", markeredgecolor="black")


def trueDikeLength(lines, dikeset, maxL, Lstep=2000, secondbar=False, axs=None):
    if axs is None:
        fig,axs=plt.subplots() #(2)
    axs.grid(False)
    
    a=axs
    bins=np.arange(1,maxL,2000)
    a.hist(lines['R_Length'], bins=bins, density=True, stacked=True,color='tab:blue', label='Linked Segment Length')
    a.set_ylabel('Count of Linked Dikes' ,color='tab:blue')
    a.tick_params(axis='y', labelcolor='tab:blue')
    oldavg=np.mean(dikeset['seg_length'])
    oldstd=np.std(dikeset['seg_length'])
    
    newavg=np.mean(lines['R_Length'])
    newstd=np.std(lines['R_Length'])
    if secondbar:
        a2=a.twinx()
        a2.grid(False)
        #a2=axs[1]
        bins=np.arange(1,maxL,500)
        a2.hist(dikeset['seg_length'], density=True, stacked=True,histtype="step", color='tab:red', label='Single Segment Length')
        a2.set_ylabel('Count of Segments', color='tab:red')
        a2.tick_params(axis='y', labelcolor='tab:red')
        
        
    a.set_xlabel('Length (m)')
    plt.tight_layout()
    
    oldavg=np.mean(dikeset['seg_length'])
    oldstd=np.std(dikeset['seg_length'])
    
    newavg=np.mean(lines['R_Length'])
    newstd=np.std(lines['R_Length'])
    
    a.text(.60,.80,'Dike mean:'+str(round(newavg,0)), transform=a.transAxes)
    a.text( .60, .70, 'Dike STD:'+str(round(newstd,0)),transform=a.transAxes)
    
    a.text( .60,.50 ,'Segment mean:'+str(round(oldavg,0)),transform=a.transAxes)
    a.text( .60,.40, 'Segment std:'+str(round(oldstd,0)),transform=a.transAxes)
    
    return axs

def plotbyAngle(dikeset, lines, AngleBin):
    
    colorsSegments=labelcolors(dikeset['Labels'])
    colorsDikes=labelcolors(lines['Label'])
    
    bins=int(180/AngleBin)
    fig,ax=plt.subplots(2,bins)
    start=-90 
    
    
    for i in range(bins): 
        stop=start+AngleBin
        mask1= (dikeset['theta']>start) & (dikeset['theta']<stop)
        mask2= (lines['AvgTheta']>start) & (lines['AvgTheta']<stop)
        
        ax1=ax[1][i]
        ax1.hist(lines[mask2]['AvgRho'], bins=30)
        
        ax2=ax[0][i]

        plotlines(lines, 'grey', ax2,alpha=0.2)
        
        plotlines(dikeset, colorsSegments, ax2, maskar=mask1)
        plotlines(lines, colorsDikes, ax2, alpha=0.4, maskar=mask2)
        
        
        start=stop
        
    return fig, ax

def HThist(rho,theta, rstep, tstep, weights=None, ax=None, rbins=None, tbins=None):
    if ax is None: 
        fig,ax=plt.subplots()
    if rbins is None: 
        rbins=np.arange(min(rho), max(rho), rstep)
    if tbins is None:
        tbins=np.arange(-90, 90, tstep)
    h=ax.hist2d(theta, rho, bins=[tbins, rbins], weights=weights, cmap=cm.magma, norm=mcolors.PowerNorm(0.8))
    
    
    return ax, h

def LinesPlusHT(dikeset,lines): 
    fig,ax=plt.subplots(1,2)
    plotlines(lines, 'grey', ax[0] )
    plotlines(dikeset, 'k', ax[0],  center=True)
    
    c1=ax[1].scatter(lines['AvgTheta'], lines['AvgRho'], c=lines['Xstart'], cmap=cm.plasma, edgecolor='black')
    ax[1].set_title('Hough Space')
    ax[0].set_title('Cartesian Space')
    ax[1].set_xlabel('Theta (degrees)')
    ax[1].set_ylabel('Rho (m)')
    cbar=fig.colorbar(c1, ax=ax[1])
    cbar.set_label('Segment Length (m)')
    return fig,ax
    
def roseDiagram(df):
    ax=plt.subplot(projection='polar')
    ax.bar(df['AvgTheta'], df['AvgRho'],alpha=0.5)
    
    return ax
    
def FullAlgoFigure(dikeset,lines, ax, fig): 
    if len(ax) <4: 
        print("Error: axis is not array of size 4")
            
    plotlines(dikeset, 'k', ax[0], center=True )
    plotlines(lines, 'grey', ax[2] )
    plotlines(dikeset, 'k', ax[2],  center=True)
    
    c1=ax[1].scatter(dikeset['theta'], dikeset['rho'], c=dikeset['Ystart'], cmap=cm.plasma, edgecolor='black')
    c2=ax[3].scatter(lines['AvgTheta'], lines['AvgRho'], c=lines['Ystart'], cmap=cm.plasma, edgecolor='black')
        #lines['AvgTheta'], lines['AvgRho'], c=lines['Xstart'], cmap=cm.plasma, edgecolor='black')
    ax[0].set_title('Cartesian Space (Raw Data), n='+str(len(dikeset)))
    ax[1].set_title('Hough Space (Raw Data)')
    ax[2].set_title('Cartesian Space (Linked), n='+str(len(lines)))
    ax[3].set_title('Hough Space (Linked)')
    
    ax[1].set_xlabel('Theta (degrees)')
    ax[1].set_ylabel('Rho (m)')
    ax[3].set_xlabel('Theta (degrees)')
    ax[3].set_ylabel('Rho (m)')
    cbar=fig.colorbar(c1, ax=ax[1])
    cbar.set_label('Latitude')
    
    cbar=fig.colorbar(c2, ax=ax[3])
    cbar.set_label('Latitude')
    
    return ax
    
def BA_HT(dikeset,lines,rstep=5000):
    fig,ax=plt.subplots(1,3)
     #lines['StdRho'].mean()*2
    tstep=2
    rbins=np.arange(min(dikeset['rho']), max(dikeset['rho']), rstep)

    #ax[0],h1=HThist(dikeset['rho'], dikeset['theta'],rstep, tstep, weights=dikeset['seg_length'], ax=ax[0], rbins=rbins)
    ax[0],h1=HThist(dikeset['rho'], dikeset['theta'],rstep, tstep, ax=ax[0], rbins=rbins)
    ax[0].set_title('Raw Data')
    ax[0].set_xlabel('Theta (degrees)')
    ax[1].set_xlabel('Theta (degrees)') 
    ax[2].set_xlabel('Theta (degrees)')
    ax[0].set_ylabel('Rho (m)')
    #ax[1], h2=HThist(lines['AvgRho'], lines['AvgTheta'], rstep, tstep, weights=lines['R_Length'], ax=ax[1],rbins=rbins)
    ax[1], h2=HThist(lines['AvgRho'], lines['AvgTheta'], rstep, tstep, ax=ax[1],rbins=rbins)
    ax[1].set_title('Clustered Data')
    fig.colorbar(h1[3], ax=ax[0])
    fig.colorbar(h2[3], ax=ax[1])
    hdiff= h2[0] - h1[0]
    x,y=np.meshgrid(h1[1], h1[2])
    divnorm = mcolors.TwoSlopeNorm(vcenter=0)
    c2=ax[2].pcolormesh(x,y,hdiff.T, cmap=cm.RdBu, norm=divnorm)
    ax[2].set_title('Change')
    fig.colorbar(c2, ax=ax[2])
    
    plt.tight_layout()
    
    return fig,ax, h1, h2

def HT3D(lines):
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    c=ax.scatter(lines['KNN2'], lines['AvgRho'], lines['AvgTheta'] , marker='o', c=np.log(lines['Ystart']), cmap="turbo")
    ax.set_xlabel('Xstart')
    ax.set_ylabel('Rho')
    ax.set_zlabel('Theta')
    cbar=fig.colorbar(c)
    cbar.set_label('YStart point (m)')

    return fig, ax


def DotsHT(lines, ColorBy="seg_length"):
    
    #plt.rcParams.update({'font.size': 50, 'font.weight': 'normal'})
    #sns.set_context("talk")
    fig,ax=plt.subplots(1,2)    #lines['StdRho'].mean()*2
    plotlines(lines, 'k', ax[0], center=True)
    ax[1].set_xlabel('Theta (degrees)')
    ax[1].set_ylabel('Rho (m)')

    #ax[1], h2=HThist(lines['AvgRho'], lines['AvgTheta'], rstep, tstep, weights=lines['R_Length'], ax=ax[1],rbins=rbins)
    c2=ax[1].scatter(lines['AvgTheta'], lines['AvgRho'], c=(lines[ColorBy]), cmap=cm.turbo, edgecolor='black')
    ax[1].set_title('HT')
    
    
    #cbar=fig.colorbar(c1, ax=ax[0])
    #cbar.set_label('Segment Length (m)')
    
    cbar=fig.colorbar(c2, ax=ax[1])
    cbar.set_label(ColorBy)
    
    
    plt.tight_layout()
    
    return fig,ax


def DotsHT2(dikeset,lines, ColorBy="seg_length"):
    
    #plt.rcParams.update({'font.size': 50, 'font.weight': 'normal'})
    #sns.set_context("talk")
    fig,ax=plt.subplots(1,3)    #lines['StdRho'].mean()*2
            
    
    #ax[0],h1=HThist(dikeset['rho'], dikeset['theta'],rstep, tstep, weights=dikeset['seg_length'], ax=ax[0], rbins=rbins)
    c1=ax[0].scatter(dikeset['theta'], dikeset['rho'], c=dikeset["seg_length"], cmap=cm.plasma, edgecolor='black', norm=mcolors.PowerNorm(0.3))
    ax[0].set_title('Raw Data')
    ax[0].set_xlabel('Theta (degrees)')
    ax[1].set_xlabel('Theta (degrees)') 
    ax[0].set_ylabel('Rho (m)')
    
    
    ax[0].set_ylim([min(dikeset['rho'])-1000, max(dikeset['rho'])+1000])
    ax[1].set_ylim([min(dikeset['rho'])-1000, max(dikeset['rho'])+1000])
    #ax[1], h2=HThist(lines['AvgRho'], lines['AvgTheta'], rstep, tstep, weights=lines['R_Length'], ax=ax[1],rbins=rbins)
    c2=ax[1].scatter(lines['AvgTheta'], lines['AvgRho'], c=(lines['KNN2']), cmap=cm.turbo, edgecolor='black')
    ax[1].set_title('Clustered Data')
    
    
    #cbar=fig.colorbar(c1, ax=ax[0])
    #cbar.set_label('Segment Length (m)')
    
    cbar=fig.colorbar(c2, ax=ax[1])
    cbar.set_label('Distance to K Nearest Neighbor (m)')
    
    
    plt.tight_layout()
    
    return fig,ax

def densityPlot(lines, binSize=1000, sampling=10, ax=None, fig=None):
    
    xlist,ylist=allpoints(lines)
    minX=min(xlist)
    minY=min(ylist)
    
    maxX=max(xlist)
    maxY=max(ylist)
    
    binX=np.arange(minX,maxX, binSize)
    binY=np.arange(minY,maxY, binSize)
    
    xlist,ylist=allpoints(lines)
    
    if ax is None: 
        fig,ax=plt.subplots()

    h=ax.hist2d(xlist,ylist, bins=[binX, binY], cmap=cm.magma, norm=mcolors.PowerNorm(0.3) )
    cbar=fig.colorbar(h[3], ax=ax)
    
    return h
    
def BA_Density(dikeset,lines): 
    fig,a=plt.subplots(1,2)
    
    h1=densityPlot(dikeset, ax=a[0], fig=fig)
    h1=densityPlot(lines,ax=a[1], fig=fig)
    
