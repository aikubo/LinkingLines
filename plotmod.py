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
from matplotlib import cm
#from skimage.transform import probabilistic_hough_line as ProbHough
#from sklearn.cluster import DBSCAN
import matplotlib.pyplot as plt

from pyproj import Proj
from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap, Normalize
import seaborn as sns
import matplotlib.colors as mcolors
from htMOD import HT_center
from fitRectangle import *
from PrePostProcess import whichForm

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
    """
    plots the rectangle defined by the center, a, and the lines
    """

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

def plotlines(data, col, ax, alpha=1, myProj=None, maskar=None, linewidth=1, ColorBy=None, center=False, xc=None, yc=None, extend=False):

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
    
    if ColorBy is not None: 
        C=data[ColorBy]
        norm=Normalize(vmin=min(C), vmax=max(C))
        cmap=cm.turbo
        m = cm.ScalarMappable(norm=norm, cmap=cmap)   
        col=m.to_rgba(C)
        
        
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
        
        if extend: 
            m=(lat1-lat2)/(lon1-lon2)
            b=lat1-lon1*m
            
            l=np.sqrt((lat1-lat2)**2+(lon1-lon2)**2)
            LONG=[lon1-l*4, lon2+l*4]
            LAT=np.multiply(LONG,m)+b
        
        #color handling 
        
        if ColorBy is None: 
            if isinstance(col,list):
                colo=col[i]
            elif isinstance(col, str):
                colo=col
        else: 
            colo=col[i]
        
        
        ax.plot(LONG,LAT, c=colo, alpha=alpha, linewidth=linewidth)
        
    if center: 
        if xc is None or yc is None:
            xc,yc=HT_center(data)
        ax.plot(xc,yc, "*r", markeredgecolor="black")
        ax.axis('equal')
    #if ColorBy is not None: 
    #   fig.colorbar(label=ColorBy)

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

def plotbyAngle(lines, AngleBin, absValue=False):
    
    #colorsSegments=labelcolors(dikeset['Labels'],cm.turbo)
    #colorsDikes=labelcolors(lines['Label'],cm.turbo)
    
    bins=int(180/AngleBin)
    start=-90 
    if absValue:
        lines['AvgTheta']=abs( lines['AvgTheta'])
        bins=int(90/AngleBin)
        start=0
    
    if bins > 9: 
        fig,ax=plt.subplots(2, int(bins/2))
        fig2,ax2=plt.subplots(2, int(bins/2))
    else:
        fig,ax=plt.subplots(2,bins)
    
    colors= cm.get_cmap('viridis', bins)
    j=0
    for i in range(bins): 
        stop=start+AngleBin
        #mask1= (dikeset['theta']>start) & (dikeset['theta']<stop)
        mask= (lines['AvgTheta']>start) & (lines['AvgTheta']<stop)
        

        if start >= 0 and bins > 9:
            
            ax1=ax2[0][j]
            ax22=ax2[1][j]
            j=j+1
        else:
            ax1=ax[0][i]
            ax22=ax[1][i]
            
        ax1.hist(lines[mask]['AvgTheta'], bins=30, color=colors(i))
        ax1.set_title(str(start)+"-"+str(stop))
        ax1.set_xlabel('Theta (deg)')
       
        ax22.hist(lines[mask]['AvgRho'],bins=30, color=colors(i))
        ax22.set_xlabel('Rho (m)')
        #plotlines(lines, 'grey', ax2,alpha=0.2)
        
        #plotlines(dikeset, colorsSegments, ax2, maskar=mask1)
        #plotlines(lines, colorsDikes, ax2, alpha=0.4, maskar=mask2)
        
        
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


def DotsHT(fig,ax,lines, ColorBy="R_length", label=None, cmap=cm.turbo, marker='o'):
    
    #plt.rcParams.update({'font.size': 50, 'font.weight': 'normal'})
    #sns.set_context("talk")
    t,r=whichForm(lines)
    print(t,r)
    ax.set_xlabel('Theta ($^\circ$)')
    ax.set_ylabel('Rho (m)')

    #ax[1], h2=HThist(lines['AvgRho'], lines['AvgTheta'], rstep, tstep, weights=lines['R_Length'], ax=ax[1],rbins=rbins)
    c2=ax.scatter(lines[t], lines[r], c=(lines[ColorBy]), cmap=cmap, edgecolor='black', marker=marker)
    ax.set_title('Hough Transform')
    
    
    #cbar=fig.colorbar(c1, ax=ax[0])
    #cbar.set_label('Segment Length (m)')
    
    cbar=fig.colorbar(c2, ax=ax)
    if label is None: 
        cbar.set_label(ColorBy)
    else: 
        cbar.set_label(label)
    
    
    plt.tight_layout()
    
    return fig,ax


def DotsLinesHT(lines, ColorBy="seg_length",cmap=cm.turbo):
    t
    #plt.rcParams.update({'font.size': 50, 'font.weight': 'normal'})
    #sns.set_context("talk")
    fig,ax=plt.subplots(1,2)    #lines['StdRho'].mean()*2
    plotlines(lines, 'k', ax[0], center=True)
    ax[1].set_xlabel('Theta (degrees)')
    ax[1].set_ylabel('Rho (m)')

    #ax[1], h2=HThist(lines['AvgRho'], lines['AvgTheta'], rstep, tstep, weights=lines['R_Length'], ax=ax[1],rbins=rbins)
    c2=ax[1].scatter(lines['AvgTheta'], lines['AvgRho'], c=(lines[ColorBy]), cmap=cmap, edgecolor='black')
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
    
def HT3D(X,Y,Z, C):
    
    sns.set(rc={'axes.facecolor':'white', 'figure.facecolor':'white'})
    fig = plt.figure(facecolor='w')
    ax = fig.add_subplot(111, projection='3d')
    
    norm=Normalize(vmin=min(C), vmax=max(C))
    cmap=cm.turbo
    m = cm.ScalarMappable(norm=norm, cmap=cmap)   
    c=m.to_rgba(C)
    
    ax.scatter(X,Y/1000,Z/1000, c=c, cmap=cm.turbo, edgecolor=c, alpha=0.3)
    ax.set_xlabel("\nTheta ($^\circ$)")
    ax.set_ylabel("\nRho (km)")
    ax.set_zlabel("\nPerp to Mid (km)", linespacing=3.1)
    plt.tight_layout()
    return fig, ax
    
def rotateHT3D(fig,ax,name):
    i=0
    for angle in range(0, 180, 10):
        ax.view_init(30, angle)
        plt.draw()
        #plt.pause(.001)
        
        fig.savefig(name+str(i)+".png",dpi=300)
        i=i+1
