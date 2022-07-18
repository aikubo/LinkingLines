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
import scipy.cluster.hierarchy as sch
import labellines
import matplotlib.gridspec as gridspec

np.random.seed(5)

import matplotlib.ticker as mticker

class FixCartesianLabels():
    """ 
        Moves Offset axis tick label to axis label 
        based on: https://stackoverflow.com/questions/45760763/how-to-move-the-y-axis-scale-factor-to-the-position-next-to-the-y-axis-label
    """
    def __init__(self,  ax):
        # self.axis = {"y":ax.yaxis, "x":ax.xaxis}[axis]
        # self.labelx=""
        ax.callbacks.connect('Cartesian Plots Updated', self.update)
        ax.figure.canvas.draw()
        self.update(ax,None)

    def update(self, ax, lim):
        
        for i,l in zip ( [ax.yaxis, ax.xaxis], ['Y', 'X']):
            fmt = i.get_major_formatter()
            i.offsetText.set_visible(False)
            i.set_label_text(l + " ("+ fmt.get_offset()+" m )" )
            # formatter = mticker.ScalarFormatter(useMathText=True)
            # formatter.set_powerlimits((-1,9))
            # i.set_major_formatter(formatter)

class Labeloffset():
    """ 
        Moves Offset axis tick label to axis label 
        based on: https://stackoverflow.com/questions/45760763/how-to-move-the-y-axis-scale-factor-to-the-position-next-to-the-y-axis-label
    """
    def __init__(self,  ax, label="", axis="y"):
        self.axis = {"y":ax.yaxis, "x":ax.xaxis}[axis]
        self.label=label
        ax.callbacks.connect(axis+'lim_changed', self.update)
        ax.figure.canvas.draw()
        self.update(None)

    def update(self, lim):
        fmt = self.axis.get_major_formatter()
        self.axis.offsetText.set_visible(False)
        self.axis.set_label_text(self.label + " ("+ fmt.get_offset()+")" )

from operator import sub
def get_aspect(ax):
    """ Returns plot aspect ratio
    based on: https://stackoverflow.com/questions/41597177/get-aspect-ratio-of-axes
    answer by @Mad Physicist
    
    """
    # Total figure size
    figW, figH = ax.get_figure().get_size_inches()
    # Axis size on figure
    _, _, w, h = ax.get_position().bounds
    # Ratio of display units
    disp_ratio = (figH * h) / (figW * w)
    # Ratio of data units
    # Negative over negative because of the order of subtraction
    data_ratio = sub(*ax.get_ylim()) / sub(*ax.get_xlim())

    return disp_ratio / data_ratio

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

def RGBArraytoHexArray(c):
    return [RGBtoHex(i) for i in c]

def StringColors(values, palette="turbo"):
    if type(values[0]) is not str:
        raise Exception("Must pass list of stings")
    if len(values)<2:
        raise Exception("Must pass list of greater than 2 ")
        
    labels=np.unique(values)
    n_colors=len(labels)
    
    colors=[RGBtoHex(x) for x in sns.color_palette(palette)]
    
    cm=LinearSegmentedColormap.from_list("StringCM", colors, N=n_colors)
    color_idx=np.array([np.where(i==labels)[0][0] for i in values])
    #[hash(label) for label in values] significantly faster than np.where
    
    return color_idx, cm
    
def StringCbar(c, fig, ax, values):
    labels=np.unique(values)
    n_colors=len(labels)
    c_ticks = np.arange(n_colors) #* (n_colors / (n_colors + 1)) + (2 / n_colors)
    #c_ticks=[hash(x) for x in labels]
    cbar = fig.colorbar(c, ticks=c_ticks, ax=ax)
    cbar.ax.set_yticklabels(labels)
    return cbar 

def clustered_lines(xs, ys, theta, length, xmid=None, ymid=None):
    xstart=np.max(xs)
    ystart=np.max(ys)
    
    xend=np.min(xs)
    yend=np.min(ys)
    
    
    
    if xmid is None or ymid is None:
        print('calculating xmid ymid')
        xmid=(xstart+xend)/2
        ymid=(ystart+yend)/2
        
    
    a = np.cos(np.deg2rad(theta))
    b = np.sin(np.deg2rad(theta))
    
    x0 = xmid
    y0 = ymid
    x1 = int(x0 + length/2 * (-b))
    y1 = int(y0 + length/2 * (a))
    x2 = int(x0 - length/2 * (-b))
    y2 = int(y0 - length/2 * (a))
    
    return x1, x2, y1, y2

def clustered_lines2(xs, ys, theta, rho, length, xc, yc):
    xstart=np.max(xs)
    ystart=np.max(ys)
    
    xend=np.min(xs)
    yend=np.min(ys)
    
    xmid=(xstart+xend)/2
    ymid=(ystart+yend)/2

    
    slope = -1/np.tan(np.deg2rad(theta)+0.00000000001)
    b=rho/np.sin((np.deg2rad(theta))+0.00000000001)+yc-xc*slope
    
    x1=xstart
    x2=xend 
    
    y1=x1*slope+b
    y2=x2*slope+b
    
    
    return x1, x2, y1, y2

def clusteredLinesComplete(lines):
    xs,ys=endpoints2(lines)
    avgtheta=np.mean(lines['theta'].values)
    rotation_angle=-1*avgtheta #+20
    rotatedLines=rotateData2(lines, rotation_angle)
    #print(avgtheta, rotation_angle)
    w,l=fit_Rec(rotatedLines, xc, yc)
    xstart=np.max(xs)
    ystart=np.max(ys)
    
    xend=np.min(xs)
    yend=np.min(ys)
    
    xmid=(xstart+xend)/2
    ymid=(ystart+yend)/2
    
    
    a = np.cos(np.deg2rad(theta))
    b = np.sin(np.deg2rad(theta))
    
    x0 = xmid
    y0 = ymid
    x1 = int(x0 + length/2 * (-b))
    y1 = int(y0 + length/2 * (a))
    x2 = int(x0 - length/2 * (-b))
    y2 = int(y0 - length/2 * (a))
    
    
    return 

def pltRec(lines, xc, yc, fig=None, ax=None): 
    """
    plots the rectangle defined by the center, a, and the lines
    """
    if fig is None or ax is None:
        fig,ax=plt.subplots()
    xi,yi=endpoints2(lines)
    x0=xc
    y0=yc
    
    ang=-np.mean(np.deg2rad(lines['theta'].values))
    
    
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

    Xmid=(np.max(xs)+np.min(xs))/2
    Ymid=(np.max(ys)+np.min(ys))/2

    xs,ys=unrotate(ang,Xedges,Yedges,x0,y0)
    Xmid, Ymid=unrotate(ang, Xmid, Ymid, x0, y0)
    
    ax.plot(xs, ys, 'k-.', alpha=0.7)

    xstart, xend, ystart, yend=clustered_lines(xi, yi, np.mean(lines['theta'].values), length, xmid=Xmid, ymid=Ymid)
    
    
    ax.plot( [xstart, xend], [ystart, yend], 'g.-')
    for i in range(0,len(lines)):
        ax.plot( [xi[i], xi[i+len(lines)]],  [yi[i], yi[i+len(lines)]], 'r-')
    ax.plot(Xmid, Ymid, 'yp')
    ax.set_aspect('equal')
    
    return fig,ax, length, width

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

def plotlines(data, col, ax, alpha=1, myProj=None, maskar=None, linewidth=1, 
              ColorBy=None, center=False, xc=None, yc=None, extend=False, 
              cmap=cm.turbo, cbarStatus=False):

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
        C=data[ColorBy].values
        
        if type(C[0]) is str:
            C,cmap=StringColors(C)
            col=[RGBtoHex(a) for a in cmap(C)]
            print("in plotlines, colorby is a str")
        
        else:
            norm=Normalize(vmin=min(C), vmax=max(C))

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
    FixCartesianLabels(ax)
    ax.set_aspect('equal')
    #if ColorBy is not None: 
    #   fig.colorbar(label=ColorBy)

def trueDikeLength(lines, dikeset, maxL, Lstep=2000, secondbar=False, axs=None):

    if axs is None:
        fig,axs=plt.subplots() #(2)
    axs.grid(False)
    
    a=axs
    bins=np.arange(1,maxL,2000)
    a.hist(lines['R_Length'], bins=bins, stacked=True,color='tab:blue', label='Linked Segment Length')
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
    
    oldavg=np.median(dikeset['seg_length'])
    oldstd=np.std(dikeset['seg_length'])
    
    newavg=np.median(lines['R_Length'])
    newstd=np.std(lines['R_Length'])
    
    a.text(.60,.80,'Dike median:'+str(round(newavg,0)), transform=a.transAxes)
    a.text( .60, .70, 'Dike STD:'+str(round(newstd,0)),transform=a.transAxes)
    
    a.text( .60,.50 ,'Segment median:'+str(round(oldavg,0)),transform=a.transAxes)
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


    
    
def HThist(lines, rstep, tstep, weights=None,fig=None, ax=None, rbins=None, tbins=None, cmap=cm.cividis):
    t,r=whichForm(lines)
    
    
    
    if fig is None or ax is None: 
        fig,ax=plt.subplots()
    if rbins is None: 
        rbins=np.arange(min(lines[r]), max(lines[r]), rstep)
    if tbins is None:
        tbins=np.arange(-90, 90, tstep)
    h,xe,ye, c=ax.hist2d(lines[t], lines[r], bins=[tbins, rbins], weights=weights, cmap=cmap, norm=mcolors.PowerNorm(0.8))
    fig.colorbar(c, label='Counts', ax=ax)
    ax.set_xlabel('Theta ($^\circ$)')
    ax.set_ylabel('Rho (m)')
    ax.set_title("HT histogram")
    
    
    return fig, ax, [h, xe,ye,c]

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



def DotsHT(fig,ax,lines, ColorBy="R_length", label=None, cmap=cm.turbo, marker='o', title='Hough Transform', CbarLabels=True, StrOn=True):
    
    #plt.rcParams.update({'font.size': 50, 'font.weight': 'normal'})
    #sns.set_context("talk")
    t,r=whichForm(lines)
    
    ax.set_xlabel('Theta ($^\circ$)')
    ax.set_ylabel('Rho (m)')
    
    if ColorBy==None:
        c='grey'
    elif type(lines[ColorBy].values[0]) is str:
        c,cmap=StringColors(lines[ColorBy].values)
        print("colorby is string")
        
    else:
        c=lines[ColorBy].values
        
        
    
    #ax[1], h2=HThist(lines['AvgRho'], lines['AvgTheta'], rstep, tstep, weights=lines['R_Length'], ax=ax[1],rbins=rbins)
    c2=ax.scatter(lines[t], lines[r], c=c, cmap=cmap, edgecolor='black', marker=marker, alpha=0.5)
    #ax.set_title(title)
    
    
 
    
    
    if ColorBy is not None: 
        
        if type(lines[ColorBy].values[0]) is str:
            cbar=StringCbar(c2, fig, ax, lines[ColorBy].values)
        if StrOn:
            cbar=StringCbar(c2, fig, ax, lines[ColorBy].values.astype(str))
            
        else:
            cbar=fig.colorbar(c2, ax=ax)
        
    if ColorBy is not None: 
        if label is None: 
            cbar.set_label(ColorBy)
        elif ColorBy is not None: 
            cbar.set_label(label)
            
    if not CbarLabels:
        cbar.ax.set_yticklabels([])
    ax.set_xlim([-90,90])
    plt.tight_layout()
    
    return fig,ax


def DotsLines(lines, ColorBy="seg_length",cmap=cm.turbo, fig=None, ax=None, CbarLabels=True, StrOn=False):
    t,r=whichForm(lines)
    #plt.rcParams.update({'font.size': 50, 'font.weight': 'normal'})
    #sns.set_context("talk")
    if fig is None:
        fig,ax=plt.subplots(1,2)    #lines['StdRho'].mean()*2

    plotlines(lines, 'k', ax[0], ColorBy=ColorBy, cmap=cmap)
   

    #ax[1], h2=HThist(lines['AvgRho'], lines['AvgTheta'], rstep, tstep, weights=lines['R_Length'], ax=ax[1],rbins=rbins)
    DotsHT(fig, ax[1], lines, ColorBy=ColorBy, cmap=cmap,CbarLabels=CbarLabels, StrOn=StrOn)
    #ax[1].set_title('HT')
    
    

    plt.tight_layout()
    
    return fig,ax

def DotsLinesHist(lines, rstep, tstep, cmap1=cm.turbo, cmap2=cm.gray, ColorBy=None):
    t,r=whichForm(lines)
    if ColorBy is None:
        ColorBy=t
    fig,ax=plt.subplots(1,3)    #lines['StdRho'].mean()*2
    plotlines(lines, 'k', ax[0], ColorBy=ColorBy, cmap=cmap1, center=True)
    ax[0].set_title('Cartesian')
    #ax[1], h2=HThist(lines['AvgRho'], lines['AvgTheta'], rstep, tstep, weights=lines['R_Length'], ax=ax[1],rbins=rbins)
    DotsHT(fig, ax[1], lines, ColorBy=ColorBy, cmap=cmap1)
    ax[1].set_title('HT')
    
    
    fig,ax, h=HThist(lines, rstep, tstep, cmap=cmap2, fig=fig, ax=ax[2])

    plt.tight_layout()

    return fig, ax
    
def BA_Density(dikeset,lines): 
    fig,a=plt.subplots(1,2)
    
    h1=densityPlot(dikeset, ax=a[0], fig=fig)
    h1=densityPlot(lines,ax=a[1], fig=fig)
    
def HT3D(lines, ColorBy='PerpOffsetDist'):
    t,r=whichForm(lines)
    X=lines[t].values
    Y=lines[r].values
    Z=lines['PerpOffsetDist'].values
    
    C=lines[ColorBy].values
    
    sns.set(rc={'axes.facecolor':'white', 'figure.facecolor':'white'})
    fig = plt.figure(facecolor='w')
    ax = fig.add_subplot(111, projection='3d')
    
    norm=Normalize(vmin=min(C), vmax=max(C))
    cmap=cm.turbo
    m = cm.ScalarMappable(norm=norm, cmap=cmap)   
    c=m.to_rgba(C)
    
    ax.scatter(X,Y/1000,Z/1000, c=c, cmap=cm.turbo, edgecolor='grey', alpha=0.3)
    ax.set_xlabel("\nTheta ($^\circ$)")
    ax.set_ylabel("\nRho (km)")
    ax.set_zlabel("\nPerp to Mid (km)", linespacing=3.1)
    plt.tight_layout()
    
    return fig, ax
    
def rotateHT3D(fig,ax,name):
    i=0
    for angle in range(0, 270, 10):
        ax.view_init(30, angle)
        plt.draw()
        #plt.pause(.001)
        
        fig.savefig(name+str(i)+".png",dpi=300)
        i=i+1


#plot results of linking algorithm
def plotResults(data):
    
    fig,ax=plt.subplots(1,5)
    
    
    #plot lines
    plotlines(data, 'k', ax[0], alpha=1, ColorBy='AvgTheta')
    
    ax[0].axis('equal')

    #plot histogram of length 
    ax[1].hist(data['R_Length']/1000, bins=np.arange(1, 200, 5))
    ax[1].set_xlabel('Dike Length (km)')

    #plot histogram of width
    ax[2].hist(data['R_Width'], bins=np.arange(100, 2000, 200))
    ax[2].set_xlabel('Width (m)')

    #plot histogram of AvgTheta 
    ax[3].hist(data['AvgTheta'], bins=np.arange(-90,90, 5))
    ax[3].set_xlabel('Theta $^\circ$')

    #plot histogram of AvgRho
    ax[4].hist(data['AvgRho']/1000, bins=np.arange(min(data['AvgRho'])/1000, max(data['AvgRho'])/1000, 20))
    ax[4].set_xlabel('Rho (km)')

# Helper function used for visualization in the following examples
def identify_axes(ax_dict, fontsize=48):
    """
    Helper to identify the Axes in the examples below.

    Draws the label in a large font in the center of the Axes.

    Parameters
    ----------
    ax_dict : dict[str, Axes]
        Mapping between the title / label and the Axes.
    fontsize : int, optional
        How big the label should be.
    """
    kw = dict(ha="center", va="center", fontsize=fontsize, color="darkgrey")
    for k, ax in ax_dict.items():
        ax.text(0.5, 0.5, k, transform=ax.transAxes, **kw)
        
def breakXaxis(xlim, numAxes=1):
    """
    function to break x axis into based in xlim 
    based on matplotlib example 
    https://matplotlib.org/stable/gallery/subplots_axes_and_figures/broken_axis.html

    num axes cannot be greater than 13

    input: 
        xlim: tuple of x limits
        nAxes: number of axes you wish to make with the same breakpoints

    output:
        fig: figure object
        ax: list of axes objects
    """
    # f, axes = plt.subplots(numAxes,2)
    # ax=axes[:,0]
    # ax2=axes[:,1]
    try :
        numAxes>13
    except ValueError:
        print('You can not have a numAxes greater than 13 ')
    
        
    mosaic="""AAAB"""
    ax1labels=["A"]
    ax2labels=["B"]
    from string import ascii_uppercase
    j=2
    if numAxes>1: 
        for i in range((numAxes-1)):
            
            letter1=ascii_uppercase[j]
            ax1labels.append(letter1)
            j=j+1
            letter2=ascii_uppercase[j]
            ax2labels.append(letter2)
            newline="\n"+letter1*3+letter2
            mosaic=mosaic+newline
            j=j+1
            
    print(mosaic)
    f = plt.figure(constrained_layout=True)
    ax_dict = f.subplot_mosaic(mosaic)
    #identify_axes(ax_dict)
    
    
    ax=[ax_dict[i] for i in ax1labels]
    ax2=[ax_dict[i] for i in ax2labels]

    
    d = .5  # proportion of vertical to horizontal extent of the slanted line
    kwargs = dict(marker=[(-1, -d), (1, d)], markersize=12,
              linestyle="none", color='k', mec='k', mew=1, clip_on=False)
    
    for i in range(numAxes):
        if numAxes == 1:

            ax.set_xlim(xlim[0])
            ax2.set_xlim(xlim[1])

            ax.spines['right'].set_visible(False)
            ax2.spines['left'].set_visible(False)
            #ax.yaxis.tick_top()
            ax2.tick_params(labeleft=False, left=False)  # don't put tick labels at the top
            #ax2.yaxis.tick_bottom()
            
            #plots break symbols

            ax.plot([1, 1], [0, 1], transform=ax.transAxes, **kwargs)
            ax2.plot([0, 0], [0, 1], transform=ax2.transAxes, **kwargs)
            continue
        
        ax[i].set_xlim(xlim[0])
        ax2[i].set_xlim(xlim[1])

        ax[i].spines['right'].set_visible(False)
        ax2[i].spines['left'].set_visible(False)
        #ax[i].yaxis.tick_top()
        ax2[i].tick_params(labelleft=False, left=False)  # don't put tick labels at the top
        #ax2[i].yaxis.tick_bottom()
        
        #plots break symbols
        ax[i].plot([1, 1], [0, 1], transform=ax[i].transAxes, **kwargs)
        ax2[i].plot([0, 0], [0, 1], transform=ax2[i].transAxes, **kwargs)

    return f, ax, ax2

def splitData(xlim, x):
    """
    function to split data into two groups based on xlim
    assume only one breakpoint
    """
    x1=[]
    x2=[]
    for i in range(len(x)):
        if x[i]<max(xlim[0]):
            x1.append(x[i])
        else:
            x2.append(x[i])
    return x1, x2

def plotBreak(xlim, x, y, ax1, ax2,marker, **kwargs):
    """
    function to plot breakpoints
    """
    x1, x2=splitData(xlim, x)
    # ax1.plot(x1, y[0:len(x1)],marker,   **kwargs)
    # ax2.plot(x2, y[len(x1):], marker,  **kwargs)
    ax1.plot(x,y, marker, **kwargs)
    ax2.plot(x,y, marker, **kwargs)
    
def NumtoStringCoord(x,y):
    return "("+str(int(x))+","+str(int(y))+")"
    
def plotRadialOver(fig,ax1, ax2,xc,yc,Crange=50000,n=4, step=None, color='gray', color2="red", colorLines=False):
    ys=np.linspace(yc-Crange, yc+Crange,n)
    xs=np.linspace(xc-Crange, xc+Crange,n)
    xr,yr=np.meshgrid( xs, ys)
    cos=np.cos(np.linspace(-np.pi,np.pi))
    sin=np.sin(np.linspace(-np.pi,np.pi))
    r=[(x-xc)*cos+(y-yc)*sin for x,y in zip(xr.flatten(), yr.flatten())]
    i=0
    labelsAll=[str(x) for x in np.arange(0,len(xr.flatten()))] #[ NumtoStringCoord(x, y) for x,y in zip(xr.flatten(), yr.flatten())]
    colors = cm.rainbow(np.linspace(0, 1, len(ys)))
    # ax1.plot(xr, yr, '*', c=color2)
    # for i, label in zip(r, labels):
    #     ax2.plot( np.rad2deg(np.linspace(-np.pi,np.pi)), i, "-.",label=label, c=color, alpha=0.3)
    
    for y,c in zip(ys, colors):
        ax1.plot(xs, [y]*n, "*", c=c)
        
        r=[(x-xc)*cos+(y-yc)*sin for x in xs]
        labels=labelsAll[(i)*n :n*(i+1)]
     
        for rho, label in zip(r, labels):
            ax2.plot( np.rad2deg(np.linspace(-np.pi,np.pi)), rho, "-.",label=label, c=c, alpha=0.7)
        i=i+1

    #labellines.labelLines(ax2.get_lines())
        
    colors = cm.rainbow(np.linspace(0, 1, len(ys)))
    
    return r
def persistancePlot(Z, fig=None,ax=None, log=False):
    
    if fig is None and ax is None: 
        fig, ax=plt.subplots()
        
    Z1 = sch.dendrogram(Z, orientation='left', no_plot=True)
    
    dcoord=np.array(Z1['dcoord'])
    icoord=np.array(Z1['icoord'])
    c=Z1['color_list']
    idx=Z1['leaves']
    
    
    #scaling for persistance
    #a1=(np.max(icoord)+np.max(dcoord))/2
    #a0=(np.min(icoord)+np.min(dcoord))/2
    
    #dcoord=(dcoord-a0)/(a1-a0)
    #icoord=(icoord-a0)/(a1-a0)
    x=np.max(dcoord)
    
    #ax[0].plot([1,1], [x,x], 'k-', linewidth=10)
    p=np.append(dcoord[:,1]-dcoord[:,0], dcoord[:,2]-dcoord[:,3])
    birth=np.array([ dcoord[:,0], dcoord[:,3]])+1
    death=np.array([ dcoord[:,1], dcoord[:,2]])+1
    
    ax.scatter(birth, p, s=p+5, alpha=0.6, edgecolors="k")
    
    if log:
        
        ax.set_yscale('log')
        ax.set_xscale('log')

    #
    ax.set_xlabel('Birth')
    ax.set_ylabel('Persistance')
    #ax.set_xlim((-2,np.max(birth)))
    #ax.set_ylim((-2,np.max(p)))

    return fig, ax, p, birth
from matplotlib.patches import Arc
from matplotlib.transforms import IdentityTransform, TransformedBbox, Bbox


class AngleAnnotation(Arc):
    """
    Draws an arc between two vectors which appears circular in display space.
    """
    def __init__(self, xy, p1, p2, size=75, unit="points", ax=None,
                 text="", textposition="inside", text_kw=None, **kwargs):
        """
        Parameters
        ----------
        xy, p1, p2 : tuple or array of two floats
            Center position and two points. Angle annotation is drawn between
            the two vectors connecting *p1* and *p2* with *xy*, respectively.
            Units are data coordinates.

        size : float
            Diameter of the angle annotation in units specified by *unit*.

        unit : str
            One of the following strings to specify the unit of *size*:

            * "pixels": pixels
            * "points": points, use points instead of pixels to not have a
              dependence on the DPI
            * "axes width", "axes height": relative units of Axes width, height
            * "axes min", "axes max": minimum or maximum of relative Axes
              width, height

        ax : `matplotlib.axes.Axes`
            The Axes to add the angle annotation to.

        text : str
            The text to mark the angle with.

        textposition : {"inside", "outside", "edge"}
            Whether to show the text in- or outside the arc. "edge" can be used
            for custom positions anchored at the arc's edge.

        text_kw : dict
            Dictionary of arguments passed to the Annotation.

        **kwargs
            Further parameters are passed to `matplotlib.patches.Arc`. Use this
            to specify, color, linewidth etc. of the arc.

        """
        self.ax = ax or plt.gca()
        self._xydata = xy  # in data coordinates
        self.vec1 = p1
        self.vec2 = p2
        self.size = size
        self.unit = unit
        self.textposition = textposition

        super().__init__(self._xydata, size, size, angle=0.0,
                         theta1=self.theta1, theta2=self.theta2, **kwargs)

        self.set_transform(IdentityTransform())
        self.ax.add_patch(self)

        self.kw = dict(ha="center", va="center",
                       xycoords=IdentityTransform(),
                       xytext=(0, 0), textcoords="offset points",
                       annotation_clip=True)
        self.kw.update(text_kw or {})
        self.text = ax.annotate(text, xy=self._center, **self.kw)

    def get_size(self):
        factor = 1.
        if self.unit == "points":
            factor = self.ax.figure.dpi / 72.
        elif self.unit[:4] == "axes":
            b = TransformedBbox(Bbox.from_bounds(0, 0, 1, 1),
                                self.ax.transAxes)
            dic = {"max": max(b.width, b.height),
                   "min": min(b.width, b.height),
                   "width": b.width, "height": b.height}
            factor = dic[self.unit[5:]]
        return self.size * factor

    def set_size(self, size):
        self.size = size

    def get_center_in_pixels(self):
        """return center in pixels"""
        return self.ax.transData.transform(self._xydata)

    def set_center(self, xy):
        """set center in data coordinates"""
        self._xydata = xy

    def get_theta(self, vec):
        vec_in_pixels = self.ax.transData.transform(vec) - self._center
        return np.rad2deg(np.arctan2(vec_in_pixels[1], vec_in_pixels[0]))

    def get_theta1(self):
        return self.get_theta(self.vec1)

    def get_theta2(self):
        return self.get_theta(self.vec2)

    def set_theta(self, angle):
        pass

    # Redefine attributes of the Arc to always give values in pixel space
    _center = property(get_center_in_pixels, set_center)
    theta1 = property(get_theta1, set_theta)
    theta2 = property(get_theta2, set_theta)
    width = property(get_size, set_size)
    height = property(get_size, set_size)

    # The following two methods are needed to update the text position.
    def draw(self, renderer):
        self.update_text()
        super().draw(renderer)

    def update_text(self):
        c = self._center
        s = self.get_size()
        angle_span = (self.theta2 - self.theta1) % 360
        angle = np.deg2rad(self.theta1 + angle_span / 2)
        r = s / 2
        if self.textposition == "inside":
            r = s / np.interp(angle_span, [60, 90, 135, 180],
                                          [3.3, 3.5, 3.8, 4])
        self.text.xy = c + r * np.array([np.cos(angle), np.sin(angle)])
        if self.textposition == "outside":
            def R90(a, r, w, h):
                if a < np.arctan(h/2/(r+w/2)):
                    return np.sqrt((r+w/2)**2 + (np.tan(a)*(r+w/2))**2)
                else:
                    c = np.sqrt((w/2)**2+(h/2)**2)
                    T = np.arcsin(c * np.cos(np.pi/2 - a + np.arcsin(h/2/c))/r)
                    xy = r * np.array([np.cos(a + T), np.sin(a + T)])
                    xy += np.array([w/2, h/2])
                    return np.sqrt(np.sum(xy**2))

            def R(a, r, w, h):
                aa = (a % (np.pi/4))*((a % (np.pi/2)) <= np.pi/4) + \
                     (np.pi/4 - (a % (np.pi/4)))*((a % (np.pi/2)) >= np.pi/4)
                return R90(aa, r, *[w, h][::int(np.sign(np.cos(2*a)))])

            bbox = self.text.get_window_extent()
            X = R(angle, r, bbox.width, bbox.height)
            trans = self.ax.figure.dpi_scale_trans.inverted()
            offs = trans.transform(((X-s/2), 0))[0] * 72
            self.text.set_position([offs*np.cos(angle), offs*np.sin(angle)])
            
def DrawHTRays(df, xc=None, yc=None):
    fig,ax=plt.subplots()
    
    if xc is None and yc is None:
        xc,yc=HT_center(data)
        
    t,r=whichForm(df)
    x1=df[r].values*np.cos(np.deg2rad(df[t].values))
    y1=df[r].values*np.sin(np.deg2rad(df[t].values))
    
    for x,y in zip(x1,y1):
        
        ax.arrow( xc,yc, x,y, head_width=10, head_length=20, fc='k', ec='k', length_includes_head=True)
        
    ax.plot(xc,yc, 'r*', markersize=5)
    ax.set_aspect('equal')
    return x1,y1

def dilationPlot(df, binWidth=1700, EWDilation=None, NSDilation=None, **kwargs):

    fig = plt.figure(figsize=(8,8))
    gs = gridspec.GridSpec(3, 3)
    ax_main = plt.subplot(gs[1:3, :2])
    ax_xDist = plt.subplot(gs[0, :2],sharex=ax_main)
    ax_yDist = plt.subplot(gs[1:3, 2],sharey=ax_main)
    
    xs=[ min(df['Xstart'].min(), df['Xend'].min() ), max( df['Xstart'].max(), df['Xend'].max() )]
    ys=[ min(df['Ystart'].min(), df['Yend'].min() ), max( df['Ystart'].max(), df['Yend'].max() )]
    
    if np.ptp(xs) < binWidth or np.ptp(ys) < binWidth:
        binWidth=np.min( [np.ptp(xs)/10, np.ptp(ys)/10])
    
    binx=np.arange(xs[0]-binWidth, xs[1]+binWidth, binWidth)
    biny=np.arange(ys[0]-binWidth, ys[1]+binWidth, binWidth)
    
    if EWDilation is None or NSDilation is None:
        from dilationCalculation import dilation
        EWDilation, NSDilation=dilation(df, binWidth=binWidth, **kwargs)
        
    plotlines(df, 'k', ax_main)
   
    #ax_xDist.plot(binx[:-1], NSDilation[:-1])
    ax_xDist.fill_between(binx,0, NSDilation)
        
        
    #ax_yDist.plot(EWDilation[:-1], biny[:-1])
    ax_yDist.fill_between(EWDilation,0,biny )
    
    ax_xDist.set(ylabel='NS Dilaton (m)')
    ax_yDist.set(xlabel='EW Dilaton (m)')
    ax_xDist.set_ylim([0,np.max(NSDilation)+100])
    ax_yDist.set_xlim([0,np.max(EWDilation)+100])
    ax_yDist.set_ylim([ys[0], ys[1]])
    
    
    
    return EWDilation, NSDilation, fig, [ax_main, ax_xDist, ax_yDist]