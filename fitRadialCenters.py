#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov  8 12:44:20 2021

@author: akh
"""
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit 
import numpy as np 
import pandas as pd
from plotmod import plotlines, DotsLines
from PrePostProcess import whichForm
from htMOD import HT_center
from clusterMod import CyclicAngleDist
import numpy as np

from scipy.spatial.distance import pdist, squareform


def CenterFunc(theta,xr,yr, xc, yc):
    rhoRadial=(xr-xc)*np.cos(np.deg2rad(theta))+(yr-yc)*np.sin(np.deg2rad(theta))
    return rhoRadial

def RadialFitLabels(lines, labels, xc,yc, plot=False):
    theta=lines['AvgTheta'].values
    rho=lines['AvgRho'].values
    
    if plot:
        fig,ax=plt.subplots(1,2)
        
        xdata=np.linspace(-90,90,200)
        ax[1].scatter(theta,rho,c=labels, label="data", cmap='turbo')
        ax[1].set_ylabel('Rho (m)')
        ax[1].set_ylabel('Theta($^\circ$)')
            
    Centers=pd.DataFrame()
    for i in np.unique(labels):
        t=theta[labels==i]
        p=rho[labels==i]
        popt, pcov=curve_fit( lambda angle, xr,yr: CenterFunc(angle, xr, yr, xc, yc), t,p )
        perr = np.sqrt(np.diag(pcov))
        Centers=Centers.append(pd.DataFrame({ "Label":i,"Center":[popt], "Std Error": [perr]}).astype(object), ignore_index=True)
        
        if plot:
           
            ax.plot(xdata, CenterFunc(xdata, *popt, xc, yc), 'r-',
                     label='fit: xr=%5.3f, yr=%5.3f' % tuple(popt))
            
        plt.legend()
    return Centers
    
def RadialFit(lines,plot=False, ColorBy=None, weight='LayerNumber', ThetaRange=[-90,90]):
    t,r=whichForm(lines)
    theta=lines[t].values
    rho=lines[r].values
    
    m=(theta>ThetaRange[0]) & (theta<ThetaRange[1])
    theta=theta[m]
    rho=rho[m]    
    
    if plot:
        fig,ax=DotsLines(lines, ColorBy=ColorBy, cmap='turbo')
        xdata=np.linspace(-90,90,200)

            
    Centers=pd.DataFrame()
    
    if 'xc' in lines.columns:
        xc=lines['xc'].values[0]
        yc=lines['yc'].values[0]
    else: 
        xc,yc=HT_center(lines)
    popt, pcov=curve_fit( lambda angle, xr,yr: CenterFunc(angle, xr, yr, xc, yc), theta,rho )
    perr = np.sqrt(np.diag(pcov))
    
    residuals=rho-CenterFunc(theta, *popt, xc,yc)
    ss_res=np.sum(residuals**2)
    ss_tot=np.sum( (rho-np.mean(rho))**2)
    r_sq=1-(ss_res/ss_tot)
    Centers=pd.DataFrame({ "Center":[popt], "Std Error": [perr], 'RSq':r_sq})
    if plot:
       
        ax[1].plot(xdata, CenterFunc(xdata, *popt, xc, yc), 'y-',
                 label='fit: xr=%5.3f, yr=%5.3f' % tuple(popt), linewidth=3)
        ax[0].plot( popt[0], popt[1], '*g', markersize=10)
        
    plt.legend()
    return Centers

def RadialAzimuthal(lines, Center):
    
    xdist=lines['Xmid'].values-Center['Center'][0][0]
    ydist=lines['Ymid'].values-Center['Center'][0][1]
    
    rAngle=np.rad2deg(np.arctan2(xdist,ydist))+180

    return rAngle
def CyclicAngle360(u,v):
    return (u-v)%360

def AngleSpacing(rAngle):
    
    SrAngle=np.sort(rAngle)
    
    spacing=np.abs(SrAngle[0:-2]-SrAngle[1:-1])
    
    return np.mean(spacing), np.median(spacing), np.min(spacing), np.max(spacing)

def RipleyRadial(rAngle, plot=False):
    
    """
    Measure of radial dike swarm angle "clumpiness"
    Inpsired by the Ripley K funciton https://stats.stackexchange.com/questions/122668/is-there-a-measure-of-evenness-of-spread
    
    Input: theta, array 
        array of angles in radial swarm
    
    Output: R, float 
    
     """
    tRange=CyclicAngle360( np.max(rAngle), np.min(rAngle))
    theta=rAngle[:,None]
    steps=np.arange(0,tRange,5)
    
    n=len(rAngle)
    l=n/tRange
    
    d=pdist(theta, metric=CyclicAngle360)
    counts,_=np.histogram(d, bins=steps)
    K=l*counts/n
    
    K_Pure=np.ones(len(K))*l/n
    
    K_se=np.sum( (K-K_Pure)**2)/n
    
    if plot:
        fg,ax=plt.subplots()
        ax.plot(steps[:-1], K_Pure, 'k.-')
        ax.plot(steps[:-1], K, 'r.-')
        return K, K_se, fg, ax
    else:
        return K, K_se
   
def NearCenters(lines, Center, tol=10000):
    t,r=whichForm(lines)
    
    xdata=lines[t].values
    xc=lines['xc'].values[0]
    yc=lines['yc'].values[0]
    rhoPerfect=CenterFunc(xdata, Center['Center'][0][0], Center['Center'][0][1], xc, yc)
    
    close=abs(lines[r].values-rhoPerfect)<tol
    Close=lines[close]
    
    print("For Center", Center['Center'][0][0], Center['Center'][0][1])
    rAngle=RadialAzimuthal(Close, Center)
    maxt=np.max(rAngle)
    mint=np.min(rAngle)
    print( "Max angle Range is ", CyclicAngle360(maxt, mint))
    
    print( "Angle Spacing is ")
    spacing=AngleSpacing(rAngle)
    print("Mean: Median: Min: Max:")
    print(spacing)
    
    K,K_se=RipleyRadial(rAngle)
    
    print("Deviation from perfect radial angle spread")
    print("K_se:", K_se)
    print("Perfect K_se would be 0")

    
    return Close
     
     # for i in theta:
     #     temp=np.abs(theta-i)
     #     for s in steps: 
     #         count=np.sum(temp<s)/n
             
     
     
     
    
    
# xdata = np.linspace(-10, 10, 100)
# rho = CenterFunc(xdata, 2500,2500, 0,0)
# rng = np.random.default_rng()
# rho_noise = 100 * rng.normal(size=xdata.size)
# ydata = rho + rho_noise
# plt.plot(xdata, ydata, 'b-', label='data, a=2500, b=2500')

# a,b=curve_fit(CenterFunc, xdata, ydata)
# xc=0
# yc=0
# a, b=curve_fit( lambda t, xr,yr: CenterFunc(t, xr, yr, xc, yc), xdata,ydata )
# xdata = np.linspace(-90, 90, 100)
# plt.plot(xdata, CenterFunc(xdata, *a, xc, yc), 'r-',
#           label='fit: a=%5.3f, b=%5.3f' % tuple(a))

# plt.ylabel('rho')
# plt.xlabel('theta')
# plt.legend()
# plt.show()