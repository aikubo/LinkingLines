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
from plotmod import plotlines, DotsLinesHT
from htMOD import HT_center
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
    
def RadialFit(lines,plot=False, ColorBy='AvgTheta', weight='LayerNumber'):
    theta=lines['AvgTheta'].values
    rho=lines['AvgRho'].values
    if plot:
        fig,ax=DotsLinesHT(lines, ColorBy=ColorBy, cmap='turbo')
        xdata=np.linspace(-90,90,200)

            
    Centers=pd.DataFrame()
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