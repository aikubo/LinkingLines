#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 28 11:24:29 2021

@author: akh
"""

from htMOD import HT_center, AKH_HT, rotateData2
from plotmod import plotlines
import pandas as pd 
import matplotlib.pyplot as plt
import numpy as np 


def rotateHT(dikeset,deg):
    xci,yci=HT_center(dikeset)
    p=2
    fig,ax=plt.subplots(1,4)
    fig2,ax2=plt.subplots(2,3)
    
    plotlines(dikeset, 'k', ax[0])
    theta,rho,xc,yc=AKH_HT(dikeset)
    theta0=theta
    rho0=rho
    ax[0].plot(xc, yc, "*", color='grey' , markeredgecolor="black", markersize=20)
    ax[1].set_title('not Rotated (black)')
    ax[1].scatter(theta,rho, c='grey')
    mrho=max(rho)
    minrho=min(rho)
    colors=['red', 'green', 'blue', 'yellow', 'orange', 'purple', 'magenta', 'lime', 'tab:cyan', 'turquoise']
    
    for i in [deg, 2*deg]:
        dfr=rotateData2(dikeset, i)
        theta,rho,xc,yc=AKH_HT(dfr)
        dist=np.sqrt((theta0-theta)**2+(rho0-rho)**2)
        print("jittering mean:", np.mean(dist), "max:", max(dist))
        ax[p].scatter(theta0,rho0, c='grey', alpha=0.5)
        ax[p].scatter(theta,rho, c=dist, cmap="plasma")
        ax[p].set_title("Rotated "+str(i)+colors[p])
        print(colors[p])
        plotlines(dfr, colors[p], ax[0])
    
        
        ax2[0,p-1].scatter(theta0, dist)
        ax2[0,p-1].set_xlabel("theta")
        ax2[0,p-1].set_ylabel("jitter")
        ax2[1,p-1].scatter((rho0), dist)
        ax2[1,p-1].set_xlabel("rho")
        ax2[1,p-1].set_xlim([-90,90])
        ax2[1,p-1].set_ylabel("jitter")

        p=p+1
        print("                         ")

    for i in range(p-1):
        ax[i+1].set_ylim(minrho,mrho)
        ax[i+1].set_xlim([-90,90])

def moveHTcenter(dikeset,r):
    xci,yci=HT_center(dikeset)
    p=2
    fig,ax=plt.subplots(1,10)
    fig2,ax2=plt.subplots(2,9)
    
    plotlines(dikeset, 'k', ax[0])
    theta,rho,xc,yc=AKH_HT(dikeset)
    theta0=theta
    rho0=rho
    ax[0].plot(xc, yc, "*", color='grey' , markeredgecolor="black", markersize=20)
    ax[1].set_title('grey')
    ax[1].scatter(theta,rho, c='grey')
    mrho=max(rho)
    minrho=min(rho)
    colors=['red', 'green', 'blue', 'yellow', 'orange', 'purple', 'magenta', 'lime', 'tab:cyan', 'turquoise']
    for i in [-r,0, r]:
        for j in [-r,0,r]:
            if i==0 and j==0:
                continue
            
            print(colors[p], i, j, np.sqrt(i**2+j**2))
            theta,rho,xc,yc=AKH_HT(dikeset, xc=xci+i, yc=yci+j)
            if max(rho) > mrho:
                mrho=max(rho)
            
            if min(rho) < minrho:
                minrho=min(rho)
                
            dist=np.sqrt((theta0-theta)**2+(rho0-rho)**2)
            print("jittering mean:", np.mean(dist), "max:", max(dist))
            ax[p].scatter(theta0,rho0, c='grey', alpha=0.5)
            ax[p].scatter(theta,rho, c=dist/np.sqrt(i**2+j**2), cmap="plasma")
            ax[p].set_title(colors[p])
            ax2[0,p-1].scatter(theta, dist)
            ax2[0,p-1].set_xlabel("theta")
            ax2[0,p-1].set_ylabel("jitter (m)")
            ax2[1,p-1].scatter((rho), dist- np.sqrt(i**2+j**2))
            ax2[1,p-1].set_xlabel("rho")
            ax2[1,p-1].set_ylabel("jitter (m)")
            
            ax[0].plot(xc, yc, "*", color=colors[p] , markeredgecolor="black", markersize=20)
            p=p+1
            print("                         ")

    for i in range(p-1):
        ax[i+1].set_ylim(minrho,mrho)