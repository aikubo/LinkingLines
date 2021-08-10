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
from matplotlib import cm
import matplotlib.animation as animation

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
        

def jitterAnimation(dikeset, frames, name):
    global Xc, Yc, HTvis 
    fig,ax=plt.subplots(1,2, figsize=[12,6])
    
    plotlines(dikeset, 'grey', ax[0])
    theta,rho,xc1,yc1=AKH_HT(dikeset)
    cmap=cm.turbo
    HTvis=ax[1].scatter(theta, rho/max(abs(rho)), c=dikeset['Xstart']-xc1, cmap=cmap, edgecolor='black')
    ax[1].set_ylim([-1, 1])
    ax[1].set_xlabel("Theta")
    ax[1].set_ylabel("Rho/max(abs(rho))")
    ax[0].plot(xc1,yc1, "k*", markersize=15)
    centerplot, =ax[0].plot([],[], 'r*', markeredgecolor='black')
    centerline, =ax[0].plot([],[], 'r--')
    #archimedian spiral
    phi=np.linspace(0,4*np.pi, frames)
    alpha=1500
    Xc=alpha*phi*np.cos(phi)+xc1
    Yc=alpha*phi*np.sin(phi)+yc1
    
    interval=500
    t=np.arange(0,frames)
    
    def UpdateJitter(t):
        global Xc, Yc, HTvis 
        theta,rho,xc,yc=AKH_HT(dikeset, xc=Xc[t], yc=Yc[t])
        HTvis.set_offsets(np.c_[theta, rho/max(abs(rho))])
        centerplot.set_xdata(Xc[t])
        centerplot.set_ydata(Yc[t])
        centerline.set_xdata(Xc[:t+1])
        centerline.set_ydata(Yc[:t+1])
        #ax[1].set_ylim(min(rho), max(rho))
        
    
    ani=animation.FuncAnimation(fig, UpdateJitter, frames=t, interval=interval)
    ani.save(name)
    
# def matchCenter(df, xc, yc, gridM):
#     global xr, yr, HTvis, A, phi,rho, theta
#     fig,ax=plt.subplots(1,2, figsize=[12,6])
#     cmap=cm.Reds_r
#     plotlines(df, 'grey', ax[0])
#     theta,rho,xc1,yc1=AKH_HT(df)
#     HTvis=ax[1].scatter(theta, rho, c=[], cmap=cmap, edgecolor='black')
    
#     xs=np.arange(min(df['Xstart']), max(df['Xstart']), gridM)
#     ys=np.arange(min(df['Ystart']), max(df['Ystart']), gridM)
    
#     xr,yr=np.meshgrid( xs, ys)
    
#     A=np.sqrt( (xc-xr)**2+ (yc-yr)**2)
#     phi=np.arctan( (yc-yr)/ (xc-xr)) 

#     for i in range(len(df)):
#         rho=A*np.sin(df['theta'].iloc[i]-phi)
#         err=(df['rho']-rho)**2
        
#     def updateRadialCenter(t):
#         global xr, yr, HTvis, A, phi,rho, theta
#         ax[0].plot(xr[t], yr[t], 'r*')
        
        
        
#         HTvis.set_array(
    
#     return xr
    
    
    
    