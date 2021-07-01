#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  1 12:49:07 2021

@author: akh
"""
import pandas as pd 
import numpy as np
import matplotlib.pyplot as plt


def HT_center(data):
    xc=np.mean( (data['Xstart'].values+data['Xend'].values)/2)
    yc=np.mean( (data['Ystart'].values+data['Yend'].values)/2)
    
    return xc,yc

def rotateData2(data, rotation_angle):
    xc,yc=HT_center(data)
    o_x1 = data['Xstart'].values - xc
    o_x2 = data['Xend'].values - xc
    o_y1 = data['Ystart'].values - yc
    o_y2 = data['Yend'].values - yc
    
    
    x1 = o_x1*np.cos(np.deg2rad(rotation_angle)) - o_y1*np.sin(np.deg2rad(rotation_angle))
    y1 = o_x1*np.sin(np.deg2rad(rotation_angle)) + o_y1*np.cos(np.deg2rad(rotation_angle))
    
    x2 = o_x2*np.cos(np.deg2rad(rotation_angle)) - o_y2*np.sin(np.deg2rad(rotation_angle))
    y2 = o_x2*np.sin(np.deg2rad(rotation_angle)) + o_y2*np.cos(np.deg2rad(rotation_angle))
    
    dataRotated=data.copy()
    
    dataRotated['Xstart']=x1
    dataRotated['Ystart']=y1
    dataRotated['Xend']=x2
    dataRotated['Yend']=y2

    return dataRotated

def transformXstart(dikeset, xc=None, yc=None):
    if xc is None or  yc is None:
        xc,yc=HT_center(dikeset)
    dist1= xc-dikeset['Xstart']
    dist2= xc-dikeset['Xend']
    switchXs=(dist1>dist2)
    
    dikeset.loc[switchXs, ['Xstart', 'Xend']]=(dikeset.loc[switchXs, ['Xend', 'Xstart']].values)
    dikeset.loc[switchXs, ['Ystart', 'Yend']]=(dikeset.loc[switchXs, ['Yend', 'Ystart']].values)
    
    return dikeset

def AKH_HT(data, xc=None, yc=None):
  
    if xc is None and yc is None:
        xc,yc=HT_center(data)

    o_x1 = data['Xstart'].values - xc
    o_x2 = data['Xend'].values - xc
    o_y1 = data['Ystart'].values - yc
    o_y2 = data['Yend'].values - yc

    A = o_y1 - o_y2 + 0.00000000001;
    B = o_x1 - o_x2+ 0.00000000001;

    m=A/B
    angle=np.arctan(-1/m)
    b1=-1*m*o_x1+o_y1
    rho=b1*np.sin(angle)
    theta=np.rad2deg(angle)
    
    return theta, rho, xc, yc

def gridSearch(df,dc):
    xc1,yc1=HT_center(df)
    t=0
    f,a=plt.subplots(4)
    for i in [-dc, dc]:
        for j in  [-dc, dc]:
            theta, rho, xc, yc=AKH_HT(df, xc=xc1+i, yc=yc1+j)
            npts=np.sum(np.abs(rho)>10000)
            a[t].scatter(theta,rho)
            title=str(i)+","+str(j)
            a[t].set_title(title)
            print(npts, i, j)
            t=t+1