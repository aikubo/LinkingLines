#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  1 12:49:07 2021

@author: akh
"""
import pandas as pd 
import numpy as np
import matplotlib.pyplot as plt
from fitRectangle import *
#from PrePostProcess import midPoint

def AKH_HT(data, xc=None, yc=None):
    """
    Calculates the Hough Transform of a dataframe of line segments.

    Parameters
    ----------
    data: pandas.Dataframe
        dataframe of the line segments
        must contain ["Xstart", "Ystart", "Xend", "Yend"]
        
    xc, yc: float, optional
        x and y location of the center of the HT.
        If none is given, the center is calculated from the dataframe.

    Returns
    -------
    HT: pandas.Dataframe
        dataframe of the Hough Transform
    """
    
  
    if xc is None or yc is None:
        xc,yc=HT_center(data)

    o_x1 = data['Xstart'].values - xc
    o_x2 = data['Xend'].values - xc
    o_y1 = data['Ystart'].values - yc
    o_y2 = data['Yend'].values - yc

    A = o_y1.astype(float) - o_y2.astype(float) + 0.0000000000000001;
    B = o_x1.astype(float) - o_x2.astype(float)+ 0.000000000000000001;

    m=A/B
    angle=np.arctan(-1/m)
    b1=-1*m*o_x1+o_y1
    rho=b1*np.sin(angle)
    theta=np.rad2deg(angle)
    
    return theta, rho, xc, yc

    
def HT_center(data):
    """
    Finds the center of a dataframe of line segments.
    
    Parameters
    ----------
    df: pandas.Dataframe
        dataframe of the line segments
        must contain ["Xstart", "Ystart", "Xend", "Yend"]

    Returns
    -------
    xc, yc: float
        x and y location of the center of the HT

    """
    xc=np.mean( (data['Xstart'].values+data['Xend'].values)/2)
    yc=np.mean( (data['Ystart'].values+data['Yend'].values)/2)
    
    return xc,yc

def rotateData2(data, rotation_angle):

    """
    Rotates a dataframe of line segments by a given angle.

    Parameters
    ----------
    x,y : numpy.Array
        array of endpoints such as output from fitRectangle/endpoints2

    rotation_angle: float
        angle of rotation in degrees

    xc, yc: float
        x and y location of the center of the HT.
        

    Returns
    -------
    df: pandas.Dataframe
        dataframe of the line segments rotated by the given angle
    """

    

    xc,yc=HT_center(data)

    rotation_angle=float(rotation_angle)
    o_x1 = data['Xstart'].values - xc
    o_x2 = data['Xend'].values - xc
    o_y1 = data['Ystart'].values - yc
    o_y2 = data['Yend'].values - yc
    #print(xc,yc)
    
    x1 = o_x1*np.cos(np.deg2rad(rotation_angle)) - o_y1*np.sin(np.deg2rad(rotation_angle))
    y1 = o_x1*np.sin(np.deg2rad(rotation_angle)) + o_y1*np.cos(np.deg2rad(rotation_angle))
    
    x2 = o_x2*np.cos(np.deg2rad(rotation_angle)) - o_y2*np.sin(np.deg2rad(rotation_angle))
    y2 = o_x2*np.sin(np.deg2rad(rotation_angle)) + o_y2*np.cos(np.deg2rad(rotation_angle))
    
    ang=np.deg2rad(rotation_angle)
    
    dataRotated=data.copy()
    
    dataRotated['Xstart']=x1+xc
    dataRotated['Ystart']=y1+yc
    dataRotated['Xend']=x2+xc
    dataRotated['Yend']=y2+yc
    
    #xcR, ycR=HT_center(dataRotated)
    theta, rho, xc, yc=AKH_HT(dataRotated, xc=xc, yc=yc)
    dataRotated['theta']=theta
    dataRotated['rho']=rho
    
    
    #print(xcR,ycR)

    return dataRotated


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
            



def MidtoPerpDistance(df, xc, yc):
    """
    Find the distance between line segment midpoint and rho line perpendicular 
    intersection. 
    
    Parameters
    ----------
    df: pandas.Dataframe 
        dataframe of the line segments
        must contain ["Xstart", "Ystart", "Xend", "Yend"]
    xc, yc: float 
        x location of HT center 

    Returns
    -------
    df: pandas.Dataframe
    with new columns of ['PerpOffDist', 'PerpIntX', 'PerpIntY']
    """
    
    try: 
        'Xmid' not in df.columns
    except: 
        print("Xmid must be calculate first")
    # if 'Xmid' not in df.columns:
    #     df=midPoint(df)
        
    
    theta,rho,xc,yc=AKH_HT(df, xc, yc)
    intx=rho*np.cos(np.deg2rad(theta))
    inty=rho*np.sin(np.deg2rad(theta))
    df['PerpOffsetDist']=np.sqrt( (df['Xmid'].values-intx)**2 +  (df['Ymid'].values-inty)**2)*np.sign((df['Ymid'].values-inty))
    df['PerpIntX']=intx
    df['PerpIntY']=inty
    
    return df

def moveHTcenter(data, xc=None, yc=None):
    if xc is None and yc is None:
        xc,yc=HT_center(data)
        
    o_x1 = data['Xstart'].values - xc
    o_x2 = data['Xend'].values - xc
    o_y1 = data['Ystart'].values - yc
    o_y2 = data['Yend'].values - yc
    
    df2=pd.DataFrame({'Xstart':o_x1, 'Ystart':o_y1, 'Xend':o_x2, 'Yend':o_y2})
    
    return df2
    
    