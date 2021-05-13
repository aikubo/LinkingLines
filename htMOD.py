#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  1 12:49:07 2021

@author: akh
"""
import pandas as pd 
import numpy as np


def HT_center(data):
    xc=np.mean( (data['Xstart'].values+data['Xend'].values)/2)
    yc=np.mean( (data['Ystart'].values+data['Yend'].values)/2)
    
    return xc,yc

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