#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul  2 11:54:30 2021

@author: akh
"""
import pandas as pd
from htMOD import AKH_HT as HT
from htMOD import HT_center
from sklearn.preprocessing import scale
import numpy as np 
from clusterMod import *
import matplotlib.pyplot as plt
from plotmod import plotlines, labelcolors, plotbyAngle, BA_HT, HThist
from examineMod import examineClusters
import seaborn as sns
from jitteringHTcenter import moveHTcenter, rotateHT


df1=pd.read_csv('/home/akh/myprojects/Linking-and-Clustering-Dikes/synthetic.csv')
df1=pd.read_csv('/home/akh/myprojects/Linking-and-Clustering-Dikes/dikedata/dikeset_ptheta.csv')



def testRotateDikeset(df1, Rot):
    
    dfR=rotateData2(df1, Rot)
    #fig,ax=plt.subplots(1,3)
    
    #plotlines(df1, 'grey', ax[0], center=True) 
    #plotlines(dfR, 'yellow', ax[0], center=True)
    
    t1,r1,xc1, yc1=AKH_HT(df1)
    tR, rR, xcR,ycR=AKH_HT(dfR)
    dfR['theta']=tR
    dfR['rho']=rR
    print("Regular Center", xc1, yc1)
    print("Rotated Center", xcR, ycR)
    if abs(xc1-xcR)<10**-6 and abs(yc1-ycR)<10**-6:
        print("test1 passed, centers equal")
    else:
        print("test1 failed")
        
        
    if np.mean(t1-tR) == -1*Rot: 
        print("test 2 passed", Rot, np.mean(t1-tR))
    else: 
        print("test2 failed", Rot, np.mean(t1-tR))
        fig,ax=plt.subplots(3)
        ax[0].scatter((t1-tR)+Rot, r1)
        mask=((t1-tR+Rot)>0.1)
        lines=df1.iloc[ (t1-tR+Rot)>0.1]
        linesR=dfR.iloc[ (t1-tR+Rot)>0.1]
        print("don't match up", len(lines))
        print(np.mean(t1[~mask]))
        print(np.mean(tR[~mask]))
        
        plotlines(lines, 'grey', ax[1])
        plotlines(linesR, 'r', ax[1])
        ax[1].plot(xc1,yc1, "r*", markersize=20)
        print(np.mean(tR[mask]), np.mean(t1[mask]))
        ax[2].scatter(tR,t1)
        ax[2].scatter(tR[mask],t1[mask])
        ax[2].set_xlabel('Rotated Theta')
        ax[2].set_ylabel('Original Theta')
        lines=df1.iloc[ (t1-tR+Rot)<0.1]
        linesR=dfR.iloc[ (t1-tR+Rot)<0.1]
        print("take bad ones out", Rot, np.mean(lines['theta']-linesR['theta']))
        
    dfRBack=rotateData2(dfR, (-1*Rot))
    
    tRB, rRB, xcRB,ycRB=AKH_HT(dfRBack)
    
    
    if abs(xc1-xcRB)<10**-6 and abs(yc1-ycRB)<10**-6:
        print("test 3 passed, centers equal")
    else:
        print("test 3 failed")

        
        
    if np.mean(t1-tRB) < 10**-6: 
        print("test 4 passed")
    else: 
        print("test 4 failed", np.mean(t1-tRB))
        
    #plotlines(dfRBack, 'r', ax[0], center=True)

    return dfR

dfR=testRotateDikeset(df1, 20)

