#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul  7 11:19:15 2021

@author: akh
"""

import pandas as pd
from htMOD import AKH_HT as HT
from sklearn.preprocessing import scale
import numpy as np 
from clusterMod import *
import matplotlib.pyplot as plt
from plotmod import *
from examineMod import examineClusters, plotlabel
from PrePostProcess import *
from htMOD import rotateData2

""" Unrotated data :"""
dikeset=pd.read_csv('/home/akh/myprojects/Linking-and-Clustering-Dikes/syntheticFrag_large.csv')
true=pd.read_csv('/home/akh/myprojects/Linking-and-Clustering-Dikes/syntheticTrue_large.csv')

def RotationTest(dikeset):
    dikeset=preprocess(dikeset)
    theta, rho, xc, yc= HT(dikeset)
    dikeset['rho']=rho
    dikeset['theta']=theta
    
    trange=2 
    rrange=5000 
    
    stdT=np.std(theta)
    stdRho=np.std(rho)
    d2=trange/stdT
    print(d2)
    
    clustering=HT_AGG(dikeset,d2)
    dikeset['Labels']=clustering.labels_
    lines,IC=examineClusters(dikeset)
    
    
    """ Rotate 20 deg """
    
    rotated20=rotateData2(dikeset, 80)
    theta, rho, xcR, ycR= HT(rotated20)
    rotated20['rho']=rho
    rotated20['theta']=theta
    
    
    stdT=np.std(theta)
    stdRho=np.std(rho)
    d2=trange/stdT
    print(d2)
    
    clustering=HT_AGG(rotated20,d2)
    rotated20['Labels']=clustering.labels_
    linesR,IC=examineClusters(dikeset)

    ## Plot it up 
    fig,ax=plt.subplots(1,3)
    #plotlines(lines, 'grey', ax[0])
    plotlines(dikeset, 'k', ax[0], center=True)
    
    #plotlines(linesR, 'red', ax[0], alpha=0.6)
    plotlines(rotated20, 'red', ax[0], center=True)
    
    ax[1].scatter(dikeset['theta'], dikeset['rho'])
    ax[2].scatter(rotated20['theta'], rotated20['rho'])
    for i in [1,2]:
        ax[i].set_xlim([-90,90])
    if abs(xc-xcR) > 10**-6 or abs(yc-ycR)>10**-6 :
        print ("test 1 failed: centers not equal")
    else:
        print("test1, passed. centers ==")
    
    if np.mean(dikeset['theta']) == np.mean(rotated20['theta']):
        print("means the same")
    else: 
        print("means different. normal:",np.mean(dikeset['theta']), "rotated:", np.mean(rotated20['theta']))
        
    if abs(np.std(dikeset['theta'])-np.std(rotated20['theta']))<10**-3:
        print("STD the same")
    else: 
        print("STD different. normal:",np.std(dikeset['theta']), "rotated:", np.std(rotated20['theta']))
    
    if len(lines) == len(linesR): 
        print("test3 passing, equal number of clusters")
        checkClusterChange(lines, linesR)
    else: 
        print("test3 failed, non-equal clusters")
        checkClusterChange(lines, linesR)
    
    
""" Makes figure of true, fragmented, and linked then does length histograms"""
# fig,ax=plt.subplots(2)
# plotlines(true, 'yellow', ax[0])
# plotlines(lines, 'r', ax[0])
# plotlines(dikeset, 'k', ax[0])

# ax=trueDikeLength( lines, dikeset, 50000, axs=ax[1], secondbar=True)
# ax.text( .60,.90, 'True Length:'+str(20000),transform=ax.transAxes, color='yellow')
# ax.axvline(x=20000, color='yellow')

