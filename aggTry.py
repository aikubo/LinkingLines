#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 17 11:23:39 2021

@author: akh
"""
import pandas as pd
from htMOD import AKH_HT as HT
from sklearn.preprocessing import scale
import numpy as np 
from clusterMod import *
import matplotlib.pyplot as plt
from plotmod import plotlines, labelcolors, plotbyAngle, BA_HT
from examineMod import examineClusters

dikeset=pd.read_csv('/home/akh/myprojects/Linking-and-Clustering-Dikes/dikedata/dikeset_ptheta.csv')
theta, rho, xc, yc= HT(dikeset)
dikeset['rho']=rho
dikeset['theta']=theta

trange=2 
rrange=5000 

stdT=np.std(theta)
stdRho=np.std(rho)
d2=trange/stdT

clustering=HT_AGG(dikeset,d2)
dikeset['Labels']=clustering.labels_
lines,IC=examineClusters(dikeset)

#generate a kmz
colorsSegments=labelcolors(dikeset['Labels'])
colorsDikes=labelcolors(lines['Label'])

errorAnalysis(lines)

TopHTSection(lines, 5000, 3)
#fig,ax=plotbyAngle(dikeset, lines, 20)

## Length and Width plots
# f,a=plt.subplots(2)
# a[0].hist(lines['R_Length'], bins=np.arange(0,100000,5000))
# a[0].set_xlabel('Dike Cluster Length')

# a[1].hist(lines['R_Width'],bins=np.arange(0,2000,100))
# a[1].set_xlabel('Dike Cluster Width')
# plt.tight_layout()

# f2,a2=plt.subplots(3)
# a2[0].scatter(lines['R_Length'], lines['R_Width'])
# a2[0].set_xlabel('Length')
# a2[0].set_ylabel('Width')

# a2[2].scatter(lines['Size'], lines['R_Width'])
# a2[2].set_xlabel('Size')
# a2[2].set_ylabel('Width')

# a2[1].scatter(lines['Size'], lines['R_Length'])
# a2[1].set_xlabel('Size')
# a2[1].set_ylabel('Length')
# fig,ax, h1, h2=BA_HT(dikeset, lines, rstep=2000)

plt.tight_layout()


