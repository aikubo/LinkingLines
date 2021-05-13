#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 20 12:46:10 2021

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

dikeset=pd.read_csv('/home/akh/myprojects/Linking-and-Clustering-Dikes/dikedata/Deccan_Central.csv')

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
dikeset['p']=rho
dikeset['theta']=theta
Deccanlines,IC=examineClusters(dikeset)

#generate a kmz
colorsSegments=labelcolors(dikeset['Labels'])
colorsDikes=labelcolors(Deccanlines['Label'])

#fig,ax=plotbyAngle(dikeset, lines, 20)

#)# Length and Width plots
f,a=plt.subplots(2)
a[0].hist(Deccanlines['R_Length'], bins=np.arange(0,100000,5000))
a[0].set_xlabel('Dike Cluster Length')



a[1].hist(Deccanlines['R_Width'],bins=np.arange(0,2000,100))
a[1].set_xlabel('Dike Cluster Width')
plt.tight_layout()

f2,a2=plt.subplots(3)
a2[0].scatter(Deccanlines['R_Length'], Deccanlines['R_Width'])
a2[0].set_xlabel('Length')
a2[0].set_ylabel('Width')

a2[2].scatter(Deccanlines['Size'], Deccanlines['R_Width'])
a2[2].set_xlabel('Size')
a2[2].set_ylabel('Width')

a2[1].scatter(Deccanlines['Size'], Deccanlines['R_Length'])
a2[1].set_xlabel('Size')
a2[1].set_ylabel('Length')



#fig,ax=plotbyAngle(dikeset, lines, 20)

fig,ax, h1, h2=BA_HT(dikeset, Deccanlines, rstep=15000)
errorAnalysis(Deccanlines)
TopHTSection(lines, 15000, 3)