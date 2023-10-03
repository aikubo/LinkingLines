#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 24 10:12:59 2021

@author: akh
"""
import pandas as pd
from htMOD import AKH_HT as HT
from htMOD import transformXstart
from sklearn.preprocessing import scale
import numpy as np 
from clusterMod import *
import matplotlib.pyplot as plt
from plotmod import *
from PrePostProcess import *


SPl=pd.read_csv('/home/akh/myprojects/Linking-and-Clustering-Dikes/dikedata/SpanishSilverLinked0621.csv')
SP=pd.read_csv('/home/akh/myprojects/Linking-and-Clustering-Dikes/dikedata/DikeMountain_Peaks_3857.csv')


DCl=pd.read_csv('/home/akh/myprojects/Linking-and-Clustering-Dikes/dikedata/DeccanCostalLinked0621.csv')
DC=pd.read_csv('/home/akh/myprojects/Linking-and-Clustering-Dikes/dikedata/Deccan_Central.csv')

CRBl=pd.read_csv('/home/akh/myprojects/Linking-and-Clustering-Dikes/dikedata/CRBLinked0621.csv')
CRB=pd.read_csv('/home/akh/myprojects/Linking-and-Clustering-Dikes/dikedata/dikeset_ptheta.csv')

fig,ax=plt.subplots(3,4)

FullAlgoFigure(SP, SPl, ax[0], fig)
FullAlgoFigure(CRB, CRBl, ax[1], fig)
FullAlgoFigure(DC, DCl, ax[2], fig)

fig,ax=plt.subplots(3)

#trueDikeLength(SPl, SP, 100000, secondbar=True, axs=ax[0])
#trueDikeLength(CRBl, CRB, 100000, secondbar=True, axs=ax[1])
#trueDikeLength(DCl, DC, 100000, secondbar=True, axs=ax[2])
ax[0].hist(SPl['AvgTheta'])
ax[1].hist(CRBl['AvgTheta'])
ax[2].hist(DCl['AvgTheta'])
ax[2].set_xlabel("Theta")
ax[0].set_title("Spanish Peaks")
ax[1].set_title("CRB")
ax[2].set_title("Deccan Coastal")

# ax[0].hist(SPl['R_Length'].values/(SPl['KNN2'].values+1), bins=20)
# ax[1].hist(CRBl['R_Length'].values/(CRBl['KNN2'].values+1), bins=20)
# ax[2].hist(DCl['R_Length'].values/(DCl['KNN2'].values+1), bins=20)

# print(np.mean(SPl['R_Length'].values/(SPl['KNN2'].values+1)))
# print(np.mean(DCl['R_Length'].values/(DCl['KNN2'].values+1)))
# print(np.mean(CRBl['R_Length'].values/(CRBl['KNN2'].values+1)))