#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 10 13:45:32 2021

@author: akh
"""
import pandas as pd
from htMOD import AKH_HT as HT
from htMOD import transformXstart
from sklearn.preprocessing import scale
import numpy as np 
from clusterMod import *
import matplotlib.pyplot as plt
from plotmod import plotlines, labelcolors, plotbyAngle, BA_HT, DotsHT
from examineMod import examineClusters, plotlabel
from PrePostProcess import *
from fitRectangle import endpoints2

dikeset=pd.read_csv('/home/akh/myprojects/Linking-and-Clustering-Dikes/dikedata/DikeMountain_Peaks_3857WKT.csv')
#dikeset=pd.read_csv('/home/akh/myprojects/Linking-and-Clustering-Dikes/dikedata/Peaks_3857.csv')
dikeset=WKTtoArray(dikeset)
dikeset=giveID(dikeset)
theta, rho, xc, yc= HT(dikeset)
dikeset['rho']=rho
dikeset['theta']=theta

dikeset.to_csv('/home/akh/myprojects/Linking-and-Clustering-Dikes/dikedata/DikeMountain_Peaks_3857.csv',index=False)

trange=2 
rrange=5000 

stdT=np.std(theta)
stdRho=np.std(rho)
d2=trange/stdT

clustering=HT_AGG(dikeset,d2)
dikeset['Labels']=clustering.labels_
dikeset['p']=rho
dikeset['theta']=theta
Splines,IC=examineClusters(dikeset)

Splines=transformXstart(Splines)

fig, ax= DotsHT(dikeset,Splines)

Splines.to_csv('/home/akh/myprojects/Linking-and-Clustering-Dikes/dikedata/SpanishSilverLinked0621.csv')