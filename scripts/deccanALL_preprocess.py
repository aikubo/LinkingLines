#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 13 15:11:18 2022

@author: akh
"""
from PrePostProcess import * 
import pandas as pd 
import numpy as np 
from htMOD import AKH_HT as HT 
from htMOD import HT_center, MidtoPerpDistance
from plotmod import DotsHT
import matplotlib.pyplot as plt 
from clusterMod import HT_AGG_custom
import scipy.cluster.hierarchy as sch
from examineMod import examineClusters

data=pd.read_csv('/home/akh/myprojects/Linking-and-Clustering-Dikes/dikedata/deccandata/ALLDECCAN_041322_EPSG3857.csv')

dikeset=completePreProcess(data)

"""
1010 dropped for not being straight
"""


theta, rho, xc,  yc=HT(dikeset)
print(xc,yc)
dikeset= MidtoPerpDistance(dikeset, xc, yc)
dikeset['theta']=theta
dikeset['rho']=rho
fig,ax=plt.subplots()
DotsHT(fig, ax, dikeset, ColorBy='SwarmID')

dtheta=3 
drho=2000
        
dikeset.to_csv('/home/akh/myprojects/Linking-and-Clustering-Dikes/dikedata/deccandata/proc_ALLDECCAN_041322_EPSG3857.csv')
path='/home/akh/myprojects/Linking-and-Clustering-Dikes/dikedata/deccandata/'
file1="_preprocesed.csv"
file2="_3_2000_LINKED_complete.csv"
file3="_3_2000_CompleteLinkageMat"

""" link coastal and central"""
ncoastal=np.sum(dikeset['SwarmID']=='Coastal')
dikeset.loc[dikeset['SwarmID']=='Coastal', 'SwarmID']= ['Central']*ncoastal


for i in np.unique(dikeset['SwarmID']):
    temp=dikeset[dikeset['SwarmID']==i]
    temp=DikesetReProcess(temp)
    #print(temp['xc'].iloc[0])
    
    
    #DotsLinesHist(temp, 10000, 3, ColorBy='PerpOffsetDist')
    name1=path+i+file1
    name2=path+i+file2
    name3=path+i+file3
    
    
    temp, Z=HT_AGG_custom(temp, dtheta, drho, linkage='complete')
    np.save(name3, Z)
    temp.to_csv(name1)
    lines,evaluation=examineClusters(temp)
    lines.to_csv(name2)
        
    