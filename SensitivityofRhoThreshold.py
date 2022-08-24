#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 11 15:33:14 2022

@author: akh
"""

import pandas as pd 
import matplotlib.pyplot as plt 
import numpy as np 
from htMOD import AKH_HT as HT
from htMOD import MidtoPerpDistance
from clusterMod import HT_AGG_custom as AggHT
from examineMod import *
from plotmod import *

dikeset=pd.read_csv('/home/akh/myprojects/Linking-and-Clustering-Dikes/dikedata/crb/CJDS_FebStraightened.csv')
# already straightened dikes and did preprecessing

evals=pd.read_csv('/home/akh/myprojects/Linking-and-Clustering-Dikes/CJDS_sensitivity.csv')

def RhoSensitivityRun(dikeset, rhos, evalPath1, rotate=False):
    evals=pd.DataFrame()
    
    dtheta=2

    for drho in rhos:
        print("Running ", dtheta, drho)
        #Run the clustering algorithm
        
        if rotate: 
            temp=RotateAndCluster(dikeset, dtheta, drho, linkage='complete')
        else:
            temp, clusters, M=HT_AGG_custom(dikeset, dtheta, drho, linkage='complete')
    
        #Run examinecluster
        temp,evaluation=examineClusters(temp)
        evaluation['Theta_Threshold']=dtheta
        evaluation['Rho_Threshold']=drho
        evaluation['Linkage']='Complete'
    
        evals=evals.append(evaluation, ignore_index=True)

    evals.to_csv(evalPath)
    
    return evals
