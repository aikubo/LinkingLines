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
from datetime import datetime
dikeset=pd.read_csv('/home/akh/myprojects/Linking-and-Clustering-Dikes/dikedata/crb/CJDS_FebStraightened.csv')
# already straightened dikes and did preprecessing

# '/home/akh/myprojects/Linking-and-Clustering-Dikes/dikedata/crb/CJDS_sensitivity'
def RhoSensitivityRun(dikeset, rhos, evalPath, rotate=False, saveTemp=True, overideSave=False):
    """
    Run sensitivity analysis based on changing Rho threshold in clustering algo
    Outputs evaluation of clusters

    Parameters
    ----------
    dikeset : pandas dataframe
        DESCRIPTION.
    rhos : list
        list of rho values as threshold.
    evalPath : string
        path to save evaluations.
    rotate : bool, optional
        Use the rotateAndCluster method. The default is False.
    saveTemp: bool, default True
    overideSave: bool, default False

    Returns
    -------
    evals : pandas dataframe
        evaluation of clusters including 
        max length, cluster sizes, cluster max size.

    """
    evals=pd.DataFrame()
    
    dtheta=2

    for drho in rhos:
        print("Running ", dtheta, drho)
        
        savePath=evalPath+"_Rho"+str(drho)
        #Run the clustering algorithm
        if ~overideSave:
            isExist = os.path.exists(savePath)
            if isExist:
                temp=pd.read_csv(savePath)
                evaluation=evaluationOnClusters(temp)
                
            else:
                if rotate: 
                    temp=RotateAndCluster(dikeset, dtheta, drho)
                else:
                    temp, M=HT_AGG_custom(dikeset, dtheta, drho, linkage='complete')
                    
                if saveTemp:
                    temp.to_csv(savePath)
                
                #Run examinecluster
                temp,evaluation=examineClusters(temp)
                evaluation['Theta_Threshold']=dtheta
                evaluation['Rho_Threshold']=drho
                evaluation['Linkage']='Complete'
    
        evals=evals.append(evaluation, ignore_index=True)

    now = datetime.now() 
    date = now.strftime("%d %b, %Y")
    evals=evals.assign(Date_Changed=date)
    evals.to_csv(evalPath)
    
    return evals

def LinkageSensitivityRun(dikeset, evalPath, linkageTypes=['complete', 'single', 'average'], drho=1000, dtheta=2, saveTemp=True, overideSave=False):
    evals=pd.DataFrame()

    for l in linkageTypes:
        print("Running ", l)
        
        savePath=evalPath+"_Linkage"+l
        
        if ~overideSave:
            isExist = os.path.exists(savePath)
            if isExist:
                continue
        #Run the clustering algorithm
        
        temp, M=HT_AGG_custom(dikeset, dtheta, drho, linkage=l)
    
        #Run examinecluster
        temp,evaluation=examineClusters(temp)
        evaluation['Theta_Threshold']=dtheta
        evaluation['Rho_Threshold']=drho
        evaluation['Linkage']=l
        
    
        evals=evals.append(evaluation, ignore_index=True)
        
        savePath=evalPath+"_Linkage"+l
        if saveTemp:
            temp.to_csv(savePath)
            
    now = datetime.now() 
    date = now.strftime("%d %b, %Y")
    evals=evals.assign(Date_Changed=date)
    
    evals.to_csv(evalPath)
    
    return evals

def runAllSensitivity(dikeset, path):
    evalPath=path+"RhoSensitivity"
    rhos=[100, 200,400,1000,5000]
    evals1= RhoSensitivityRun(dikeset, rhos, evalPath)
    
    evalPath=path+"LinkageSensitivity"
    evals2=LinkageSensitivityRun(dikeset, evalPath)
    
    evalPath=path+"RhoSensitivityRotated"
    evals1= RhoSensitivityRun(dikeset, rhos, evalPath, rotate=True)
    