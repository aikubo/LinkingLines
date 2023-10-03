#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  7 16:20:24 2022

@author: akh
"""
import pandas as pd 
import matplotlib.pyplot as plt 
import numpy as np 
from htMOD import AKH_HT as HT
from htMOD import MidtoPerpDistance
from clusterMod import HT_AGG_custom as AggHT
from examineMod import examineClusters
from plotmod import *
from PrePostProcess import * 

dikeset=pd.read_csv('/home/akh/myprojects/Linking-and-Clustering-Dikes/dikedata/deccandata/deccan_maincontientaldikes_32244_processed.csv')
#dikeset=completePreProcess(dikeset)

# already straightened dikes and did preprecessing

#do the hough transform 
theta,rho,xc,yc=HT(dikeset)
dikeset['theta']=theta
dikeset['rho']=rho
dikeset=MidtoPerpDistance(dikeset, xc, yc)


#DotsLinesHT(dikeset, ColorBy="theta")
#evals=pd.read_csv('/home/akh/myprojects/Linking-and-Clustering-Dikes/CJDS_sensitivity.csv')

evals=pd.DataFrame()

rhos=[200,500,1000,5000]
#Run sensitivity analysis first
thetas=[1,2,3,5,6]
for dtheta in thetas:
    lengths=[]
    needikes=[]
    nclusters=[]
    for drho in rhos:
        print("Running ", dtheta, drho)
        #Run the clustering algorithm
        temp, clusters, M=AggHT(dikeset, dtheta, drho)
        
        #Run examinecluster
        temp,evaluation,EEdikes=examineClusters(temp, ifEE=True)
        evaluation['Theta_Threshold']=dtheta
        evaluation['Rho_Threshold']=drho
        
        evals=evals.append(evaluation, ignore_index=True)

evals=pd.read_csv('/home/akh/myprojects/Linking-and-Clustering-Dikes/deccan_sensitivity2.csv')
def plotSensitivityRhoColumns(evals, columns):
    xlim=[(150,1050), (4950,5005)]
    f1,ax11,ax12=breakXaxis(xlim, numAxes=len(columns))
    rhos=np.unique(evals['Rho_Threshold'])
    theta=np.unique(evals['Theta_Threshold'])
    markers= ["rs-", "yv-", "g*-", "bp-", "kD-" ]
    for i,j in zip(columns, range(len(columns))):
        
        for t,m in zip(theta, markers): 
            y1=evals.loc[evals['Theta_Threshold']==t][i]
            name=str(t)+"$^\circ$"
            kwargs=dict( label=name, markersize=7, alpha=0.5)
            plotBreak(xlim, rhos, y1, ax11[j], ax12[j], m, **kwargs)
            ax11[j].set_ylabel(i)
    plt.tight_layout()
    return f1, (ax11, ax12)


def plotSensitivityThetaColumns(evals, columns):
    
    f1,ax=plt.subplots(len(columns), sharex=True)
    rhos=np.unique(evals['Rho_Threshold'])
    theta=np.unique(evals['Theta_Threshold'])
    markers= ["rs-", "yv-", "g*-", "bp-" ]
    
    for i,j in zip(columns, range(len(columns))):
        
        for p,m in zip(rhos, markers): 
            y1=evals.loc[evals['Rho_Threshold']==p][i]
            name=str(p)+"m"
            kwargs=dict( label=name, markersize=7, alpha=0.5)
            ax[j].plot(theta, y1, m, **kwargs)
            ax[j].set_ylabel(i)
    
    ax[j].set_xlabel("Theta ($^\circ$) ")
    plt.tight_layout()
    
    return f1, ax