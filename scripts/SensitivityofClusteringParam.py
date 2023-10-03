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
import seaborn as sns
import matplotlib.pyplot as plt

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
                temp, M=HT_AGG_custom(dikeset, dtheta, drho, linkage='complete', rotate=True, metric='Euclidean')
                    

                
                #Run examinecluster
                lines,evaluation=examineClusters(temp)
                lines=lines.assign(Rho_Threshold=drho, Theta_Threshold=dtheta, Linkage='Complete')
                if saveTemp:
                    lines.to_csv(savePath)
                evaluation['Theta_Threshold']=dtheta
                evaluation['Rho_Threshold']=drho
                evaluation['Linkage']='Complete'
    
        evals=evals.append(evaluation, ignore_index=True)

    now = datetime.now() 
    date = now.strftime("%d %b, %Y")
    evals=evals.assign(Date_Changed=date)
    evals.to_csv(evalPath)
    
    return evals

def thetaSensitivityRun(dikeset, thetas, evalPath, overideSave=False, saveTemp=True):
    """
    Run sensitivity analysis based on changing theta threshold in clustering algo
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
    
    drho=400

    for dtheta in thetas:
        print("Running ", dtheta, drho)
        
        savePath=evalPath+"_Theta"+str(dtheta)
        #Run the clustering algorithm
        if ~overideSave:
            isExist = os.path.exists(savePath)
            if isExist:
                temp=pd.read_csv(savePath)
                evaluation=evaluationOnClusters(temp)
                
            else:
                temp, M=HT_AGG_custom(dikeset, dtheta, drho, linkage='complete', rotate=True, metric='Euclidean')
                    

                
                #Run examinecluster
                now = datetime.now() 
                date = now.strftime("%d %b, %Y")
                temp,evaluation=examineClusters(temp)
                
                temp=temp.assign(Rho_Threshold=drho, Theta_Threshold=dtheta, Linkage='Complete', dateChanged=date)
                if saveTemp or overideSave:
                    temp.to_csv(savePath)
                evaluation['Theta_Threshold']=dtheta
                evaluation['Rho_Threshold']=drho
                evaluation['Linkage']='Complete'



def LinkageSensitivityRun(dikeset, evalPath, linkageTypes=['complete', 'single', 'average'], drho=400, dtheta=2, saveTemp=True, overideSave=False):
    evals=pd.DataFrame()
    now = datetime.now() 
    date = now.strftime("%d %b, %Y")
    for l in linkageTypes:
        print("Running ", l)
        
        savePath=evalPath+"_Linkage"+l
        
        if ~overideSave:
            isExist = os.path.exists(savePath)
            if isExist:
                continue
        #Run the clustering algorithm
        
        temp, M=HT_AGG_custom(dikeset, dtheta, drho, linkage=l,  rotate=True, metric='Euclidean')
    
        #Run examinecluster
        temp,evaluation=examineClusters(temp)
        evaluation['Theta_Threshold']=dtheta
        evaluation['Rho_Threshold']=drho
        evaluation['Linkage']=l
        
    
        evals=evals.append(evaluation, ignore_index=True)
        temp=temp.assign(Rho_Threshold=drho, Theta_Threshold=dtheta, Linkage=l, dateChanged=date)
        savePath=evalPath+"_Linkage"+l
        if saveTemp:
            temp.to_csv(savePath)
            

    evals=evals.assign(Date_Changed=date)
    
    evals.to_csv(evalPath)
    
    return evals


    
def sensitivityFigures(path):
    SMALL_SIZE = 8
    MEDIUM_SIZE = 8
    BIGGER_SIZE = 8
  
    plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
    plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
    plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
    plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
    plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
    plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
    plt.rc('figure', titlesize=MEDIUM_SIZE)  # fontsize of the figure title
    plt.rcParams['legend.title_fontsize'] = SMALL_SIZE
    
    #rho 
    evalPath=path+"RhoSensitivity"
    rhos=[200,400,1000,5000]
    rho_files=[evalPath+"_Rho"+str(drho) for drho in rhos]
    dfAll = pd.concat((pd.read_csv(f) for f in rho_files), ignore_index=True)
    dfAll['Dike Cluster Length (km)']=dfAll['R_Length']/1000
    dfAll['Dike Cluster Width (m)']=dfAll['R_Width']+2
    dikeset['Segment Length (km)']=dikeset['seg_length']/1000
    dfAll=dfAll.assign(Status=[ 'Filtered Cluster' if i==1 else 'Cluster' for i in dfAll['TrustFilter']])
    
    fig,ax=plt.subplots(3,3)
    fig.set_size_inches(190/25.4, 190/25.4)
    
    sns.boxplot(data=dfAll, x='Rho_Threshold', y='Size', ax=ax[0,0], hue='Status')

    g=sns.boxplot(data=dfAll, x='Rho_Threshold', y='Dike Cluster Length (km)', ax=ax[1,0], hue='Status', )
    g.legend_.remove()
    g=sns.boxplot(data=dfAll, x='Rho_Threshold', y='Dike Cluster Width (m)', ax=ax[2,0], hue='Status')
    g.legend_.remove()

    evalPath=path+"LinkageSensitivity"
    linkages=['complete', 'single', 'average']
    files=[evalPath+"_Linkage"+l for l in linkages]
    dfAll = pd.concat((pd.read_csv(f) for f in files), ignore_index=True)
    dfAll['Dike Cluster Length (km)']=dfAll['R_Length']/1000
    dfAll['Dike Cluster Width (m)']=dfAll['R_Width']+2
    dikeset['Segment Length (km)']=dikeset['seg_length']/1000
    dfAll=dfAll.assign(Status=[ 'Filtered Cluster' if i==1 else 'Cluster' for i in dfAll['TrustFilter']])
    # fig,ax=plt.subplots(3,1)
    # fig.set_size_inches(75/25.4, 225/25.4)
    g=sns.boxplot(data=dfAll, x='Linkage', y='Size', ax=ax[0,1], hue='Status')
    g.legend_.remove()
    g=sns.boxplot(data=dfAll, x='Linkage', y='Dike Cluster Length (km)', ax=ax[1,1], hue='Status')
    g.legend_.remove()
    g=sns.boxplot(data=dfAll, x='Linkage', y='Dike Cluster Width (m)', ax=ax[2,1], hue='Status')
    g.legend_.remove()

    
    evalPath=path+"ThetaSensitivity"
    thetas=[1,2,4,10]
    files=[evalPath+"_Theta"+str(l) for l in thetas]
    dfAll = pd.concat((pd.read_csv(f) for f in files), ignore_index=True)
    dfAll['Dike Cluster Length (km)']=dfAll['R_Length']/1000
    dfAll['Dike Cluster Width (m)']=dfAll['R_Width']+2
    dikeset['Segment Length (km)']=dikeset['seg_length']/1000
    dfAll=dfAll.assign(Status=[ 'Filtered Cluster' if i==1 else 'Cluster' for i in dfAll['TrustFilter']])
    # fig,ax=plt.subplots(3,1)
    # fig.set_size_inches(75/25.4, 225/25.4)
    g=sns.boxplot(data=dfAll, x='Theta_Threshold', y='Size', ax=ax[0,2], hue='Status')
    g.legend_.remove()
    g=sns.boxplot(data=dfAll, x='Theta_Threshold', y='Dike Cluster Length (km)', ax=ax[1,2], hue='Status')
    g.legend_.remove()

    g=sns.boxplot(data=dfAll, x='Theta_Threshold', y='Dike Cluster Width (m)', ax=ax[2,2], hue='Status')

    fig.savefig(path+"SensitivityFigure.pdf", dpi=600)
    g=sns.boxplot(data=dfAll, x='Theta_Threshold', y='Dike Cluster Width (m)', ax=ax[2,2], hue='Status')
    g.legend_.remove()
    
    for a in [ax[1,0], ax[2,0], ax[1,1], ax[2,1], ax[2,2], ax[1,2], ax[0,1]]:
        a.set_yscale('log')
        
    for a in ax[:,0]:
        a.set_xlabel('Rho Threshold (m)')
    for a in ax[:,2]:
        a.set_xlabel('Theta Threshold ($^\circ$)')
        
    
    labelSubplots([a for a in ax.flatten()])
    plt.tight_layout()
    
    fig.savefig(path+'SenseFigurelog.pdf', dpi=600)
    
def runAllSensitivity(dikeset, path):
    evalPath=path+"RhoSensitivity"
    rhos=[200,400,1000,5000]
    evals1= RhoSensitivityRun(dikeset, rhos, evalPath)
    
    evalPath=path+"LinkageSensitivity"
    evals2=LinkageSensitivityRun(dikeset, evalPath)
    
    evalPath=path+"ThetaSensitivity"
    thetas=[1,2,4,10]
    thetaSensitivityRun(dikeset, thetas, evalPath)
    
    sensitivityFigures(path)
    
dikeset=pd.read_csv('/home/akh/myprojects/Linking-and-Clustering-Dikes/dikedata/crb/CJDS_FebStraightened.csv')
runAllSensitivity(dikeset, '/home/akh/myprojects/Linking-and-Clustering-Dikes/dikedata/crb/CJDS_sensitivity/')