#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 28 10:28:12 2022

@author: akh
"""
import pandas as pd 
import matplotlib.pyplot as plt 
import numpy as np 
from htMOD import AKH_HT as HT
from htMOD import MidtoPerpDistance
from clusterMod import HT_AGG_custom as AggHT
from clusterMod import HT_AGG_custom
from examineMod import examineClusters
from plotmod import *
import scipy.cluster.hierarchy as sch
from PrePostProcess import * 
import pandas as pd 
import numpy as np 
from htMOD import AKH_HT as HT 
from htMOD import HT_center, MidtoPerpDistance
from plotmod import *
import matplotlib.pyplot as plt 
from clusterMod import HT_AGG_custom
import scipy.cluster.hierarchy as sch
from examineMod import examineClusters, checkoutCluster, TopHTSection
import os
import labellines
from ResultsPlots import *
import seaborn as sns
from datetime import datetime

d=['dikedata/deccandata/Central_preprocesed.csv',
   'dikedata/deccandata/NarmadaTapi_preprocesed.csv',
   'dikedata/deccandata/Saurashtra_preprocesed.csv',
   'dikedata/crb/CJDS_FebStraightened.csv',
   'dikedata/crb/Steens_dikes_FebStraightened.csv',
   'dikedata/crb/Monument_Dikes_FebStraightened.csv',
   'dikedata/crb/Ice_Harbor_FebStraightened.csv',
   'dikedata/spanish peaks/SpanishPeaks_3857_preprocessed.csv'
   ]

l=['dikedata/deccandata/Central_Complete_euc_2_', 
   'dikedata/deccandata/NarmadaTapi_Complete_euc_2_',
   'dikedata/deccandata/Saurashtra_Complete_euc_2_',
   'dikedata/crb/CJDS_Complete_euc_2_',
   'dikedata/crb/Steens_Complete_euc_2',
   'dikedata/crb/Monument_Complete_euc_2',
   'dikedata/crb/IceHarbor_Complete_euc_2',
   'dikedata/spanish peaks/SpanishPeaks_Complete_euc_3_']
names=['Deccan:Central',
       'Deccan:Narmada-Tapi',
       'Deccan:Saurashtra',
       'CRBG:CJDS',
       'CRBG:Steens',
       'CRBG:Monument',
       'CRBG:Ice Harbor',
       'Spanish Peaks']
def printInfo(i,j):
    print( "Number of Segments", len(i))
    
    print('Rho threshold:', j['Rho_Threshold'].unique())
    print( "Number of Clusters", np.sum(j['Linked']==1))
    print( "Number of Trusted Clusters", np.sum(j['TrustFilter']==1))
    print( 'Center', HT_center(i))
    
overide=False
examineOveride=True
width=[8,8,8,8,8,8,8,2]
allDeccanLinked=pd.DataFrame()
allCRBLinked=pd.DataFrame()
changedClustering=False
for i,j,w, n in zip(d,l,width, names): 
    print("Loading...")
    print(i)
    dikeset=pd.read_csv(i)

    
    dtheta=2 
    drho=np.floor(dikeset['seg_length'].mean())
    
    name=j+str(int(drho))+".csv"
    npname=j+str(int(drho))+"LINKAGE"+".npy"
    now = datetime.now() 
    date = now.strftime("%d %b, %Y")
    date2=now.strftime('%d_%m_%Y')
    
    if os.path.exists(name) and os.path.exists(npname) :
        print("opening")
        lines=pd.read_csv(name)

    else:
        print("Re-running clustering...")
        dikeset=DikesetReProcess(dikeset)
        dikeset, Z=HT_AGG_custom(dikeset, dtheta, drho, linkage='complete', rotate=True, metric='Euclidean')
        lines,evaluation=examineClusters(dikeset)

        
        dikeset=dikeset.assign(Date_Changed=date)
        lines=lines.assign(Date_Changed=date)
        lines=lines.assign(Rho_Threshold=drho)
        lines=lines.assign(Theta_Threshold=dtheta)

        np.save(npname, Z)
        lines.to_csv(name, index=False)
        dikeset.to_csv(i, index=False)
        changedClustering=True

    if overide:
        print("Re-running clustering...")
        dikeset=DikesetReProcess(dikeset)
        dikeset, Z=HT_AGG_custom(dikeset, dtheta, drho, linkage='complete', rotate=True, metric='Euclidean')
        lines,evaluation=examineClusters(dikeset)

        
        dikeset=dikeset.assign(Date_Changed=date)
        lines=lines.assign(Date_Changed=date)
        lines=lines.assign(Rho_Threshold=drho)
        lines=lines.assign(Theta_Threshold=dtheta)

        np.save(npname, Z)
        lines.to_csv(name, index=False)
        dikeset.to_csv(i, index=False)
        changedClustering=True
        
    if examineOveride:
        lines,evaluation=examineClusters(dikeset)
        lines=lines.assign(Date_Changed=date)
        lines=lines.assign(Rho_Threshold=drho)
        lines=lines.assign(Theta_Threshold=dtheta)
        lines.to_csv(name, index=False)
        
        changedClustering=True


    lines=lines.assign(SwarmID=n)
    printInfo(dikeset,lines)
    print("Plotting Figures")
    #allFigures(i, name, w, Large=False, overide=False)
    
    if 'deccandata' in j:
        allDeccanLinked=allDeccanLinked.append(lines)
    if 'crb' in j: 
        allCRBLinked=allCRBLinked.append(lines)
changedClustering=True


if changedClustering:      
    allCRBLinked=LinesReProcess(allCRBLinked, HTredo=True)
    
    allCRBLinked.to_csv('dikedata/crb/AllCRBLinked_euc_'+date2+'.csv', index=False)
    
    allFigures('dikedata/crb/allCRB_dikes_PreProcessed.csv', 'dikedata/crb/AllCRBLinked_euc_'+date2+'.csv', 8, Large=False)
    
    allDeccanLinked=LinesReProcess(allDeccanLinked, HTredo=True)
    allDeccanLinked.to_csv('dikedata/deccandata/AllDeccanLinked_euc_18_10_2022.csv', index=False)
        
    allFigures('/home/akh/myprojects/Linking-and-Clustering-Dikes/dikedata/deccandata/AllDeccan_PreProcessed.csv', 'dikedata/deccandata/AllDeccanLinked_euc_'+date2+'.csv',8, Large=False)

    Region=['Deccan', 'Deccan', 'Deccan', 'CRBG', 'CRBG', 'CRBG', 'CRBG', 'SpanishPeaks']
    ALLDATA=pd.DataFrame()
    ALL_LINES=pd.DataFrame()
    
    for i,j,n,r in zip(d,l,names,Region):
        df=pd.read_csv(i)
        dtheta=2 
        drho=np.floor(df['seg_length'].mean())
    
        name=j+str(int(drho))+".csv"
        lines=pd.read_csv(name)
        df=df.assign(SwarmID=n)
        lines=lines.assign(SwarmID=n)
        df=df.assign(Region=r)
        lines=lines.assign(Region=r)
        ALLDATA=ALLDATA.append(df)
        ALL_LINES=ALL_LINES.append(lines)
        
        ALLDATA.to_csv('dikedata/AllDatasets_18_10_2022.csv', index=False)
        ALL_LINES.to_csv('dikedata/AllDatasetsLinked_18_10_2022.csv', index=False)
else:
    allFigures('dikedata/crb/allCRB_dikes_PreProcessed.csv', 'dikedata/crb/AllCRBLinked_euc_18_10_2022.csv', 8, Large=False)
    allFigures('/home/akh/myprojects/Linking-and-Clustering-Dikes/dikedata/deccandata/AllDeccan_PreProcessed.csv', 'dikedata/deccandata/AllDeccanLinked_euc_18_10_2022.csv',8, Large=False)
    
        