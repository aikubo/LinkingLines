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

l=['dikedata/deccandata/Central_Complete_2_', 
   'dikedata/deccandata/NarmadaTapi_Complete_2_',
   'dikedata/deccandata/Saurashtra_Complete_2_',
   'dikedata/crb/CJDS_Complete_2_',
   'dikedata/crb/Steens_Complete_2',
   'dikedata/crb/Monument_Complete_2',
   'dikedata/crb/IceHarbor_Complete_2',
   'dikedata/spanish peaks/SpanishPeaks_Complete_3_']


width=[8,8,8,8,8,8,8,2]
allDeccanLinked=pd.DataFrame()
allCRBLinked=pd.DataFrame()
for i,j,w in zip(d,l,width): 
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
    
    if os.path.exists(name) and os.path.exists(npname):
        print("opening")
        lines=pd.read_csv(name)
        
    else:
        print("Re-running clustering...")
        dikeset=DikesetReProcess(dikeset)
        dikeset, Z=HT_AGG_custom(dikeset, dtheta, drho, linkage='complete')
        lines,evaluation=examineClusters(dikeset)

        
        dikeset=dikeset.assign(Date_Changed=date)
        lines=lines.assign(Date_Changed=date)
        lines=lines.assign(Rho_Threshold=drho)
        lines=lines.assign(Theta_Threshold=dtheta)

        np.save(npname, Z)
        lines.to_csv(name, index=False)
        dikeset.to_csv(i, index=False)

        
    print("Plotting Figures")
    #allFigures(i, name, w)
    
    if 'deccandata' in j:
        allDeccanLinked=allDeccanLinked.append(lines)
    if 'crb' in j: 
        allCRBLinked=allCRBLinked.append(lines)
        
allCRBLinked.to_csv('dikedata/crb/AllCRBLinked_'+date2+'.csv', index=False)
allFigures('dikedata/crb/allCRB_dikes_PreProcessed.csv', 'dikedata/crb/AllCRBLinked_'+date2+'.csv', 8, Large=False)

allDeccanLinked.to_csv('dikedata/deccandata/AllDeccanLinked_'+date2+'.csv', index=False)
    
allFigures('/home/akh/myprojects/Linking-and-Clustering-Dikes/dikedata/deccandata/AllDeccan_PreProcessed.csv', 'dikedata/deccandata/AllDeccanLinked_'+date2+'.csv',8, Large=False)