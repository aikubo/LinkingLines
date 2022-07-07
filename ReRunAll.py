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

d=[#'dikedata/deccandata/Central_preprocesed.csv',
   #'dikedata/deccandata/NarmadaTapi_preprocesed.csv',
   #'dikedata/deccandata/Saurashtra_preprocesed.csv',
   #'dikedata/crb/CJDS_FebStraightened.csv',
   'dikedata/spanish peaks/SpanishPeaks_3857_preprocessed.csv'
   ]

l=[#'dikedata/deccandata/Central_Complete3_', 
   # 'dikedata/deccandata/NarmadaTapi_Complete_3_',
   # 'dikedata/deccandata/Saurashtra_Complete_3_',
   # 'dikedata/crb/CJDS_Complete_3',
    'dikedata/spanish peaks/SpanishPeaks_Complete_3_']

for i,j in zip(d,l): 
    dikeset=pd.read_csv(i)

    dikeset=DikesetReProcess(dikeset)
    dtheta=3 
    drho=np.floor(dikeset['seg_length'].mean())
    dikeset, Z=HT_AGG_custom(dikeset, dtheta, drho, linkage='complete')
    lines,evaluation=examineClusters(dikeset)
    now = datetime.now() 
    date = now.strftime("%d %b, %Y")

    dikeset=dikeset.assign(Date_Changed=date)
    lines=lines.assign(Date_Changed=date)
    
    name=j+str(int(drho))+"_.csv"
    npname=j+str(int(drho))+"LINKAGE"
    np.save(npname, Z)
    lines.to_csv(name)
    dikeset.to_csv(i)
    
    allFigures(i, name)
