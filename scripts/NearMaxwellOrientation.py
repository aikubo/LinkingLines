#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 29 14:05:34 2022

@author: akh
"""
from PrePostProcess import * 
import pandas as pd 
import numpy as np 
from htMOD import AKH_HT as HT 
from htMOD import HT_center, MidtoPerpDistance
from plotmod import *
import matplotlib.pyplot as plt 
from clusterMod import HT_AGG_custom
import scipy.cluster.hierarchy as sch
from examineMod import examineClusters, checkoutCluster, TopHTSection, CheckoutBy
import os
import labellines
from scipy import stats
import seaborn as sns
from scipy.stats import lognorm
import time
import dateutil.parser
from datetime import datetime 
from datetime import timedelta
from fitRadialCenters import RadialFit, CenterFunc, NearCenters


lines=pd.read_csv('/home/akh/myprojects/Linking-and-Clustering-Dikes/dikedata/crb/CJDS_Complete_2_433.csv')
dikeset=pd.read_csv("/home/akh/myprojects/Linking-and-Clustering-Dikes/dikedata/crb/CJDS_FebStraightened.csv")

limits=[ (-40,-0), (-40,-20)]

for i in limits: 
    linesTemp=lines[ (lines['AvgTheta']>i[0]) & (lines['AvgTheta']<i[1])]
    dfTemp=dikeset[  (dikeset['theta']>i[0]) & (dikeset['theta']<i[1])]
    
    print("Number of Segments for Orientation", i )
    print(len(dfTemp))
    
    print("Number of Clusters for Orientation", i )
    print(len(linesTemp[linesTemp['Linked']==1]))
    
    print("Number of TrustedCluster for Orientation", i )
    print(len(linesTemp[linesTemp['TrustFilter']==1]))
    