#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 28 15:06:20 2022

@author: akh
"""

import pandas as pd 
from PrePostProcess import *

dikeset=pd.read_csv("/home/akh/myprojects/Linking-and-Clustering-Dikes/dikedata/crb/allCRB_dikes.csv")
dikeset=completePreProcess(dikeset)

swarms=['CJDS', 'Ice_Harbor', 'Monument_Dikes',
        'Steens_dikes']

for layer in swarms:
    name=layer+"_FebStraightened.csv"
    dikeset.loc[dikeset['layer']==layer].to_csv(name)