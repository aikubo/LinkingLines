#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 23 00:16:47 2021

@author: akh
"""

from pyproj import Proj
import pandas as pd 

lines=pd.read_csv('/home/akh/myprojects/Linking-and-Clustering-Dikes/CRB_linked.csv')

myProj = Proj("+proj=utm +zone=11, +ellps=WGS84 +datum=WGS84 +units=m")
lon1, lat1 = myProj(lines['Xstart'], lines['Ystart'], inverse = True)
lon2, lat2 = myProj(lines['Xend'], lines['Yend'], inverse = True)

lines['Lon1']=lon1
lines['Lon2']=lon2
lines['Lat1']=lat1
lines['Lat2']=lat2

linesDict={name: col.values for name, col in lines.items()}

from scipy.io import savemat

savemat('CJDS_linked.mat', linesDict)