#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 15 12:28:04 2022

@author: akh
"""
from examineMod import * 
from PrePostProcess import * 
from plotmod import * 

import pandas as pd
import numpy as np 
from htMOD import rotateData2

df=pd.read_csv('/home/akh/myprojects/Linking-and-Clustering-Dikes/dikedata/test_overlap.csv')
df=DikesetReProcess(df)
fig,ax=plt.subplots()

plotlines(df, 'k', ax, ColorBy='label')

# for i in df['label'].unique():
#     lines=df[ df['label']==i]
#     o=overlapSegments(lines)
    
#     if ~(np.isclose(o, lines['trueoverlap'].iloc[0])):
        
#         print('Test', str(i),' Failed')
#         print("Overlap calcuation error!")
#         print(o, lines['trueoverlap'].iloc[0])
#         print('                    ')
#     else: 
#         print('Test', str(i), 'passed')
#         print('                    ')
        
def RotateOverlapTest(lines):
    theta=np.mean(lines['theta'].values)
    
    dfRotated=rotateData2(lines, (90-theta))
    dfRotated=transformXstart(dfRotated)
    fig,ax=plt.subplots()

    plotlines(lines, 'k', ax)
    plotlines(dfRotated, 'r', ax)
    
    Xstart=dfRotated['Xstart'].to_numpy()
    Ystart=dfRotated['Ystart'].to_numpy()
    Xend=dfRotated['Xend'].to_numpy()
    Yend=dfRotated['Yend'].to_numpy()
    step=1
    totalL=np.sum( np.sqrt( (Xstart-Xend)**2 + (Ystart-Yend)**2))
    xs=[np.arange(min(x,y), max(x,y), step) for x,y in zip(np.floor(Xstart),np.ceil(Xend))]
    
    l=np.max([len(xi) for xi in xs])
    xs_sameLength=[np.append(xi,[np.nan]*(l-len(xi))) for xi in xs]

    arr=np.vstack(xs_sameLength) #better
    u,xcounts=np.unique(arr[~np.isnan(arr)], return_counts=True)
    
    overlapx=np.float64(np.sum(xcounts[xcounts>1]-1)*step)
    
    print(overlapx)
    print(totalL)
    overlap=overlapx/(totalL-overlapx+0.00001)
    if overlapx > totalL:
        overlap=overlapx/totalL
    
    return overlap

print('                    ')

for i in df['label'].unique():
    lines=df[ df['label']==i]
    o=RotateOverlapTest(lines)

    if ~(np.isclose(o, lines['trueoverlap'].iloc[0])):
        print('Test', str(i),' Failed')
        print("Overlap calcuation error!")
        print(o, lines['trueoverlap'].iloc[0])
        print('                    ')
        


    else: 
        print('Test', str(i), 'passed')
        print('                    ')
        
        

