#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 15 12:28:04 2022

@author: akh
"""
from examineMod import * 
from PrePostProcess import * 
from plotmod import * 
from dilationCalculation import dilation
import pandas as pd
import numpy as np 
from htMOD import rotateData2
from synthetic import makeLinear2


def SimpleTestFunc(Value, Expected, TestName):
    if type(Value) is list and type(Expected) is list and type(TestName) is list: 
        passed=0
        for Value1, Expected1, TestName1 in zip(Value, Expected, TestName):
            
            if np.isclose(Expected1, Value1):
                print("Pass", TestName1)
                passed=passed+1
            else:
                print("Fail", TestName1)
                print("Expected value", Expected1, "Returned", Value1)
        return passed
    else:
        if np.isclose(Expected, Value):
            print("Pass", TestName)
        else:
            print("Fail", TestName)
            print("Expected value", Expected, "Returned", Value)
        return np.isclose(Expected, Value)*1
    
    
def TestDilation():
    #df=makeLinear2(1000, 1, 0, 5000, 500, ndikes=5)
    df2=makeLinear2(1000, 90, 0, 5000, 500, ndikes=5)
    df3=makeLinear2(1000, 45, 0, 5000, 500, ndikes=5)
    
    
    #DotsLines(df, ColorBy='Labels')
    fig, ax=DotsLines(df2, ColorBy='Labels')
    fig, ax=DotsLines(df3, ColorBy='Labels')

    EWDilation, NSDilation, _, _=dilationPlot(df2, binWidth=1)
    EWDilation2, NSDilation2, _, _=dilationPlot(df3, binWidth=1)

    print("Testing 90 Degrees")
    passed=0
    if np.isclose(np.max(EWDilation), 0):
        print("Pass EW 1")
        passed=passed+1
    else:
        print("Fail EW 1")
        print("Expected value", 0, "Returned", np.max(EWDilation))

    if np.isclose(np.max(NSDilation), len(df2)):
        print("Pass NS 1")
        passed=passed+1
    else:
        print("Fail NS 1")
        print("Expected value", len(df2), "Returned", np.max(NSDilation))

    print("Testing 45 Degrees")
    EW45=np.cos(np.deg2rad(45))*2 #acounts for where tips interact
    if np.isclose(np.max(EWDilation2), EW45):
        print("Pass EW 2")
        passed=passed+1
    else:
        print("Fail EW 2")
        print("Expected value", EW45, "Returned", np.max(EWDilation2))

    NS45=np.sin(np.deg2rad(45))*5
    if np.isclose(np.max(NSDilation2), NS45):
        print("Pass NS 2")
        passed=passed+1
    else:
        print("Fail NS 2")
        print("Expected value", NS45, "Returned", np.max(NSDilation2))
        
        
        
    print('Testing non-1 averageWidth')
    EWDilation, NSDilation, _, _=dilationPlot(df2, binWidth=1, averageWidth=10)
    Value=[np.max(NSDilation), np.max(EWDilation)]
    Expected=[10*len(df2), 0]
    TestName=['NS 3', 'EW 3']
    p=SimpleTestFunc(Value, Expected, TestName)
        
    passed=passed+p
    
    EWDilation, NSDilation, _, _=dilationPlot(df3, binWidth=1, averageWidth=10)
    Value=[np.max(NSDilation), np.max(EWDilation)]
    Expected=[NS45*10, EW45*10]
    TestName=['NS 4', 'EW 4']
    p=SimpleTestFunc(Value, Expected, TestName)
        
    passed=passed+p
    
    
    print('Testing Expanded Method')
    EWDilation, NSDilation, _, _=dilationPlot(df2, binWidth=1, averageWidth=10, method='Expanded')
    Value=[np.max(NSDilation), np.max(EWDilation)]
    Expected=[10*len(df2), 0]
    TestName=['NS 5', 'EW 5']
    p=SimpleTestFunc(Value, Expected, TestName)
        
    passed=passed+p
    
    EWDilation, NSDilation, _, _=dilationPlot(df3, binWidth=1, averageWidth=10, method='Expanded')
    Value=[np.max(NSDilation), np.max(EWDilation)]
    Expected=[NS45*10, EW45*10]
    TestName=['NS 6', 'EW 6']
    p=SimpleTestFunc(Value, Expected, TestName)
    passed=passed+p
    total=12
    print( passed, "of total tests", total)

''' Test Overlap Calculation'''
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

def AllOverlapTests():
    print("Entering OverlapTests")
    
    df=pd.read_csv('/home/akh/myprojects/Linking-and-Clustering-Dikes/dikedata/test_overlap.csv')
    df=DikesetReProcess(df)
    fig,ax=plt.subplots()
    
    plotlines(df, 'k', ax, ColorBy='label')
    
    for i in df['label'].unique():
        lines=df[ df['label']==i]
        o=overlapSegments(lines)
        
        if ~(np.isclose(o, lines['trueoverlap'].iloc[0])):
            
            print('Test', str(i),' Failed')
            print("Overlap calcuation error!")
            print(o, lines['trueoverlap'].iloc[0])
            print('                    ')
        else: 
            print('Test', str(i), 'passed')
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
            
            


