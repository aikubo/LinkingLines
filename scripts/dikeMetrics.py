#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 24 11:52:12 2021

@author: akh

Calculates dike metrics
"""

import numpy as np 
import pandas as pd 

def addSlopeIntercept(df):
    df['Slope']=(df['Ystart']-df['Yend'])/(df['Xstart']-df['Xend'])
    
    df['Int']=df['Ystart'].values-df['Xstart'].values*df['Slope'].values
    
    return df

def dikeDensity(df,gridSize, searchRadius ):
    '''
    calculates the density of dikes over the grid size
    after the arcgis LineDensity function

    Parameters
    ----------
    df : pandas data frame
        
    gridSize: float
        DESCRIPTION.

    searchRadius: 
        search radiusl around the grid
        
    Returns
    -------
    density: n x m ndarray
        

    '''
    if "Slope" not in df.columns:
        df=addSlopeIntercept(df)
        
    Xextent=[min(df['Xstart'], df['Xend']), max(df['Xstart'], df['Xend'])]
    Yextent=[min(df['Ystart'], df['Yend']), max(df['Ystart'], df['Yend'])]
    
    gridx, gridy=np.meshgrid( np.arange(Xextent[0], Xextent[1], gridSize),
                             np.arange(Yextent[0], Yextent[1], gridSize))
    
    