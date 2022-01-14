#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul  1 11:04:53 2021

@author: akh

Contains various preprocessing data and post processing 

    writeToWQGIS makes the dataframe into a WKT (well known text) string and writes as CSV (comma seperated values)
    WKTtoArray processes from a WKT CSV to a pandas Dataframe with columns Xstart,Ystart,Xend,Yend,seg_length
    giveID gives a numeric ID to data
"""
import pandas as pd
import numpy as np 
def writeToQGIS(df,name, myProj=None):
    """
    
    Writes a dataframe to a CSV file in the format of a QGIS layer.
    Uses well known text (WKT) to write the data as readable line vectors in QGIS.
    Uses the WKTtoArray function to process the dataframe.
    
    Input:
        df: a pandas dataframe with columns Xstart,Ystart,Xend,Yend,seg_length
        name: the name of the file to be written
        myProj: the projection of the dataframe. If none is given, the projection is set to WGS84
        Output:
            A CSV file with the dataframe in the format of a QGIS layer
            """
    
    #if myProj not None: 
        
    front="LINESTRING("
    linestring=[]
    for i in range(len(df)):
        line=front+str(df['Xstart'].iloc[i])+" "+str(df['Ystart'].iloc[i])+","+str(df['Xend'].iloc[i])+" "+str(df['Yend'].iloc[i])+")"
        linestring.append(line)
    
    df['Linestring']=linestring 
    df['Linestring']=df['Linestring'].astype(str)
    df.to_csv(name)
    
    return df

def writeToQGISLong(df,name, myProj=None):
    
    """
    Writes a dataframe to a CSV file in the format of a QGIS layer.
    Uses well known text (WKT) to write the data as readable line vectors in QGIS.

    Input:
        df: a pandas dataframe with columns Xstart,Ystart,Xend,Yend,seg_length

        name: the name of the file to be written

        myProj: the projection of the dataframe. If none is given, the projection is set to WGS84

    Output:
        A CSV file with the dataframe in the format of a QGIS layer
    """
        
    front="LINESTRING("
    linestring=[]
    for i in range(len(df)):
        line=front+str(df['XstartL'].iloc[i])+" "+str(df['YstartL'].iloc[i])+","+str(df['XendL'].iloc[i])+" "+str(df['YendL'].iloc[i])+")"
        linestring.append(line)
    
    df['Linestring']=linestring 
    df['Linestring']=df['Linestring'].astype(str)
    df.to_csv(name)
    
    return df


def WKTtoArray(df):

    '''
    Processes a dataframe with columns Xstart,Ystart,Xend,Yend,seg_length to a pandas dataframe with columns Xstart,Ystart,Xend,Yend,seg_length

    Input:
        df: a pandas dataframe with a Linestring column containing WKT strings
    Output:
        df: a pandas dataframe with columns Xstart,Ystart,Xend,Yend,seg_length
    '''
    import re
    xstart=[]
    ystart=[]
    
    xend=[]
    yend=[]
    drop=[]
    for i in range(len(df)):
        temp=df["WKT"].iloc[i]
        temp=re.split(r'[,\s]+', temp[20:-2])
        if len(temp) < 4 :
            drop.append(i)
            continue
        
       
        if len(temp)%3 == 0 and '0' in temp: 
            xstart.append(float(temp[0]))
            ystart.append(float(temp[1]))
            xend.append(float(temp[-3]))
            yend.append(float(temp[-2]))
        else: 
            xstart.append(float(temp[0]))
            ystart.append(float(temp[1]))
            xend.append(float(temp[-2]))
            yend.append(float(temp[-1]))
        #print(temp)
            
            
    length=np.sqrt((np.array(xstart)-np.array(xend))**2+(np.array(ystart)-np.array(yend))**2)
    
    if len(drop) >= 0:
        df=df.drop(drop)

    print(len(df), len(xstart))
    df['Xstart']=xstart
    df['Ystart']=ystart
    df['Xend']=xend
    df['Yend']=yend
    df['seg_length']=length
    
    return df

def giveID(df):
    """
    
    Gives a numeric ID to a dataframe."""
    col=[x.upper() for x in df.columns]
    
    if "ID" not in col:
        newID=np.arange(0,len(df))
    else: 
        idCol = [i for i, s in enumerate(col) if 'ID' in s]
        newID=df[df.columns[idCol][0]]
        if len(np.unique(newID)) < len(df):
            newID=np.arange(0,len(df))
    df['ID']=newID
    
    return df

def giveHashID(df):

    """
    Gives a hash ID to a dataframe based on line endpoints.

    Input:
        df: a pandas dataframe with columns Xstart,Ystart,Xend,Yend,seg_length
    Output:
        df: a pandas dataframe with columns HashID
    """

    ids=[]
    for i in range(len(df)):
        h=pd.util.hash_pandas_object(df.iloc[i])
        ids.append(hash( (h['Xstart'], h['Ystart'], h['Xend'], h['Yend'])))
    df['HashID']=ids
    
    return df

def segLength(df):
    length=np.sqrt( (df['Xstart']-df['Xend'])**2 + (df['Ystart']-df['Yend'])**2)
    df['seg_length']=length
    return df

    
def preprocess(df):
    col=[x.upper() for x in df.columns]
    df=giveID(df)
    if "SEG_LENGTH" not in col:
        df=segLength(df)
        
    return df

def whichForm(lines):
    '''
    Returns the form of the dataframe column names 

    Input:
        lines: a dataframe with columns containing theta, rho values 
    Output: 
        form: a string with the form of the dataframe column names 
    
    '''


    col=lines.columns 
    
    post=['Theta', 'AvgTheta', 'theta']
    posr=['Rho', 'AvgRho', 'rho']

    for p in post: 
        if p in col:
            t=p 
    
    for p in posr: 
        if p in col:
            r=p
            
    
    return t,r
    