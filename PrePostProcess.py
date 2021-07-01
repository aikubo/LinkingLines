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

import numpy as np 
def writeToQGIS(df,name):
    front="LINESTRING("
    linestring=[]
    for i in range(len(df)):
        line=front+str(df['Xstart'].iloc[i])+" "+str(df['Ystart'].iloc[i])+","+str(df['Xend'].iloc[i])+" "+str(df['Yend'].iloc[i])+")"
        linestring.append(line)
    
    df['Linestring']=linestring 
    df['Linestring']=df['Linestring'].astype(str)
    df.to_csv(name)
    
    return df


def WKTtoArray(df):
    import re
    xstart=[]
    ystart=[]
    
    xend=[]
    yend=[]
    drop=[]
    for i in range(len(df)):
        temp=df["WKT"].iloc[i]
        temp=re.split(r'[,\s]+', temp[18:-2])
        if len(temp) < 4 :
            drop.append(i)
            continue

        xstart.append(float(temp[0]))
        ystart.append(float(temp[1]))
        xend.append(float(temp[-2]))
        yend.append(float(temp[-1]))
    length=np.sqrt((np.array(xstart)-np.array(xend))**2+(np.array(ystart)-np.array(yend))**2)
    df=df.drop(drop)
    df['Xstart']=xstart
    df['Ystart']=ystart
    df['Xend']=xend
    df['Yend']=yend
    df['seg_length']=length
    
    return df

def giveID(df):
    col=[x.upper() for x in df.columns]
    
    if "id" not in col:
        newID=np.arange(0,len(df))
    else: 
        idCol = [i for i, s in enumerate(col) if 'ID' in s]
        newID=df[df.columns[idCol][0]]
        if np.unique(newID) < len(df):
            newID=np.arange(0,len(df))
    df['ID']=newID
    
    return df

