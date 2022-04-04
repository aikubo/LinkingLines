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
from htMOD import MidtoPerpDistance, HT_center

def midPoint(df):
    """
    Finds the midpoint of a dataframe of line segments.
    
    Parameters
    ----------
    df: pandas.Dataframe 
        dataframe of the line segments
        must contain ["Xstart", "Ystart", "Xend", "Yend"]

    Returns
    -------
    df: pandas.Dataframe
    with new columns of ['Xmid', 'Ymid']
    """
    df['Xmid']=(df['Xstart']+df['Xend'])/2
    df['Ymid']=(df['Ystart']+df['Yend'])/2
    
    return df


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

def WKTtoArray(df, plot=True):

    '''
    Processes a dataframe with columns Xstart,Ystart,Xend,Yend,seg_length to a pandas dataframe with columns Xstart,Ystart,Xend,Yend,seg_length

    Input:
        df: a pandas dataframe with a Linestring column containing WKT strings
    Output:
        df: a pandas dataframe with columns Xstart,Ystart,Xend,Yend,seg_length
    '''
    import matplotlib.pyplot as plt
    import re
    from scipy import stats
            
    xstart=[]
    ystart=[]
    
    xend=[]
    yend=[]
    drop=[]
    if plot: 
        fig,ax=plt.subplots()
    for i in range(len(df)):
        temp=df["WKT"].iloc[i]
        temp=re.split(r'[(|)]', temp)
        t1=temp[0]
        temp=re.split(r'[,\s]+', temp[2])
        
        if "Z" in t1:
            tempx=np.array(temp[::3]).astype(float)
            tempy=np.array(temp[1::3]).astype(float)
        else:
            tempx=np.array(temp[::2]).astype(float)
            tempy=np.array(temp[1::2]).astype(float)
               
        # #print(tempx, tempy)
        # if interp: 
        #     xvals=np.linspace(min(tempx), max(tempx), 50)
        #     yvals=np.interp(xvals, tempx, tempy)
            
        #     tempx=xvals
        #     tempy=yvals

        
        slope, intercept, r_value, p_value, std_err = stats.linregress(tempx, tempy)
        #for x,y in zip(tempx, tempy):
        if any(np.isnan( [slope, intercept])):
            drop.append(i)
            continue
        
        
        
        if p_value > 0.05 and len(tempx)>3: 
            if plot:
                ax.plot(tempx, tempy)
            drop.append(i)
            continue
        
        
        x=np.array( [ np.min(tempx), np.max(tempx)])
        y=x*slope+intercept
        
        if np.sum(np.isnan( [x, y]))>0:
            drop.append(i)
            continue
        if np.sqrt( (x[1]-x[0])**2 + (y[1]-y[0])**2)<1:
            drop.append(i)
            continue
        
        
        xstart.append(x[0])
        ystart.append(y[0])
        xend.append(x[1])
        yend.append(y[1])
        if plot:
            ax.plot(x, x*slope+intercept, '*-' )
        


    length=np.sqrt((np.array(xstart)-np.array(xend))**2+(np.array(ystart)-np.array(yend))**2)
    
    if len(drop) >= 0:
         df=df.drop(drop)

    #print(len(df), len(xstart))
    df['Xstart']=xstart
    df['Ystart']=ystart
    df['Xend']=xend
    df['Yend']=yend
    df['seg_length']=length
    print( len(drop), "dropped for not being straight")
    
    
    return df
        
# def WKTtoArray(df):

#     '''
#     Processes a dataframe with columns Xstart,Ystart,Xend,Yend,seg_length to a pandas dataframe with columns Xstart,Ystart,Xend,Yend,seg_length

#     Input:
#         df: a pandas dataframe with a Linestring column containing WKT strings
#     Output:
#         df: a pandas dataframe with columns Xstart,Ystart,Xend,Yend,seg_length
#     '''
#     import re
#     xstart=[]
#     ystart=[]
    
#     xend=[]
#     yend=[]
#     drop=[]
#     for i in range(len(df)):
#         temp=df["WKT"].iloc[i]
#         temp=re.split(r'[(|)]', temp)[2]
#         temp=re.split(r'[,\s]+', temp)
        
#         if len(temp) < 4 :
#             drop.append(i)
#             continue
        
#         tempx=np.array(temp[::2]).astype(float)
#         tempy=np.array(temp[1::2]).astype(float)
        
        
#         from scipy import stats
#         slope, intercept, r_value, p_value, std_err = stats.linregress(tempx, tempy)
#         #for x,y in zip(tempx, tempy):
#         #    m= 
       
#         if len(temp)%3 == 0 and '0' in temp: 
#             xstart.append(float(temp[0]))
#             ystart.append(float(temp[1]))
#             xend.append(float(temp[-3]))
#             yend.append(float(temp[-2]))
#         else: 
#             xstart.append(float(temp[0]))
#             ystart.append(float(temp[1]))
#             xend.append(float(temp[-2]))
#             yend.append(float(temp[-1]))
#         #print(temp)
            
            
#     length=np.sqrt((np.array(xstart)-np.array(xend))**2+(np.array(ystart)-np.array(yend))**2)
    
#     if len(drop) >= 0:
#         df=df.drop(drop)

#     print(len(df), len(xstart))
#     df['Xstart']=xstart
#     df['Ystart']=ystart
#     df['Xend']=xend
#     df['Yend']=yend
#     df['seg_length']=length
    
#     return df

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

def completePreProcess(df):
    
    if 'Xstart' not in df.columns:
        df=WKTtoArray(df)

    if 'seg_length' not in df.columns:        
        df=segLength(df)
    if 'HashID' not in df.columns: 
        df=giveHashID(df)
    if 'Xmid' not in df.columns: 
        df=midPoint(df)
    if 'PerpOffsetDist' not in df.columns: 
        xc,yc=HT_center(df)
        df=MidtoPerpDistance(df, xc, yc)
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
    