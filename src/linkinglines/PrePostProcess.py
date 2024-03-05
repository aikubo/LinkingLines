# LinkingLines Package
 # Written by aikubo
 # Version: 2.1.0
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul  1 11:04:53 2021

@author: akh

Contains various preprocessing data and post processing including reading in
WKT files and exporting WKT files to use in GIS programs

    writeFile makes the dataframe into a variety of file types including .csv, .txt, .shp, .geojson, and .json
    readFile reads in a file and returns a pandas dataframe
    writeToWKT writes a dataframe to a CSV file using WKT
    writetoGeoData writes a dataframe to a shapefile, geopackage, or geojson file
    WKTtoArray processes from a WKT CSV to a pandas Dataframe with columns Xstart,Ystart,Xend,Yend,seg_length
    giveID gives a numeric ID to data
    midPoint: Finds the midpoint of a dataframe of line segments.
    giveHashID: assigns hash ID to a dataframe based on line endpoints.
    segLength: calculates segment length
    transformXstart: reorders dataframe so that Xstart is always < Xend
    DikesetReprocess: Reprocesses a dataframe containing dike line data to ensure it has essential attributes and is properly formatted.
    LinesReprocess: Reprocesses a dataframe containing line data to ensure it has essential attributes and is properly formatted.
    CompletePreprocess:     Fully preprocesses a dataframe containing line data to ensure it has essential attributes and is properly formatted.
    whichForm: Returns the form of the dataframe column names
    MaskArea: Returns dataframe masked by bounds
    getCartLimits: Computes the Cartesian limits (x and y) of a set of lines.


"""
import pandas as pd
import numpy as np
from .HT import MidtoPerpDistance, HT_center
from .HT import HoughTransform, segLength
from datetime import datetime
import re
from scipy import stats

import matplotlib.pyplot as plt

import geopandas

def readFile(path):

    """
    Reads in a file and returns a pandas dataframe

    Parameters:
        path: (string) the path to the file to be read in
    
    Returns:
        data: (pandas.DataFrame) a pandas or geopandas dataframe
    """

    # if not a valid path, return error
    if not os.path.exists(path):
        raise ValueError("Invalid path")
    # if file is not .csv, .txt, or .shp, return error
    if path.endswith('.csv') or path.endswith('.txt') or path.endswith('.shp') or path.endswith('.geojson') or path.endswith('.json'):
        raise ValueError("Invalid file type")
    
    # identify the type of file
    # read in .csv 
    if path.endswith('.csv'):
        data=pd.read_csv(path)
    elif path.endswith('.txt'):
        data=pd.read_csv(path, delimiter='\t')
    else:
        data=geopandas.read_file(path)
        data=data.to_wkt()

    return data


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


def writeToWKT(df,name, myProj=None):
    """

    Writes a dataframe to a CSV file in the format of a QGIS layer.
    Uses well known text (WKT) to write the data as readable line vectors in QGIS.
    Uses the WKTtoArray function to process the dataframe.

    Parameters:
        df: a pandas dataframe with columns Xstart,Ystart,Xend,Yend,seg_length
        name: the name of the file to be written
        myProj: the projection of the dataframe. If none is given, the projection is set to WGS84
        Returns:
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


def writetoGeoData(df, name, driver, myProj=None):
    """
    Writes a dataframe to a shapefile.

    Parameters:
        df: a pandas dataframe with columns Xstart,Ystart,Xend,Yend,seg_length
        name: the name of the file to be written
        myProj: the projection of the dataframe. If none is given, the projection is set to WGS84

    Returns:
        A shapefile with the dataframe
    """

    if myProj is None:
        myProj='WGS84'

    gdf = geopandas.GeoDataFrame(df, geometry=geopandas.points_from_xy(df.Xstart, df.Ystart))
    gdf.to_file(name, driver=driver, crs=myProj)

    return df


def WKTtoArray(df, plot=False):

    '''
    Processes a dataframe with columns Xstart,Ystart,Xend,Yend,seg_length to a pandas dataframe with columns Xstart,Ystart,Xend,Yend,seg_length

    Parameters:
        df: a pandas dataframe with a Linestring column containing WKT strings must contain ["WKT"]
    Returns:
        df: a pandas dataframe with columns Xstart,Ystart,Xend,Yend,seg_length
    '''
    if not isinstance(df, pd.DataFrame):
        raise ValueError("Input 'data' must be a pandas DataFrame.")

    if len(df) < 1:
        raise ValueError("DataFrame is empty")

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

        #print("dike #:",i)
        print(temp)
        if 'EMPTY' in temp[0]:
            drop.append(i)
            continue
        temp=re.split(r'[,\s]+', temp[2])

        if "Z" in t1:
            tempx=np.array(temp[::3]).astype(float)
            tempy=np.array(temp[1::3]).astype(float)
        else:
            tempx=np.array(temp[::2]).astype(float)
            tempy=np.array(temp[1::2]).astype(float)

        print(tempx, tempy)
        slope, intercept, r_value, p_value, std_err = stats.linregress(tempx, tempy)
        print(p_value)
        #for x,y in zip(tempx, tempy):
        if any(np.isnan( [slope, intercept])):
            drop.append(i)
            continue



        if p_value > 0.05 and len(tempx)>=3:
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

    Parameters:
        df: a pandas dataframe with columns Xstart,Ystart,Xend,Yend,seg_length
    Returns:
        df: a pandas dataframe with columns HashID
    """

    ids=[]
    for i in range(len(df)):
        h=pd.util.hash_pandas_object(df.iloc[i])
        ids.append(hash( (h['Xstart'], h['Ystart'], h['Xend'], h['Yend'])))
    df['HashID']=ids

    return df




def transformXstart(dikeset, HTredo=True):
    """
    Transforms the Xstart column of a dataframe of line segments.
    Sets the Start point as the smallest x value of the line segment.

    Parameters
    ----------
    df: pandas.Dataframe
        dataframe of the line segments
        must contain ["Xstart", "Ystart", "Xend", "Yend"]

    Returns
    -------
    df: pandas.Dataframe
        dataframe of the line segments with transformed Xstart column
    """


    dist1= dikeset['Xstart'].values
    dist2= dikeset['Xend'].values
    switchXs=(dist1>dist2)

    dikeset.loc[switchXs, ['Xstart', 'Xend']]=(dikeset.loc[switchXs, ['Xend', 'Xstart']].values)
    dikeset.loc[switchXs, ['Ystart', 'Yend']]=(dikeset.loc[switchXs, ['Yend', 'Ystart']].values)

    return dikeset

def dikesetReProcess(df, HTredo=True, xc=None, yc=None):
    """
    Reprocesses a dataframe containing dike line data to ensure it has essential attributes and is properly formatted.

    Parameters:
        df (DataFrame): The input dataframe containing dike line data.
        HTredo (bool, optional): Whether to recalculate Hough Transform attributes (default is True).
        xc (float, optional): The x-coordinate of the center point for the Hough Transform (default is None).
        yc (float, optional): The y-coordinate of the center point for the Hough Transform (default is None).

    Returns:
        DataFrame: The processed dataframe with added or updated attributes.
    """

    # Check and transform dataframe columns if necessary
    if 'Xstart' not in df.columns:
        df = WKTtoArray(df)
    if 'seg_length' not in df.columns:
        df = segLength(df)
    df = transformXstart(df)

    # Calculate Hough Transform center coordinates if not provided
    if xc is None or yc is None:
        xc, yc = HT_center(df)

    # Assign unique hash IDs to the dataframe
    df = giveHashID(df)

    # Remove duplicate entries and report any found duplicates
    l = len(df)
    df = df.drop_duplicates(subset=['HashID'])
    if l != len(df):
        print("Found", l - len(df), "duplicates")

    # Calculate midpoints if not present
    if 'Xmid' not in df.columns:
        df = midPoint(df)

    # Calculate Hough Transform attributes (theta, rho) if not present
    if 'theta' not in df.columns or 'rho' not in df.columns:
        df, _, _ = HoughTransform(df, xc=xc, yc=yc)

    # Assign or update Hough Transform center coordinates
    if 'xc' not in df.columns:
        df, xc, yc = HoughTransform(df, xc=xc, yc=yc)
        df=MidtoPerpDistance(df, xc, yc)
    elif xc is not df['xc'].iloc[0] and HTredo:
        df, xc, yc=HoughTransform(df, xc=xc, yc=yc)
        df=MidtoPerpDistance(df, xc, yc)

    elif HTredo:
        df,xc,yc=HoughTransform(df, xc=xc, yc=yc)
        df=MidtoPerpDistance(df, xc, yc)

    if 'PerpOffsetDist' not in df.columns:
        df=MidtoPerpDistance(df, xc, yc)

    now = datetime.now()
    d = now.strftime("%d %b, %Y")

    df=df.assign(Date_Changed=d)


    return df

def LinesReProcess(df, HTredo=True):
    """
    Reprocesses a dataframe containing line data to ensure it has essential attributes and is properly formatted.

    Parameters:
        df (DataFrame): The input dataframe containing line data.
        HTredo (bool, optional): Whether to recalculate Hough Transform attributes (default is True).

    Returns:
        DataFrame: The processed dataframe with added or updated attributes.
    """
    # Calculate Hough Transform center coordinates and transform 'Xstart' if necessary
    xc, yc = HT_center(df)
    df = transformXstart(df)

    # Assign unique hash IDs to the dataframe
    df = giveHashID(df)

    # Remove duplicate entries and report any found duplicates
    l = len(df)
    df = df.drop_duplicates(subset=['HashID'])
    if l != len(df):
        print("Found", l - len(df), "duplicates")

    # Calculate midpoints if not present
    if 'Xmid' not in df.columns:
        df = midPoint(df)

    # Calculate or recalculate Hough Transform attributes (theta, rho)
    if HTredo:
        df, xc, yc = HoughTransform(df)
        df = MidtoPerpDistance(df, xc, yc)
        df['AvgTheta'] = df['theta'].values
        df['AvgRho'] = df['rho'].values

    # Calculate perpendicular offset distances if not present
    if 'PerpOffsetDist' not in df.columns:
        df = MidtoPerpDistance(df, xc, yc)

    # Assign the processing date
    now = datetime.now()
    d = now.strftime("%d %b, %Y")
    df = df.assign(Date_Changed=d)

    return df

def preProcess(data):

    df = data.copy()
    """
    Fully preprocesses a dataframe containing line data to ensure it has essential attributes and is properly formatted.

    Parameters:
        df (DataFrame): The input dataframe containing line data.

    Returns:
        DataFrame: The fully processed dataframe with all required attributes.
    """
    # Convert WKT column to array format if present
    if "WKT" in df.columns:
        df = WKTtoArray(df)

    # Transform 'Xstart', calculate segment lengths, and assign unique hash IDs
    df = transformXstart(df)
    df = segLength(df)
    df = giveHashID(df)

    # Calculate midpoints
    df = midPoint(df)

    # Calculate Hough Transform attributes (theta, rho, xc, yc) and perpendicular offset distances
    df, xc, yc = HoughTransform(df)
    df = MidtoPerpDistance(df, xc, yc)

    # Assign the processing date
    now = datetime.now()
    d = now.strftime("%d %b, %Y")
    df['Date Changed'] = [d] * len(df)

    # Remove duplicate entries and report any found duplicates
    l = len(df)
    df = df.drop_duplicates(subset=['HashID'])
    if l != len(df):
        print("Found", l - len(df), "duplicates")

    return df




def whichForm(lines):
    '''
    Returns the form of the dataframe column names

    Parameters:
        lines: a dataframe with columns containing theta, rho values
    Returns:
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

    if 'Average Rho (m)' in col:
        r='Average Rho (m)'
        t=r'Average Theta ($^\circ$)'

    return t,r

def MaskArea(df, bounds):

    """
    Returns dataframe masked by bounds
    Parameters
        df: pandas.dataframe with columns 'Xstart' and 'YStart'
        bounds: X and Y bounds in form [x1, y1, x2, y2]
        x1<x2 and y1<y2

    Returns:
        df_masked: returns all values for dataframe within those area bounds
    """

    maskX=( df['Xstart']>bounds[0]) & (df['Xstart']<bounds[2])
    maskY=( df['Ystart']>bounds[1]) & (df['Ystart']<bounds[3])

    masklatlong= (maskX==1) & (maskY==1)

    df_masked=df.loc[masklatlong]

    return df_masked

def FilterLines(lines):
    """
    Filters lines based on a trust filter.

    Parameters:
        lines (DataFrame): The input dataframe containing line data.

    Returns:
        DataFrame: A filtered dataframe containing only the lines marked as trusted (TrustFilter == 1).
    """
    mask = lines['TrustFilter'] == 1
    return lines[mask]


def getCartLimits(lines):
    """
    Computes the Cartesian limits (x and y) of a set of lines.

    Parameters:
        lines (DataFrame): The input dataframe containing line data.

    Returns:
        tuple: A tuple containing the x and y limits (xlim, ylim) as lists [min, max].
    """
    xlim = [np.min([lines['Xstart'].min(), lines['Xend'].min()]), np.max([lines['Xstart'].max(), lines['Xend'].max()])]
    ylim = [np.min([lines['Ystart'].min(), lines['Yend'].min()]), np.max([lines['Ystart'].max(), lines['Yend'].max()])]
    return xlim, ylim

def writeFile(df, name, myProj=None):
    """
    Writes a dataframe to a file based on the file extension.

    Parameters:
        df: (pandas.DataFrame) a pandas dataframe
        name: (string) the name of the file to be written with file extension

    Returns:
        df: (pandas.DataFrame) the input dataframe
    """

    # if file is not .csv, .txt, or .shp, return error
    if path.endswith('.csv') or path.endswith('.txt') or path.endswith('.shp') or path.endswith('.geojson') or path.endswith('.json'):
        raise ValueError("Invalid file type")

   # if ends with .csv or .txt, write as csv
    if name.endswith('.csv') or name.endswith('.txt'):
        df = writeToWKT(df, name, myProj=myProj)
    # if ends with .shp, write as shapefile
    else: 
        if name.endswith('.shp'):
            driver = 'ESRI Shapefile'
        elif name.endswith('.geojson') or name.endswith('.json'):
            driver = 'GeoJSON'
        elif name.endswith('.gpkg'):
            driver = 'GPKG'
        else:
            raise ValueError("Invalid file type")
        
        df = writetoGeoData(df, name, driver, myProj=myProj)


    return df
