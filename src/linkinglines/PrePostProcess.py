# LinkingLines Package
 # Written by aikubo

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
    preProcess:     Fully preprocesses a dataframe containing line data to ensure it has essential attributes and is properly formatted.
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
import os
import geopandas

def readFile(name, preprocess=True):
    """
    Reads in a file and returns a pandas dataframe.

    Parameters
    ----------
    name : str
        The path to the file to be read in.

    preprocess : bool, optional
        Indicates whether to preprocess the data, by default True.

    Returns
    -------
    data : pandas.DataFrame
        A pandas or geopandas dataframe containing the read data.

    """
    # if not a valid path, return error
    if not os.path.exists(name):
        raise ValueError("Invalid path")
    # if file is not .csv, .txt, or .shp, return error
    valid_extensions = ['.csv', '.txt', '.shp', '.geojson', '.json']

    if not any(name.endswith(ext) for ext in valid_extensions):
        raise ValueError("Invalid file type")

    # identify the type of file
    # read in .csv
    if name.endswith('.csv'):
        data=pd.read_csv(name)
    elif name.endswith('.txt'):
        data=pd.read_csv(name, delimiter='\t')
    else:
        data=geopandas.read_file(name)
        data=data.to_wkt()

    if 'Xstart' not in data.columns:
        data = WKTtoArray(data)

    # if preprocess is True, preprocess the data
    if preprocess:
        data = preProcess(data)

    return data


def midPoint(df):
    """
    Finds the midpoint of a dataframe of line segments.

    Parameters
    ----------
    df : pandas.DataFrame
        Dataframe containing the line segments with columns ["Xstart", "Ystart", "Xend", "Yend"].

    Returns
    -------
    df : pandas.DataFrame
        Dataframe with new columns ['Xmid', 'Ymid'] indicating the midpoints of the line segments.
    """

    df['Xmid']=(df['Xstart']+df['Xend'])/2
    df['Ymid']=(df['Ystart']+df['Yend'])/2

    return df


def writeToWKT(df,name, myProj=None):
    """
    Writes a dataframe to a CSV file using Well-Known Text (WKT) format for line vectors.

    Parameters
    ----------
    df : pandas.DataFrame
        Dataframe with columns ['Xstart', 'Ystart', 'Xend', 'Yend', 'seg_length'].
    name : str
        Name of the output file.
    myProj : pyproj.CRS, optional
        Projection of the dataframe. If None, WGS84 is assumed, by default None.

    Returns
    -------
    pandas.DataFrame
        The input dataframe with an additional 'Linestring' column.
    """

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
    Writes a dataframe to a geospatial file (shapefile, geopackage, or geojson).

    Parameters
    ----------
    df : pandas.DataFrame
        Dataframe with geospatial data to be written.
    name : str
        Name of the output file.
    driver : str
        Format of the file (e.g., 'ESRI Shapefile', 'GeoJSON').
    myProj : pyproj.CRS, optional
        Projection of the dataframe. If None, WGS84 is assumed, by default None.

    Returns
    -------
    pandas.DataFrame
        The input dataframe.
    """

    gdf = geopandas.GeoDataFrame(df, geometry=geopandas.points_from_xy(df.Xstart, df.Ystart))
    gdf.to_file(name, driver=driver, crs=myProj)

    return df


def WKTtoArray(df, plot=False):
    """
    Processes a dataframe with WKT strings to a pandas dataframe with explicit geometry columns.

    Parameters
    ----------
    df : pandas.DataFrame
        Dataframe with a 'WKT' or 'geometry' column.
    plot : bool, optional
        Whether to plot the processed lines, by default False.

    Returns
    -------
    pandas.DataFrame
        Dataframe with columns ['Xstart', 'Ystart', 'Xend', 'Yend', 'seg_length'].

    """
    if not isinstance(df, pd.DataFrame):
        raise ValueError("Input 'data' must be a pandas DataFrame.")

    if len(df) < 1:
        raise ValueError("DataFrame is empty")

    #     # if neither is in columns raise value error
    if not ("WKT" in df.columns ):
        if not ("geometry" in df.columns):
         raise ValueError("No geometry present")

    xstart=[]
    ystart=[]

    xend=[]
    yend=[]
    drop=[]

    if plot:
        fig,ax=plt.subplots()

    # check for either "WKT" or "geometry" columns
    if "geometry" in df.columns:
        tag = "geometry"
    else:
        tag = "WKT"

    for i in range(len(df)):
        temp=df[tag].iloc[i]
        t1=temp[0]
        # Using regex to find all numbers in the string
        temp = re.findall(r"[-+]?\d*\.\d+|\d+", temp)


        if len(temp)<1:
            drop.append(i)
            continue

        if "Z" in t1:
            tempx=np.array(temp[::3]).astype(float)
            tempy=np.array(temp[1::3]).astype(float)
        else:
            tempx=np.array(temp[::2]).astype(float)
            tempy=np.array(temp[1::2]).astype(float)

        if np.unique(tempx).shape[0]==1 or np.unique(tempy).shape[0]==1:
            drop.append(i)
            continue


        slope, intercept, r_value, p_value, std_err = stats.linregress(tempx, tempy)

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
    Assigns a numeric ID to the dataframe rows.

    Parameters
    ----------
    df : pandas.DataFrame
        The input dataframe.

    Returns
    -------
    pandas.DataFrame
        The dataframe with an added 'ID' column.
    """
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
    Assigns a hash ID to a dataframe based on line endpoints.

    Parameters
    ----------
    df : pandas.DataFrame
        Dataframe with line segment data.

    Returns
    -------
    pandas.DataFrame
        Dataframe with an added 'HashID' column based on the hash of line endpoints.
    """

    ids=[]
    for i in range(len(df)):
        h=pd.util.hash_pandas_object(df.iloc[i])
        ids.append(hash( (h['Xstart'], h['Ystart'], h['Xend'], h['Yend'])))
    df['HashID']=ids

    return df




def transformXstart(dikeset, HTredo=True):
    """
    Ensures that 'Xstart' is always less than 'Xend' in a dataframe of line segments.

    Parameters
    ----------
    df : pandas.DataFrame
        Dataframe with line segments.

    Returns
    -------
    df : pandas.DataFrame
        Transformed dataframe where 'Xstart' < 'Xend'.
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

    Parameters
    ----------
    df : pandas.DataFrame
        Input dataframe containing dike line data.
    HTredo : bool, default True
        Whether to recalculate Hough Transform attributes, by default True.
    xc : float, optional
        X-coordinate of the center point for the Hough Transform, by default None.
    yc : float, optional
        Y-coordinate of the center point for the Hough Transform, by default None.

    Returns
    -------
    df : pandas.DataFrame
        Processed dataframe with added or updated attributes.
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

    Parameters
    ----------
    df : pandas.DataFrame
        The input dataframe containing line data.
    HTredo : bool, default True
        Indicates whether to recalculate Hough Transform attributes, by default True.

    Returns
    -------
    pandas.DataFrame
        The processed dataframe with added or updated attributes.
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
    """
    Fully preprocesses a dataframe containing line data to ensure it has essential attributes and is properly formatted.
        1. Convert WKT column to array format if present
        2. Transform Xstart, so Xstart < Xend
        3. calculate segment lengths
        4. assign unique hash IDs
        5. calculate midpoints
        6. calculate Hough Transform attributes (theta, rho, xc, yc)
        7. calculate perpendicular offset distances
        8. assign the processing date
        9. remove duplicate entries and report any found duplicates

    Parameters
    ----------
    data : pandas.DataFrame
        The input dataframe containing line data.

    Returns
    -------
    pandas.DataFrame
        The fully processed dataframe with all required attributes.
    """
    df = data.copy()

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
    """
    Identifies the form of the dataframe column names.

    Parameters
    ----------
    lines : pandas.DataFrame
        A dataframe with columns containing theta and rho values.

    Returns
    -------
    tuple
        A tuple (t, r) containing the string identifiers for theta and rho columns.
    """

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
    Masks a dataframe by specified bounds.

    Parameters
    ----------
    df : pandas.DataFrame
        A dataframe with columns 'Xstart' and 'Ystart'.
    bounds : list or tuple
        The bounding box specified as [x1, y1, x2, y2], where x1<x2 and y1<y2.

    Returns
    -------
    pandas.DataFrame
        A masked dataframe containing only the rows within the specified bounds.
    """

    maskX=( df['Xstart']>bounds[0]) & (df['Xstart']<bounds[2])
    maskY=( df['Ystart']>bounds[1]) & (df['Ystart']<bounds[3])

    masklatlong= (maskX==1) & (maskY==1)

    df_masked=df.loc[masklatlong]

    return df_masked

def FilterLines(lines):
    """
    Filters lines based on a trust filter.

    Parameters
    ----------
    lines : pandas.DataFrame
        The input dataframe containing line data.

    Returns
    -------
    pandas.DataFrame
        A filtered dataframe containing only the lines marked as trusted (TrustFilter == 1).
    """
    mask = lines['TrustFilter'] == 1
    return lines[mask]


def getCartLimits(lines):
    """
    Computes the Cartesian limits (x and y) of a set of lines.

    Parameters
    ----------
    lines : pandas.DataFrame
        The input dataframe containing line data.

    Returns
    -------
    xlim : float
        The limits of the x-axis.
    ylim : float
        The limits of the y-axis.
    """
    xlim = [np.min([lines['Xstart'].min(), lines['Xend'].min()]), np.max([lines['Xstart'].max(), lines['Xend'].max()])]
    ylim = [np.min([lines['Ystart'].min(), lines['Yend'].min()]), np.max([lines['Ystart'].max(), lines['Yend'].max()])]
    return xlim, ylim

def writeFile(df, name, myProj=None):
    """
    Writes a dataframe to a file based on the file extension.

    Parameters
    ----------
    df : pandas.DataFrame
        A pandas dataframe to be written to file.
    name : str
        The name of the file to be written, including the file extension.
    myProj : str, optional
        The projection of the dataframe, if applicable. Default is None.

    Returns
    -------
    pandas.DataFrame
        The input dataframe.
    """
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
            driver = 'GeoJSON'

        df = writetoGeoData(df, name, driver, myProj=myProj)


    return df
