# LinkingLines Package
 # Written by aikubo
 # Version: 2.1.0
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
htMOD Module

This module provides functions for working with line segments and performing Hough Transform-related calculations.

Functions:
    - CyclicAngleDistance(u,v): calculates smallest distance between angles
    - HoughTransform(data, xc=None, yc=None): calculates Hough Transform
    - HT_center(data): Calculates only HT center cartesian coordinates
    - rotateData(data, rotation_angle): rotates cartesian data by rotation_angle
    - MidtoPerpDistance(df, xc, yc): calculates mid to perp distance
    - moveHTcenter(data, xc=None, yc=None): moves HT center and recalcultes

Notes:
    - This module contains functions for calculating the Hough Transform of line segments, rotating line segments,
      finding the center of line segments, calculating distances, and moving the center of line segments.

Author:
    Created by akh on Thu Apr 1 12:49:07 2021.

Dependencies:
    - pandas
    - numpy
    - matplotlib.pyplot
"""

import pandas as pd
import numpy as np



def segLength(df):
    """
    Computes and adds a 'seg_length' column to a DataFrame, representing the length of line segments.

    Parameters:
        df (DataFrame): The input DataFrame containing line data with 'Xstart', 'Xend', 'Ystart', and 'Yend' columns.

    Returns:
        DataFrame: The input DataFrame with an additional 'seg_length' column representing the length of line segments.
    """
    length = np.sqrt((df['Xstart'] - df['Xend'])**2 + (df['Ystart'] - df['Yend'])**2)
    df['seg_length'] = length
    return df

def CyclicAngleDist(u, v):
    """
    Calculate the cyclic angle distance between two angles in degrees.

    Parameters:
        u (list or array): The first angle(s) in degrees.
        v (list or array): The second angle(s) in degrees.

    Returns:
        dist (list or array): The cyclic angle distance between the two angles, ranging from 0 to 90 degrees.

    This function calculates the cyclic angle distance between two angles in
    degrees, considering the cyclical nature of angles.
    The result represents the smallest angle difference between the two input
    angles, ranging from 0 to 90 degrees.

    Example usage:
    angle1 = [45.0]
    angle2 = [160.0]
    distance = CyclicAngleDist(angle1, angle2)
    print("Cyclic Angle Distance (degrees):", distance)
    """
    u[0]=u[0]+90
    v[0]=v[0]+90
    return min( (u[0]-v[0])%180, (v[0]-u[0])%180)



def HoughTransform(data, xc=None, yc=None):
    """
    Calculates the Hough Transform of a dataframe of line segments.

    Parameters
    ----------
    data: pandas.Dataframe
        dataframe of the line segments
        must contain ["Xstart", "Ystart", "Xend", "Yend"]

    xc, yc: float, optional
        x and y location of the center of the HT.
        If none is given, the center is calculated from the dataframe.

    Returns
    -------
    theta: HT angle in degrees
    rho: HT rho in same unit as 'Xstart'

    """

    # Check if 'data' is a pandas DataFrame
    if not isinstance(data, pd.DataFrame):
        raise ValueError("Input 'data' must be a pandas DataFrame.")

    if len(data) < 1:
        raise ValueError("DataFrame is empty")

    data=segLength(data)

    if any(np.isclose(data['seg_length'],0)):
        raise ValueError('Some lines are points')


    if xc is None or yc is None:
        xc,yc=HT_center(data)


    o_x1 = data['Xstart'].values - xc
    o_x2 = data['Xend'].values - xc
    o_y1 = data['Ystart'].values - yc
    o_y2 = data['Yend'].values - yc

    A = o_y1.astype(float) - o_y2.astype(float) +0.00000000001
    B = o_x1.astype(float) - o_x2.astype(float) +0.00000000001

    m=A/B
    angle=np.arctan(-1/m)
    b1=-1*m*o_x1+o_y1
    rho=b1*np.sin(angle)
    theta=np.rad2deg(angle)

    return theta, rho, xc, yc


def HT_center(data):
    """
    Finds the center of a dataframe of line segments.

    Parameters
    ----------
    df: pandas.Dataframe
        dataframe of the line segments
        must contain ["Xstart", "Ystart", "Xend", "Yend"]

    Returns
    -------
    xc, yc: float
        x and y location of the center of the HT

    """
    xc=np.mean( (data['Xstart'].values+data['Xend'].values)/2)
    yc=np.mean( (data['Ystart'].values+data['Yend'].values)/2)

    return xc,yc

def rotateData(data, rotation_angle):

    """
    Rotates a dataframe of line segments by a given angle.

    Parameters
    ----------
    x,y : numpy.Array
        array of endpoints such as output from fitRectangle/endpoints2

    rotation_angle: float
        angle of rotation in degrees

    xc, yc: float
        x and y location of the center of the HT.


    Returns
    -------
    df: pandas.Dataframe
        dataframe of the line segments rotated by the given angle
    """



    xc,yc=HT_center(data)

    rotation_angle=float(rotation_angle)
    o_x1 = data['Xstart'].values - xc
    o_x2 = data['Xend'].values - xc
    o_y1 = data['Ystart'].values - yc
    o_y2 = data['Yend'].values - yc
    #print(xc,yc)

    x1 = o_x1*np.cos(np.deg2rad(rotation_angle)) - o_y1*np.sin(np.deg2rad(rotation_angle))
    y1 = o_x1*np.sin(np.deg2rad(rotation_angle)) + o_y1*np.cos(np.deg2rad(rotation_angle))

    x2 = o_x2*np.cos(np.deg2rad(rotation_angle)) - o_y2*np.sin(np.deg2rad(rotation_angle))
    y2 = o_x2*np.sin(np.deg2rad(rotation_angle)) + o_y2*np.cos(np.deg2rad(rotation_angle))

    ang=np.deg2rad(rotation_angle)

    dataRotated=data.copy()

    dataRotated['Xstart']=x1+xc
    dataRotated['Ystart']=y1+yc
    dataRotated['Xend']=x2+xc
    dataRotated['Yend']=y2+yc

    #xcR, ycR=HT_center(dataRotated)
    theta, rho, xc, yc=HoughTransform(dataRotated, xc=xc, yc=yc)
    dataRotated['theta']=theta
    dataRotated['rho']=rho


    #print(xcR,ycR)

    return dataRotated





def MidtoPerpDistance(df, xc, yc):
    """
    Find the distance between line segment midpoint and rho line perpendicular
    intersection.

    Parameters
    ----------
    df: pandas.Dataframe
        dataframe of the line segments
        must contain ["Xstart", "Ystart", "Xend", "Yend"]
    xc, yc: float
        x location of HT center

    Returns
    -------
    df: pandas.Dataframe
    with new columns of ['PerpOffDist', 'PerpIntX', 'PerpIntY']
    """

    try:
        'Xmid' not in df.columns
    except:
        print("Xmid must be calculate first")
    # if 'Xmid' not in df.columns:
    #     df=midPoint(df)


    theta,rho,xc,yc=HoughTransform(df, xc, yc)
    intx=rho*np.cos(np.deg2rad(theta))
    inty=rho*np.sin(np.deg2rad(theta))
    df['PerpOffsetDist']=np.sqrt( (df['Xmid'].values-intx)**2 +  (df['Ymid'].values-inty)**2)*np.sign((df['Ymid'].values-inty))
    df['PerpIntX']=intx
    df['PerpIntY']=inty

    return df

def moveHTcenter(data, xc=None, yc=None):
    """
    Move the center of a DataFrame of line segments to new coordinates (xc, yc).

    Parameters
    ----------
    data: pandas.DataFrame
        A DataFrame containing line segments with columns ['Xstart', 'Ystart', 'Xend', 'Yend'].

    xc, yc: float, optional
        The new x and y coordinates of the center of the line segments.
        If 'xc' and 'yc' are not provided, the center is calculated from the input DataFrame.

    Returns
    -------
    pandas.DataFrame
        A new DataFrame with line segments adjusted to the new center coordinates (xc, yc).

    Notes
    -----
    - This function takes a DataFrame of line segments and moves their center to the specified coordinates (xc, yc).
    - If 'xc' and 'yc' are not provided, the center is calculated based on the input DataFrame.
    - The resulting DataFrame 'df2' contains line segments with updated coordinates.

    Example
    -------
    >>> data = pd.DataFrame({'Xstart': [0, 1, 2], 'Ystart': [0, 1, 2], 'Xend': [1, 2, 3], 'Yend': [1, 2, 3]})
    >>> new_data = moveHTcenter(data, xc=2.0, yc=2.0)
    >>> print(new_data)
       Xstart  Ystart  Xend  Yend
    0    -2.0    -2.0  -1.0  -1.0
    1    -1.0    -1.0   0.0   0.0
    2     0.0     0.0   1.0   1.0
    """

    if xc is None and yc is None:
        xc,yc=HT_center(data)

    o_x1 = data['Xstart'].values - xc
    o_x2 = data['Xend'].values - xc
    o_y1 = data['Ystart'].values - yc
    o_y2 = data['Yend'].values - yc

    df2=pd.DataFrame({'Xstart':o_x1, 'Ystart':o_y1, 'Xend':o_x2, 'Yend':o_y2})

    return df2
