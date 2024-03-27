# LinkingLines Package
 # Written by aikubo
 
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
r"""
htMOD Module

HT: provides functions for working with line segments and performing Hough Transform-related calculations.

The Hough Transform (HT) is a feature extraction technique used in image analysis and computer vision to detect lines. The Hough Transform algorithm 
works by transforming lines from Cartesian space to a "slope-intercept" space, where each line is represented by a point in the new space.
We use the algorithm from Ballard, D.H. (1981). Generalizing the Hough Transform to Detect Arbitrary Shapes. Pattern Recognition, 13(2), 111-122.

.. math :: \rho = x \cos(\theta) + y \sin(\theta)

where :math:`\rho` is the perpendicular distance from the origin to the line, and :math:`\theta` is the angle between the x-axis and the line.

.. math :: \theta = \arctan(-1/m) 

where "m" is the slope of the line.


Functions:
    - `CyclicAngleDistance(u,v)`: calculates smallest distance between angles
    - `HoughTransform(df, xc=None, yc=None)`: calculates Hough Transform
    - `HT_center(df)`: Calculates only HT center cartesian coordinates
    - `rotateData(df, rotation_angle)`: rotates cartesian dataframe by rotation_angle
    - `MidtoPerpDistance(df, xc, yc)`: calculates mid to perp distance
    - `moveHTcenter(df, xc=None, yc=None)`: moves HT center and recalcultes

"""

import pandas as pd
import numpy as np



def segLength(df):
    """
    Computes and adds a 'seg_length' column to a DataFrame, representing the length of line segments.

    Parameters
    ----------
    df : pandas.DataFrame
        The input DataFrame containing line data with 'Xstart', 'Xend', 'Ystart', and 'Yend' columns.

    Returns
    -------
    df : pandas.DataFrame
        The input DataFrame with an additional 'seg_length' column representing the length of line segments.

    """
    if not isinstance(df, pd.DataFrame):
        raise ValueError("Input 'df' must be a pandas DataFrame.")
    if not all(column in df.columns for column in ['Xstart', 'Ystart', 'Xend', 'Yend']):
        raise ValueError("Input DataFrame must contain columns 'Xstart', 'Ystart', 'Xend', and 'Yend'.")
    
    length = np.sqrt((df['Xstart'] - df['Xend'])**2 + (df['Ystart'] - df['Yend'])**2)
    df['seg_length'] = length
    return df

def CyclicAngleDist(u, v):
    """
    Calculate the cyclic angle distance between two angles in degrees.

    This function calculates the cyclic angle distance between two angles in
    degrees, considering the cyclical nature of angles.
    The result represents the smallest angle difference between the two input
    angles, ranging from 0 to 90 degrees.


    Parameters
    ----------
    u : list or array of floats
        The first angle(s) in degrees.
    v : list or array of floats
        The second angle(s) in degrees.

    Returns
    -------
    dist : list or array of floats
        The cyclic angle distance between the two angles, ranging from 0 to 90 degrees.


    Example
    -------
    >>> angle1 = [45.0]
    >>> angle2 = [160.0]
    >>> distance = CyclicAngleDist(angle1, angle2)
    >>> print("Cyclic Angle Distance (degrees):", distance)

    """

    u[0]=u[0]+90
    v[0]=v[0]+90
    return min( (u[0]-v[0])%180, (v[0]-u[0])%180)



def HoughTransform(df, xc=None, yc=None):
    """
    Calculates the Hough Transform of a dataframe of line segments.

    Parameters
    ----------
    df : pandas.Dataframe
        dataframe of the line segments must contain ["Xstart", "Ystart", "Xend", "Yend"]
    xc : float, optional
        x location of the HT center.
        If none is given, the center is calculated from the dataframe.
    yc : float, optional
        y location of the HT center.
        If none is given, the center is calculated from the dataframe.

    Returns
    -------
    newdf : pandas.Dataframe
        line segments with new columns of ['theta', 'rho', 'xc', 'yc']

    """

    # Check if 'df' is a pandas DataFrame
    if not isinstance(df, pd.DataFrame):
        raise ValueError("Input 'df' must be a pandas DataFrame.")

    if len(df) < 1:
        raise ValueError("DataFrame is empty")

    df=segLength(df)

    if any(np.isclose(df['seg_length'],0)):
        raise ValueError('Some lines are points')
    
    if not all(column in df.columns for column in ['Xstart', 'Ystart', 'Xend', 'Yend']):
        raise ValueError("Input DataFrame must contain columns 'Xstart', 'Ystart', 'Xend', and 'Yend'.")

    


    if xc is None or yc is None:
        xc,yc=HT_center(df)


    o_x1 = df['Xstart'].values - xc
    o_x2 = df['Xend'].values - xc
    o_y1 = df['Ystart'].values - yc
    o_y2 = df['Yend'].values - yc

    A = o_y1.astype(float) - o_y2.astype(float) +0.00000000001
    B = o_x1.astype(float) - o_x2.astype(float) +0.00000000001

    m=A/B
    angle=np.arctan(-1/m)
    b1=-1*m*o_x1+o_y1
    rho=b1*np.sin(angle)
    theta=np.rad2deg(angle)

    newdf= df.copy()
    newdf['theta']=theta
    newdf['rho']=rho
    newdf['xc']=xc
    newdf['yc']=yc

    return newdf, xc, yc


def HT_center(df):
    """
    Finds the center of a dataframe of line segments.

    Parameters
    ----------
    df : pandas.Dataframe
        dataframe of the line segments
        must contain ["Xstart", "Ystart", "Xend", "Yend"]

    Returns
    -------
    xc : float
        x location of the HT center.
    yc : float
        y location of the HT center.

    """
    
    if not all(column in df.columns for column in ['Xstart', 'Ystart', 'Xend', 'Yend']):
        raise ValueError("Input DataFrame must contain columns 'Xstart', 'Ystart', 'Xend', and 'Yend'.")
    

    xc=np.mean( (df['Xstart'].values+df['Xend'].values)/2)
    yc=np.mean( (df['Ystart'].values+df['Yend'].values)/2)

    return xc,yc

def rotateData(df, rotation_angle, xc = None, yc = None):

    """
    Rotates a dataframe of line segments by a given angle.

    Parameters
    ----------
    df : pandas.Dataframe
        line segments dataframe 
    rotation_angle : float
        angle of rotation in degrees
    xc : float, optional
        x location of the HT center.
        If none is given, the center is calculated from the dataframe.
    yc : float, optional
        y location of the HT center.
        If none is given, the center is calculated from the dataframe.


    Returns
    -------
    dfRotated: pandas.Dataframe
        dataframe of the line segments rotated by the given angle
    """


    if xc is None or yc is None:
        xc,yc=HT_center(df)

    rotation_angle=float(rotation_angle)
    o_x1 = df['Xstart'].values - xc
    o_x2 = df['Xend'].values - xc
    o_y1 = df['Ystart'].values - yc
    o_y2 = df['Yend'].values - yc
    #print(xc,yc)

    x1 = o_x1*np.cos(np.deg2rad(rotation_angle)) - o_y1*np.sin(np.deg2rad(rotation_angle))
    y1 = o_x1*np.sin(np.deg2rad(rotation_angle)) + o_y1*np.cos(np.deg2rad(rotation_angle))

    x2 = o_x2*np.cos(np.deg2rad(rotation_angle)) - o_y2*np.sin(np.deg2rad(rotation_angle))
    y2 = o_x2*np.sin(np.deg2rad(rotation_angle)) + o_y2*np.cos(np.deg2rad(rotation_angle))

    ang=np.deg2rad(rotation_angle)

    dfRotated = df.copy(deep=True)

    dfRotated['Xstart']=x1+xc
    dfRotated['Ystart']=y1+yc
    dfRotated['Xend']=x2+xc
    dfRotated['Yend']=y2+yc

    #xcR, ycR=HT_center(dfRotated)
    dfRotated, _, _=HoughTransform(dfRotated, xc=xc, yc=yc)


    #print(xcR,ycR)

    return dfRotated





def MidtoPerpDistance(df, xc=None, yc=None):
    """
    Find the distance between line segment midpoint and rho line perpendicular
    intersection.

    Parameters
    ----------
    df: pandas.Dataframe
        dataframe of the line segments
        must contain ["Xstart", "Ystart", "Xend", "Yend"]
    xc : float, optional
        x location of the HT center.
        If none is given, the center is calculated from the dataframe.
    yc : float, optional
        y location of the HT center.
        If none is given, the center is calculated from the dataframe.


    Returns
    -------
    df: pandas.Dataframe
    with new columns of ['PerpOffDist', 'PerpIntX', 'PerpIntY']
    """

    try:
        'Xmid' not in df.columns
    except:
        print("Xmid must be calculate first")
    
    if not all(column in df.columns for column in ['Xstart', 'Ystart', 'Xend', 'Yend']):
        raise ValueError("Input DataFrame must contain columns 'Xstart', 'Ystart', 'Xend', and 'Yend'.")


    df, xc, yc=HoughTransform(df, xc=xc, yc=yc)
    rho = df['rho'].values 
    theta = df['theta'].values

    intx=rho*np.cos(np.deg2rad(df['theta'].values))
    inty=rho*np.sin(np.deg2rad(df['theta'].values))
    df['PerpOffsetDist']=np.sqrt( (df['Xmid'].values-intx)**2 +  (df['Ymid'].values-inty)**2)*np.sign((df['Ymid'].values-inty))
    df['PerpIntX']=intx
    df['PerpIntY']=inty

    return df

def moveHTcenter(df, xc=None, yc=None):
    """
    Move the center of a DataFrame of line segments to new coordinates (xc, yc).

    Parameters
    ----------
    df: pandas.DataFrame
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
    >>> df = pd.DataFrame({'Xstart': [0, 1, 2], 'Ystart': [0, 1, 2], 'Xend': [1, 2, 3], 'Yend': [1, 2, 3]})
    >>> new_df = moveHTcenter(df, xc=2.0, yc=2.0)
    >>> print(new_df)
       Xstart  Ystart  Xend  Yend
    0    -2.0    -2.0  -1.0  -1.0
    1    -1.0    -1.0   0.0   0.0
    2     0.0     0.0   1.0   1.0
    """
    if not all(column in df.columns for column in ['Xstart', 'Ystart', 'Xend', 'Yend']):
        raise ValueError("Input DataFrame must contain columns 'Xstart', 'Ystart', 'Xend', and 'Yend'.")
    
    if xc is None and yc is None:
        xc,yc=HT_center(df)

    o_x1 = df['Xstart'].values - xc
    o_x2 = df['Xend'].values - xc
    o_y1 = df['Ystart'].values - yc
    o_y2 = df['Yend'].values - yc

    df2=pd.DataFrame({'Xstart':o_x1, 'Ystart':o_y1, 'Xend':o_x2, 'Yend':o_y2})

    return df2
