# LinkingLines Package
 # Written by aikubo
 # Version: 2.1.0
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 10 13:30:34 2021

@author: akh


fitRectangle: A Module for Fitting Rectangles to Rotated Line Segments

The 'fitRectangle' module provides functions for rotating line segments and fitting rectangles around them.
It includes utilities for calculating rectangle properties and performing operations on line segments.

Functions:
- `rotateXYShift(ang, x, y, h, k)`: Rotate and shift coordinates (x, y) about a center point (h, k) by an angle 'ang'.
- `unrotate(ang, x, y, h, k)`: Reverse the rotation and shift of coordinates (x, y) about a center point (h, k) by an angle 'ang'.
- `endpoints(lines)`: Extract and return the x and y coordinates of endpoints from a DataFrame of line segments.
- `midpoint(lines)`: Calculate and return the x and y coordinates of midpoints for line segments in a DataFrame.
- `allpoints(lines)`: Calculate and return all the x and y coordinates along the line segments in a DataFrame.
- `fit_Rec(lines, xc, yc)`: Fit a rectangle to a cluster of lines and return its width, length, correlation coefficient, and coordinates.
- `RecEdges(xi, yi, avgtheta, x0, y0)`: Calculate the coordinates of the edges of a rectangle based on parameters.
- `pltLine(lines, xc, yc, ax)`: Plot a line representing the fitted rectangle on a specified matplotlib axis.
- `W_L(Clusters)`: Calculate the widths and lengths of rectangles fitted to clusters of lines.
- `squaresError(lines, xc, yc)`: Calculate the sum of squared errors for a cluster of lines fitted to a rectangle.


"""

import numpy as np
import matplotlib.pyplot as plt


def rotateXYShift(ang,x,y,h,k):
    """
  Rotate and shift coordinates (x, y) about a center point (h, k) by an angle 'ang' (in radians).

  Parameters
  ----------
  ang : float
      The angle in radians by which to rotate the coordinates.
  x : float or numpy array
      The x-coordinate of the point to be transformed.
  y : float or numpy array
      The y-coordinate of the point to be transformed.
  h : float
      The x-coordinate of the center point about which the rotation and shift will be performed.
  k : float
      The y-coordinate of the center point about which the rotation and shift will be performed.

  Returns
  -------
  xp : float or numpy array
    The transformed x-coordinate after rotating and shifting (x, y) about (h, k).
  yp : float or numpy array
    The transformed y-coordinate after rotating and shifting (x, y) about (h, k).

  Notes
  -----
  The function rotates the point (x, y) counterclockwise by 'ang' radians around the center point (h, k).

  """

    xp= (x-h)*np.cos(ang)-(y-k)*np.sin(ang)
    yp= (x-h)*np.sin(ang)+(y-k)*np.cos(ang)
    return xp, yp


def unrotate(ang, x, y, h, k):
    """
    Reverse the rotation and shift of coordinates (x, y) about a center point (h, k) by an angle 'ang' (in radians).

    Parameters
    ----------
    ang : float
        The angle in radians by which the coordinates were previously rotated.

    x : float or numpy array
        The x-coordinate of the transformed point.

    y : float or numpy array
        The y-coordinate of the transformed point.

    h : float
        The x-coordinate of the center point about which the reverse rotation and shift will be performed.

    k : float
        The y-coordinate of the center point about which the reverse rotation and shift will be performed.

    Returns
    -------
    float, float or numpy array, numpy array
        The original coordinates (xr, yr) before the rotation and shift were applied.

    Notes
    -----
    - The function reverses the rotation and shift applied to the point (x, y) by 'ang' radians around the center point (h, k).
    - It returns the original coordinates (xr, yr) before the transformation.

    See Also
    --------
    rotateXYShift : Rotate and shift coordinates (x, y) about a center point (h, k) by an angle 'ang'.

    """
    xr = x * np.cos(ang) + y * np.sin(ang) + h
    yr = -1 * (x) * np.sin(ang) + y * np.cos(ang) + k
    return xr, yr


def endpoints(lines):
    """
    Extracts and returns the x and y coordinates of endpoints from a DataFrame of line segments.

    Parameters
    ----------
    lines : pandas.DataFrame
        A DataFrame containing line segments with columns ['Xstart', 'Ystart', 'Xend', 'Yend'].

    Returns
    -------
    xlist : numpy.ndarray
        A NumPy array containing the x coordinates of endpoints.
    ylist : numpy.ndarray
        A NumPy array containing the y coordinates of endpoints.

    """
    xlist = np.array([])
    ylist = np.array([])

    for i in range(0, len(lines)):
        x1 = lines['Xstart'].iloc[i]
        y1 = lines['Ystart'].iloc[i]
        y2 = lines['Yend'].iloc[i]
        x2 = lines['Xend'].iloc[i]

        xlist = np.append(xlist, x1)
        xlist = np.append(xlist, x2)
        ylist = np.append(ylist, y1)
        ylist = np.append(ylist, y2)

    return xlist, ylist

def midpoint(lines):
    """
    Calculate and return the x and y coordinates of midpoints for line segments in a DataFrame.

    Parameters
    ----------
    lines : pandas.DataFrame
        A DataFrame containing line segments with columns ['Xstart', 'Ystart', 'Xend', 'Yend'].

    Returns
    -------
    xlist : numpy.ndarray
        A NumPy array containing the x coordinates of endpoints.
    ylist : numpy.ndarray
        A NumPy array containing the y coordinates of endpoints.
    """

    xstart=lines['Xstart'].to_numpy()
    ystart=lines['Ystart'].to_numpy()
    xend=lines['Xend'].to_numpy()
    yend=lines['Yend'].to_numpy()
    xlist=np.array([])
    ylist=np.array([])

    for i in range(len(xstart)):
        xs=np.mean([xstart[i], xend[i]])
        ys=np.mean([ystart[i], yend[i]])

        xlist=np.append(xlist, xs)
        ylist=np.append(ylist, ys)

    return xlist,ylist

def allpoints(lines):
    """
    Calculate and return all the x and y coordinates along the line segments in a DataFrame.

    This function calculates all the x and y coordinates that lie along the line segments in the DataFrame 'lines'.
    It evenly samples points along each line segment using linear interpolation.

    Parameters
    ----------
    lines : pandas.DataFrame
        A DataFrame containing line segments with columns ['Xstart', 'Ystart', 'Xend', 'Yend'].

    Returns
    -------
    xlist : numpy.ndarray
        A NumPy array containing the x coordinates of endpoints.
    ylist : numpy.ndarray
        A NumPy array containing the y coordinates of endpoints.


    """

    xstart=lines['Xstart'].to_numpy()
    ystart=lines['Ystart'].to_numpy()
    xend=lines['Xend'].to_numpy()
    yend=lines['Yend'].to_numpy()
    xlist=np.array([])
    ylist=np.array([])

    for i in range(len(xstart)):
        xs=np.linspace(xstart[i], xend[i])
        ys=np.linspace(ystart[i], yend[i])

        xlist=np.append(xlist, xs)
        ylist=np.append(ylist, ys)

    return xlist,ylist



def fit_Rec(lines, xc, yc):
    """
    Fits a rectangle to a cluster of lines and returns its width, length, correlation coefficient, and coordinates.

    Parameters
    ----------
    lines : pandas.DataFrame
        A DataFrame containing line segments with columns ['Xstart', 'Ystart', 'Xend', 'Yend', 'theta', 'rho'].
    xc : float
        x-coordinate of the center of the rectangle.
    yc : float
        y-coordinate of the center of the rectangle.

    Returns
    -------
    width : float
        The width of the fitted rectangle.
    length : float
        The length of the fitted rectangle.
    r : float
        The correlation coefficient of the fitted rectangle.
    xs : numpy.ndarray
        Array containing x coordinates of the edges of the fitted rectangle.
    ys : numpy.ndarray
        Array containing y coordinates of the edges of the fitted rectangle.
    Xmid : float
        The x-coordinate of the midpoint of the fitted rectangle.
    Ymid : float
        The y-coordinate of the midpoint of the fitted rectangle.
  

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
        t='Average Theta ($^\circ$)'

    if t=='AvgTheta' or t=='Average Theta ($^\circ$)':
        segl='R_Length'
    else:
        segl='seg_length'

    xi,yi=endpoints(lines)
    if len(lines) == 1:
        return 0, lines[segl], 0, xi,yi, lines['Xmid'], lines['Ymid']

    x0=xc
    y0=yc

    size=len(lines)

    if abs(np.sum(np.sign(lines[t].values))) < size:
        crossZero=True
        ang=np.mean(abs(lines[t].values))
        tol=6
        if np.isclose(ang,0, atol=4):
            ang=np.mean((lines[t].values))
    else:
        crossZero=False

        ang=np.mean((lines[t].values))



    xp, yp= rotateXYShift(np.deg2rad(-1*ang), xi,yi, x0,y0)

    width=np.ptp(xp.flatten())
    length=np.ptp(yp.flatten())

    xc=(max(xp)-min(xp))/2 + min(xp)
    yc=(max(yp)-min(yp))/2 + min(yp)

    xr=xc+width/2
    xl=xc-width/2
    yu=yc+length/2
    yd=yc-length/2
    xs=np.append(xr,xl)
    ys=np.append(yu,yd)


    xpi, ypi=unrotate(np.deg2rad(-1*ang), xp, yp, x0, y0)

    Xedges=np.array([xs[0], xs[0], xs[1], xs[1], xs[0]])
    Yedges=np.array([ys[1], ys[0], ys[0], ys[1], ys[1]])

    Xmid=(np.max(xs)+np.min(xs))/2
    Ymid=(np.max(ys)+np.min(ys))/2

    xs,ys=unrotate(np.deg2rad(-1*ang),Xedges,Yedges,x0,y0)
    Xmid, Ymid=unrotate(np.deg2rad(-1*ang), Xmid, Ymid, x0, y0)



    r=np.sum((yc-yp)**2)/lines[segl].sum()


    return width, length, r, xs, ys, Xmid, Ymid

def RecEdges(xi, yi, avgtheta, x0, y0):
    """
    Calculates the coordinates of the edges of a rectangle based on parameters.

    Parameters
    ----------
    xi : numpy.ndarray
        Array of x-coordinates of line endpoints.

    yi : numpy.ndarray
        Array of y-coordinates of line endpoints.

    avgtheta : float
        The average angle of the lines.

    x0 : float
        x-coordinate of the rectangle's center.

    y0 : float
        y-coordinate of the rectangle's center.

    Returns
    -------
    xs : numpy.ndarray
        Array containing x coordinates of the edges of the rectangle.
    ys : numpy.ndarray
        Array containing y coordinates of the edges of the rectangle.

    """

    ang=-1*np.deg2rad(avgtheta)

    xp, yp= rotateXYShift(ang, xi,yi, x0,y0)
    #plotlines(lines, 'k.-', a)

    width=np.ptp(xp)
    length=np.ptp(yp)

    # if width>length :
    #     length=width
    #     width=length
    xc=(max(xp)-min(xp))/2 + min(xp)
    yc=(max(yp)-min(yp))/2 + min(yp)

    xr=xc+width/2
    xl=xc-width/2
    yu=yc+length/2
    yd=yc-length/2
    xs=np.append(xr,xl)
    ys=np.append(yu,yd)


    xpi, ypi=unrotate(ang, xp, yp, x0, y0)

    Xedges=np.array([xs[0], xs[0], xs[1], xs[1], xs[0]])
    Yedges=np.array([ys[1], ys[0], ys[0], ys[1], ys[1]])
    # a.plot(Xedges, Yedges, 'r.-')

    xs,ys=unrotate(ang,Xedges,Yedges,x0,y0)


    return xs, ys


def pltLine(lines, xc, yc, ax):
    """
    Plots a line representing the fitted rectangle on an existing matplotlib axis.

    Parameters
    ----------
    lines : pandas.DataFrame
        A DataFrame containing line segments with columns ['Xstart', 'Ystart', 'Xend', 'Yend', 'theta', 'rho'].

    xc : float
        x-coordinate of the center of the rectangle.

    yc : float
        y-coordinate of the center of the rectangle.

    ax : matplotlib.axes.Axes
        The axis on which the line will be plotted.

    Returns
    -------
    None


    """

    avgtheta=np.deg2rad(np.average(lines['theta']))
    avgrho=np.average(lines['rho'])
    xs,ys=allpoints(lines)

    m=-np.cos(avgtheta)/np.sin(avgtheta)
    b=avgrho/np.sin(avgtheta)

    ax.plot( [min(xs), max(xs)], [m*min(xs-xc)+b+yc, m*max(xs-xc)+b+yc], 'orange', '..')


def W_L(Clusters):
    """
    Calculates the widths and lengths of rectangles fitted to clusters of lines.

    Parameters
    ----------
    Clusters : pandas.DataFrame
        A DataFrame containing clusters of lines with columns ['Labels', 'theta', 'rho', 'seg_length'].

    Returns
    -------
    width : numpy.ndarray
        Array containing widths of fitted rectangles for each cluster.
    length : numpy.ndarray
        Array containing lengths of fitted rectangles for each cluster.
    """

    width=np.array([])
    length=np.array([])

    labels=np.unique(Clusters['Labels'])
    for i in labels:
        mask=(Clusters['Labels']==i)
        lines=Clusters[mask]
        w,l=fit_Rec(lines)
        width=np.append(width,w)
        length=np.append(length,l)

    return width, length

def squaresError(lines, xc, yc):
    """
    Calculates the sum of squared errors for a cluster of lines fitted to a rectangle.

    Parameters
    ----------
    lines : pandas.DataFrame
        A DataFrame containing line segments with columns ['Xstart', 'Ystart', 'Xend', 'Yend', 'theta', 'seg_length'].

    xc : float
        x-coordinate of the center of the rectangle.

    yc : float
        y-coordinate of the center of the rectangle.

    Returns
    -------
    r : float
        The sum of squared errors.

    """


    avgtheta=np.deg2rad(np.average(lines['theta']))
    #print(avgtheta)
    avgrho=np.average(lines['rho'])
    xs,ys=midpoint(lines)

    m=-np.cos(avgtheta)/np.sin(avgtheta)
    b=avgrho/np.sin(avgtheta)

    r=np.sum((ys-(m*(xs-xc)+b+yc))**2)/lines['seg_length'].sum() #len(lines)


    return r
