#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
The provided Python module appears to be a collection of functions for analyzing
 and processing dike lines data in a geospatial context. Here's a summary of the
 module's main functionalities and the purpose of each function:

- `CenterFunc(theta, xr, yr, xc, yc)`: Calculates the radial distance of a point from a given center for a specified angle.

- `RadialFit(lines, plot=False, ColorBy=None, weight='LayerNumber', ThetaRange=[-90, 90], xc=None, yc=None)`: Fits a radial model to dike lines data, optionally visualizing the fit. It returns center information such as coordinates and goodness of fit.

- `RadialAzimuthal(lines, Center)`: Calculates radial azimuthal angles of dike lines relative to a specified center.

- `CyclicAngle360(u, v)`: Computes the cyclic angle between two angles in the range [0, 360).

- `AngleSpacing(rAngle)`: Calculates statistics related to the spacing between radial azimuthal angles.

- `RipleyRadial(rAngle, plot=False)`: Measures the "clumpiness" of radial dike swarm angles using the Ripley K function, with an option to visualize the results.

- `ExpandingR(lines, Center)`: Measures the density of lines at different radial distances from a center.

- `NearCenters(lines, Center, tol=10000, printOn=False)`: Identifies lines near a specified center within a tolerance, providing detailed information about those lines.

- `writeCenterWKT(df, name)`: Writes center information to a Well-Known Text (WKT) file, allowing for easy geospatial data export.
Overall, this module is designed for the analysis of geological dike structures, 
particularly in terms of their orientation and distribution relative to a central
 point. It includes functions for fitting radial models, assessing clustering patterns,
 and exporting results in a geospatial format.
 
 written by akh
'''
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
from scipy.spatial.distance import pdist, squareform
import matplotlib.pyplot as plt
from PrePostProcess import whichForm, writeToQGIS
from plotmod import plotlines, DotsLines
from clusterMod import CyclicAngleDist
from htMOD import HT_center

def CenterFunc(theta, xr, yr, xc, yc):
    """
    Calculate the radial distance of a point from a center given an angle.

    Parameters:
        theta (float): Angle in degrees.
        xr (float): X-coordinate of the point.
        yr (float): Y-coordinate of the point.
        xc (float): X-coordinate of the center.
        yc (float): Y-coordinate of the center.

    Returns:
        float: Radial distance of the point from the center.

    """
    rhoRadial = (xr - xc) * np.cos(np.deg2rad(theta)) + (yr - yc) * np.sin(np.deg2rad(theta))
    return rhoRadial

def RadialFit(lines, plot=False, ColorBy=None, weight='LayerNumber', ThetaRange=[-90, 90], xc=None, yc=None):
    """
    Fit a radial model to dike lines data.

    Parameters:
        lines (pd.DataFrame): Dataframe containing dike lines data.
        plot (bool, optional): Whether to plot the fitted model. Defaults to False.
        ColorBy (str, optional): Column name to color the lines in the plot. Defaults to None.
        weight (str, optional): Column name for weighting in the plot. Defaults to 'LayerNumber'.
        ThetaRange (list, optional): Range of angles to consider for the fit. Defaults to [-90, 90].
        xc (float, optional): X-coordinate of the center. Defaults to None.
        yc (float, optional): Y-coordinate of the center. Defaults to None.

    Returns:
        pd.DataFrame: Dataframe containing center information.

    """
    t, r = whichForm(lines)
    theta = lines[t].values
    rho = lines[r].values

    m = (theta > ThetaRange[0]) & (theta < ThetaRange[1])
    theta = theta[m]
    rho = rho[m]

    if plot:
        fig, ax = DotsLines(lines, ColorBy=ColorBy, cmap='turbo')
        xdata = np.linspace(-90, 90, 200)

    Centers = pd.DataFrame()

    if 'xc' in lines.columns and xc is None:
        xc = lines['xc'].values[0]
        yc = lines['yc'].values[0]
    elif xc is None:
        xc, yc = HT_center(lines)

    popt, pcov = curve_fit(lambda angle, xr, yr: CenterFunc(angle, xr, yr, xc, yc), theta, rho)
    perr = np.sqrt(np.diag(pcov))

    residuals = rho - CenterFunc(theta, *popt, xc, yc)
    ss_res = np.sum(residuals**2)
    ss_tot = np.sum((rho - np.mean(rho))**2)
    r_sq = 1 - (ss_res / ss_tot)
    Centers = pd.DataFrame({"Center": [popt], "Std Error": [perr], 'RSq': r_sq})

    if plot:
        ax[1].plot(xdata, CenterFunc(xdata, *popt, xc, yc), 'y-', label='fit: xr=%5.3f, yr=%5.3f' % tuple(popt), linewidth=3)
        ax[0].plot(popt[0], popt[1], '*g', markersize=10)
        plt.legend()

    return Centers

def RadialAzimuthal(lines, Center):
    """
    Calculate the radial azimuthal angles of dike lines relative to a given center.

    Parameters:
        lines (pd.DataFrame): Dataframe containing dike lines data.
        Center (pd.DataFrame): Dataframe containing center information.

    Returns:
        np.array: Array of radial azimuthal angles.

    """
    xdist = lines['Xmid'].values - Center['Center'][0][0]
    ydist = lines['Ymid'].values - Center['Center'][0][1]

    rAngle = np.rad2deg(np.arctan2(xdist, ydist)) + 180

    return rAngle

def CyclicAngle360(u, v):
    """
    Calculate the cyclic angle between two angles in the range [0, 360).

    Parameters:
        u (float): First angle.
        v (float): Second angle.

    Returns:
        float: Cyclic angle between the two angles.

    """
    return (u - v) % 360

def AngleSpacing(rAngle):
    """
    Calculate statistics related to angle spacing in a set of radial azimuthal angles.

    Parameters:
        rAngle (np.array): Array of radial azimuthal angles.

    Returns:
        tuple: Mean, median, minimum, and maximum angle spacing.

    """
    SrAngle = np.sort(rAngle)
    spacing = np.abs(SrAngle[0:-2] - SrAngle[1:-1])

    return np.mean(spacing), np.median(spacing), np.min(spacing), np.max(spacing)

def RipleyRadial(rAngle, plot=False):
    """
    Measure the radial dike swarm angle "clumpiness" using the Ripley K function.

    Parameters:
        rAngle (np.array): Array of radial azimuthal angles.
        plot (bool, optional): Whether to plot the Ripley K function. Defaults to False.

    Returns:
        tuple: Ripley K function, K function standard error, and optional plot figure and axis.

    """
    tRange = CyclicAngle360(np.max(rAngle), np.min(rAngle))
    theta = rAngle[:, None]
    steps = np.arange(0, 360, 10)

    n = len(rAngle)
    l = n / 360

    d = squareform(pdist(theta, metric=CyclicAngle360))

    counts = [np.histogram(i, bins=steps)[0] for i in d]
    K = l * np.cumsum(counts, axis=1) / n
    L = np.sqrt(np.sum(np.cumsum(counts, axis=1), axis=1) / (np.pi * n * (n - 1)))

    K_Pure = np.ones(len(K)) * l / n
    K_se = np.sum((K - K_Pure)**2) / n

    if plot:
        fg, ax = plt.subplots()
        ax.plot(steps[:-1], K_Pure, 'k.-')
        ax.plot(steps[:-1], K, 'r.-')
        return K, K_se, fg, ax
    else:
        return K, K_se

def ExpandingR(lines, Center):
    """
    Measure the density of lines at different radial distances from a center.

    Parameters:
        lines (pd.DataFrame): Dataframe containing dike lines data.
        Center (pd.DataFrame): Dataframe containing center information.

    Returns:
        list: List of density values at different radial distances.

    """
    t, r = whichForm(lines)
    tols = np.array([0.0001, 0.005, 0.01, 0.1, 0.2, 0.5, 0.75, 1]) * np.ptp(lines[r].values)
    ntol = np.empty(len(tols))

    xdata = lines[t].values
    xc = lines['xc'].values[0]
    yc = lines['yc'].values[0]

    rhoPerfect = CenterFunc(xdata, Center['Center'][0][0], Center['Center'][0][1], xc, yc)

    ntol = [np.sum(abs(lines[r].values - rhoPerfect) < tol) / len(lines) for tol in tols]

    return ntol

def NearCenters(lines, Center, tol=10000, printOn=False):
    """
    Find lines near a given center within a specified tolerance.

    Parameters:
        lines (pd.DataFrame): Dataframe containing dike lines data.
        Center (pd.DataFrame): Dataframe containing center information.
        tol (float, optional): Tolerance for considering lines as near centers. Defaults to 10000.
        printOn (bool, optional): Whether to print information about the lines. Defaults to False.

    Returns:
        pd.DataFrame: Dataframe containing lines near centers.

    """
    t, r = whichForm(lines)

    xdata = lines[t].values
    xc = lines['xc'].values[0]
    yc = lines['yc'].values[0]
    rhoPerfect = CenterFunc(xdata, Center['Center'][0][0], Center['Center'][0][1], xc, yc)

    close = abs(lines[r].values - rhoPerfect) < tol
    ntol = ExpandingR(lines, Center)
    Close = lines[close]
    rAngle = RadialAzimuthal(Close, Center)
    maxt = np.max(rAngle)
    mint = np.min(rAngle)
    spacing = AngleSpacing(rAngle)

    CenterDist1 = np.sqrt((Center['Center'][0][0] - Close['Xstart'].values)**2 + (Center['Center'][0][1] - Close['Ystart'].values)**2)
    CenterDist2 = np.sqrt((Center['Center'][0][0] - Close['Xend'].values)**2 + (Center['Center'][0][1] - Close['Yend'].values)**2)

    MinDist = np.min(np.vstack((CenterDist1, CenterDist2)), axis=0)

    Close = Close.assign(CenterDistStart=CenterDist1, CenterDistEnd=CenterDist2, CenterDist=MinDist, CenterX=Center['Center'][0][0], CenterY=Center['Center'][0][1], RadialAngles=rAngle)

    if "T" in t:
        clusters = len(Close)
        dikes = Close['Size'].sum()
        filt = Close['TrustFilter'].sum()

        Center = Center.assign(Spacing=[spacing], ExpandingR=[ntol], AngleRange=CyclicAngle360(maxt, mint), ClustersN=clusters, DikesN=dikes, FilteredN=filt)
    else:
        n = len(Close)
        Center = Center.assign(Spacing=[spacing], ExpandingR=[ntol], AngleRange=CyclicAngle360(maxt, mint), nDikes=n)

    if printOn:
        print("For Center", Center['Center'][0][0], Center['Center'][0][1])
        print("Max angle Range is ", CyclicAngle360(maxt, mint))
        print("Angle Spacing is ")
        print("Mean: Median: Min: Max:")
        print(spacing)
        print("Deviation from perfect radial angle spread")
        print("Perfect K_se would be 0")

    return Close, Center

def writeCenterWKT(df, name):
    """
    Write center information to a WKT file.

    Parameters:
        df (pd.DataFrame): Dataframe containing center information.
        name (str): File name to save the WKT data.

    Returns:
        pd.DataFrame: Dataframe containing center information.

    """
    front = "POINT ("
    linestring = []
    for i in range(len(df)):
        line = front + str((df['Center'].iloc[i])[0]) + " " + str((df['Center'].iloc[i])[1]) + ")"
        linestring.append(line)

    df['Pointstring'] = linestring
    df['Pointstring'] = df['Pointstring'].astype(str)
    df.to_csv(name)

    return df
