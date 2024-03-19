# LinkingLines Package
 # Written by aikubo
 
#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
FitRadialCenters Module

This module provides functions for fitting a radial model to dike lines data, calculating radial azimuthal angles, measuring the "clumpiness" of radial dike swarm angles, and more.

Functions:
- `CenterFunc(theta, xr, yr, xc, yc)`: Calculates the radial distance of a point from a given center for a specified angle.
- `RadialFit(lines, plot=False, ColorBy=None, weight='LayerNumber', ThetaRange=[-90, 90], xc=None, yc=None)`: Fits a radial model to dike lines data, optionally visualizing the fit. It returns center information such as coordinates and goodness of fit.
- `RadialAzimuthal(lines, Center)`: Calculates radial azimuthal angles of dike lines relative to a specified center.
- `CyclicAngle360(u, v)`: Computes the cyclic angle between two angles in the range [0, 360).
- `AngleSpacing(rAngle)`: Calculates statistics related to the spacing between radial azimuthal angles.
- `RipleyRadial(rAngle, plot=False)`: Measures the "clumpiness" of radial dike swarm angles using the Ripley K function, with an option to visualize the results.
- `ExpandingR(lines, Center)`: Measures the density of lines at different radial distances from a center.
- `NearCenters(lines, Center, tol=10000, printOn=False)`: Identifies lines near a specified center within a tolerance, providing detailed information about those lines.
- `writeCenterWKT(df, name)`: Writes center information to a Well-Known Text (WKT) file, allowing for easy geospatial data export.

'''
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
from scipy.spatial.distance import pdist, squareform
import matplotlib.pyplot as plt
from .PrePostProcess import whichForm
from .PlotUtils import DotsLines
import os
import geopandas as gpd


def CenterFunc(theta, xr, yr, xc, yc):
    """
    Calculate the radial distance of a point from a center given an angle.

    This function computes the radial distance of a point (or points) from a specified center based on
    a given angle. The radial distance, `rhoRadial`, is determined from the cartesian coordinates of the
    point and the center, along with the angle in degrees between them. This calculation is useful for
    geometric analyses and transformations in a polar coordinate system.

    Parameters
    ----------
    theta : float or array_like
        Angle(s) in degrees between the point(s) and the horizontal axis.
    xr : float
        X-coordinate of the point.
    yr : float
        Y-coordinate of the point.
    xc : float
        X-coordinate of the center.
    yc : float
        Y-coordinate of the center.

    Returns
    -------
    rhoRadial : float or ndarray
        The radial distance(s) of the point(s) from the center. If `theta` is a single float, `rhoRadial`
        will be a single float. If `theta` is an array_like, `rhoRadial` will be an ndarray of the same shape.

    """

    rhoRadial = (xr - xc) * np.cos(np.deg2rad(theta)) + (yr - yc) * np.sin(np.deg2rad(theta))
    return rhoRadial

def RadialFit(lines, plot=False, ColorBy=None, weight='LayerNumber', ThetaRange=[-90, 90], xc=None, yc=None):
    """
    Fit a radial model to dike lines data.

    This function fits a radial model to a dataset containing dike lines, optionally plots the fitted model,
    and returns a DataFrame with center information. The fitting process considers dike line orientations and
    positions to determine the most representative central point. Users can specify parameters for plotting,
    including whether to plot, how to color the lines, what weights to use, and the range of angles to consider
    for the fit. The function also allows specifying the center coordinates explicitly.

    Parameters
    ----------
    lines : pandas.DataFrame
        A DataFrame containing the dike lines data. Expected to include columns that match the `ColorBy`
        and `weight` parameters if they are used.
    plot : bool, optional
        If True, the fitted model and the lines will be plotted. Defaults to False.
    ColorBy : str, optional
        The name of the column in `lines` DataFrame to use for coloring the lines in the plot. If None,
        no coloring is applied. Defaults to None.
    weight : str, optional
        The name of the column in `lines` DataFrame to use for weighting in the plot. This can be used
        to emphasize certain lines over others. Defaults to 'LayerNumber'.
    ThetaRange : list of float, optional
        A list specifying the range of angles in degrees to consider for the fit, in the format [min, max].
        This can be used to limit the analysis to a certain orientation of dike lines. Defaults to [-90, 90].
    xc : float, optional
        The x-coordinate of the center to be used in the fit. If None, the function attempts to calculate
        the center automatically. Defaults to None.
    yc : float, optional
        The y-coordinate of the center to be used in the fit. If None, the function attempts to calculate
        the center automatically. Defaults to None.

    Returns
    -------
    Centers : pandas.DataFrame
        A DataFrame containing the information about the calculated center(s), including coordinates and
        possibly other metrics derived from the fitting process.

    Examples
    --------
    >>> import pandas as pd
    >>> import numpy as np
    >>> from linkinglines.FitRadialCenters import RadialFit
    >>> lines = pd.DataFrame({
    ...     'Theta': np.linspace(-90, 90, 100),
    ...     'Rho': np.linspace(0, 1, 100)
    ... })
    >>> Centers = RadialFit(lines, plot=True, ColorBy='LayerNumber', weight='LayerNumber',
    ...                     ThetaRange=[-90, 90], xc=None, yc=None)
    >>> print(Centers)

    """

    t, r = whichForm(lines)
    theta = lines[t].values
    rho = lines[r].values

    m = (theta > ThetaRange[0]) & (theta < ThetaRange[1])
    theta = theta[m]
    rho = rho[m]

    if plot:
        if ColorBy not in lines.columns:
            raise ValueError(f"Column '{ColorBy}' not found in the DataFrame.")
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

    In a radial structure of lines, the hough transform does not distinguish between lines on either side of the center. This function calculates the radial azimuthal angles of the lines relative to a given center.
    Which assigns lines based on their angle relative to the center.

    Parameters
    ----------
    lines : pandas.DataFrame
        Dataframe containing dike lines data.
    Center : pandas.DataFrame
        Dataframe containing center information.

    Returns
    -------
    rAngle : numpy.array
        Array of radial azimuthal angles.

    """
    xdist = lines['Xmid'].values - Center['Center'][0][0]
    ydist = lines['Ymid'].values - Center['Center'][0][1]

    rAngle = np.rad2deg(np.arctan2(xdist, ydist)) + 180

    return rAngle

def CyclicAngle360(u, v):
    """
    Calculate the cyclic angle between two angles in the range [0, 360).

    Parameters
    ----------
    u : float or array_like
        First angle.
    v : float or array_like
        Second angle.

    Returns
    -------
    d : float or array_like
        Cyclic angle between the two angles.

    """
    return (u - v) % 360

def AngleSpacing(rAngle):
    """
    Calculate statistics related to angle spacing in a set of radial azimuthal angles.

    This function computes statistical measures related to the spacing between consecutive angles in
    a set of radial azimuthal angles. It calculates the mean, median, minimum, and maximum spacing between
    angles. These statistics provide insights into the distribution and variability of angle spacings within
    the dataset, which can be critical for understanding patterns or uniformity in radial distributions.

    Parameters
    ----------
    rAngle : numpy.ndarray
        An array of radial azimuthal angles, assumed to be in degrees. The angles should be sorted in
        ascending order for accurate spacing calculations.

    Returns
    -------
    mean : float
        The mean (average) spacing between consecutive angles in the dataset.
    median : float
        The median spacing between consecutive angles, representing the middle value when the spacings
        are sorted in ascending order.
    min : float
        The minimum spacing between any two consecutive angles in the dataset, indicating the closest
        proximity between points.
    max : float
        The maximum spacing between any two consecutive angles, indicating the farthest apart points.

    See Also:
    --------
    nearCenters: Find lines near a given center within a specified tolerance.

    """
    SrAngle = np.sort(rAngle)
    spacing = np.abs(SrAngle[0:-2] - SrAngle[1:-1])

    return np.mean(spacing), np.median(spacing), np.min(spacing), np.max(spacing)

def RipleyRadial(rAngle, plot=False):
    """
    Measure the radial dike swarm angle "clumpiness" using the Ripley K function.

    Parameters
    ----------
    rAngle : numpy.array
        Array of radial azimuthal angles.
    plot : bool, default False
        Whether to plot the Ripley K function. Defaults to False.

    Returns
    -------
    K : numpy.array
        Ripley K function values.
    K_se : float
        Standard error of the Ripley K function.

        if plot:
            fg : matplotlib.figure.Figure
                The figure object.
            ax : matplotlib.axes.Axes
                The axes object.

    See Also:
    --------
    AngleSpacing: Calculate statistics related to angle spacing in a set of radial azimuthal angles.
    NearCenters: Find lines near a given center within a specified tolerance.

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

    Parameters
    ----------
    lines : pandas.DataFrame
        Dataframe containing dike lines data.
    Center : pandas.DataFrame
        Dataframe containing center information.

    Returns
    -------
    ntol : list
        List of density values at different radial distances.

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

    This function identifies and extracts lines from a dataset that are within a certain distance (tolerance)
    from a specified center point. The lines considered "near" are those whose distance from the center does
    not exceed the given tolerance. Optionally, the function can print details about these lines. The primary
    use of this function is in geospatial analyses, where identifying features relative to a point of interest
    is necessary, such as finding geological lines or faults near a geographic center.

    Parameters
    ----------
    lines : pandas.DataFrame
        A DataFrame containing the dike lines data. The DataFrame is expected to have columns that
        allow calculating the distance of each line from a specified center point.
    Center : pandas.DataFrame
        A DataFrame containing the center information. This DataFrame should include at least the
        x and y coordinates of the center.
    tol : float, optional
        The tolerance within which lines are considered near the center, in the same units as the
        line and center coordinates. Defaults to 10000 units.
    printOn : bool, optional
        If True, prints information about the lines found near the center. This can include distances,
        line IDs, or any other relevant information contained in the `lines` DataFrame. Defaults to False.

    Returns
    -------
    Close : pandas.DataFrame
        A DataFrame containing the subset of lines from the `lines` DataFrame that are within the
        specified tolerance of the center. The structure of this DataFrame mirrors that of `lines`.
    Center : pandas.DataFrame
        The unchanged DataFrame containing center information as provided in the input. This is returned
        to maintain consistency in function outputs and facilitate further processing if needed.

    See Also:
    --------
    CenterFunc: Calculate the radial distance of a point from a given center for a specified angle.
    RadialAzimuthal: Calculate radial azimuthal angles of dike lines relative to a specified center.
    ExpandingR: Measure the density of lines at different radial distances from a center.
    RipleyRadial: Measure the "clumpiness" of radial dike swarm angles using the Ripley K function.
    

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
    Write center information to a geospatial data file.
    Can be "csv" or "shp" or "geojson".

    Parameters
    ----------
    df : pandas.DataFrame 
        Dataframe containing center information.
    name : str
        File name to save the geospacial data.

    Returns
    -------
    df : pandas.DataFrame
        Dataframe containing center information.

    """
    front = "POINT ("
    linestring = []
    for i in range(len(df)):
        line = front + str((df['Center'].iloc[i])[0]) + " " + str((df['Center'].iloc[i])[1]) + ")"
        linestring.append(line)

    df['geometry'] = linestring
    df['geometry'] = df['geometry'].astype(str)

    if os.name.endswith('.csv'):
        df.to_csv(name)
    elif os.name.endswith('.shp') or os.name.endswith('.geojson'):
        df = gpd.GeoDataFrame(df, geometry='geometry')
        df.to_file(name)

    return df
