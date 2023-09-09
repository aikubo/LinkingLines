#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  1 12:54:02 2021

@author: akh

examineMod: Module for examining and analyzing line clusters.

This module provides functions to analyze and examine line clusters, including computing bounding rectangles, 
evaluating cluster properties, and checking for changes between clusters.

Functions:
    - OutputRectangles(clusters): Compute the coordinates of bounding rectangles for each cluster.
    - examineCluster(clusters): Generate a summary of properties for each cluster in a set of line clusters.
    - examineClustersShort(clusters): Generate a summary of properties for clusters in a set of line clusters.
    - checkClusterChange(lines1, lines2): Check if two sets of line clusters are the same.
    - checkoutCluster(clusters,label): Display information and plots related to a cluster of lines.
    - checkoutClusterCart(clusters, label): Display information and plots related to a cluster of lines.in Cartesian space
    - checkoutby(dikeset, lines, col): plot information based on cluster metrics
    - RotateOverlap(lines):Calculate the overlap ratio and maximum overlap count of lines after rotation.
    - enEchelonTwistAngle:Calculate the angle twist and statistical significance of features.
    -extendLines(lines):Extends lines in dataframe up to min(x)*.5 and max(x)*1.5
        
    
Dependences: 

    fitRectangle (this package)
    htMOD (this package)
    plotMod (this package)
    statsmodels
    scipy
    matplotlib
    numpy
    pandas
"""

import numpy as np 
import pandas as pd 
from fitRectangle import *
from htMOD import HT_center
from plotmod import HThist, plotlines, DotsHT, FixAxisAspect, clustered_lines, pltRec, FixCartesianLabels

from PrePostProcess import *
from htMOD import rotateData, CyclicAngleDist
from clusterMod import *
import statsmodels.api as sm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap

from PrePostProcess import completePreProcess, whichForm, midPoint


import matplotlib.pyplot as plt 
from scipy import stats
import warnings
from matplotlib import cm
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

from datetime import datetime


def checkoutCluster(dikeset, label):
    """
    Display information and plots related to a cluster of lines.

    Parameters:
        dikeset (DataFrame): A DataFrame containing line data with columns 'xc', 'yc', 'Labels', 'theta', and 'rho'.
        label (int): The label of the cluster to be analyzed.

    Returns:
        fig, ax (tuple): A tuple containing the Figure and Axes objects for the generated plots.

    This function takes a DataFrame containing line data and a label specifying a cluster of lines. It generates two subplots:
    1. Scatter plot of lines' theta and rho values, highlighting the selected cluster in red.
    2. Rectangle plot showing the lines in the cluster along with additional information.

    The function calculates and displays the following information:
    - Mean angle (in degrees) of the lines in the cluster.
    - Mean rho (in km) of the lines in the cluster.
    - Length and width of the cluster.
    - Aspect ratio (length/width) of the cluster.
    - Size of the cluster (number of lines).

    Additionally, it handles cases where the cluster crosses the zero angle boundary and adjusts the plot accordingly.

    Note:
    - The input DataFrame 'dikeset' must contain columns 'xc', 'yc', 'Labels', 'theta', and 'rho'.
    - 'label' specifies the cluster to be analyzed.
    - The function returns the Figure and Axes objects for further customization or saving.

    Example usage:
    fig, ax = checkoutCluster(lines, 2)
    plt.show()
    """

    xc=dikeset['xc'].iloc[0]
    yc=dikeset['yc'].iloc[0]
    
    fig, ax=plt.subplots(1,2) 
    fig.set_size_inches( 12,6)
    fig.suptitle("Label "+str(int(label)))
    mask=(dikeset['Labels']==label)
    lines=dikeset[mask]
    size=len(lines)
    
    if abs(np.sum(np.sign(lines['theta'].values))) < size: 


        ang=np.mean(abs(lines['theta'].values))

        if np.isclose(ang,0, atol=4):
            ang=np.mean((lines['theta'].values))
    else:
        ang=np.average((lines['theta']))
        
        
    print(ang)
    xlim1=min(lines['theta'])-2
    xlim2=max(lines['theta'])+2
    ylim1=(min(lines['rho'])-10000)/1000
    ylim2=(max(lines['rho'])+10000)/1000
    
    ax[0].scatter(dikeset['theta'], dikeset['rho']/1000, c='grey', alpha=0.2) 
    ax[0].scatter(lines['theta'], lines['rho']/1000, c='r')
    ax[0].plot(ang, lines['rho'].mean()/1000, 'g*')
    
    ax[0].plot([xlim1, xlim1, xlim2, xlim2, xlim1], [ylim1,ylim2, ylim2, ylim1, ylim1], 'k')
    
    if np.mean(lines['rho'].values) > 0 :
        left, bottom, width, height = [0.2, 0.2, 0.2, 0.2]
    else:
        left, bottom, width, height = [0.2, 0.55, 0.2, 0.2]
        

    ax2 = fig.add_axes([left, bottom, width, height])
    ax2.scatter(lines['theta'], lines['rho']/1000, c='r')
    ax2.plot(ang, lines['rho'].mean()/1000, 'g*')
    ax2.scatter(dikeset['theta'], dikeset['rho']/1000, c='grey', alpha=0.1) 
    ax2.set_ylim([ylim1+9.5, ylim2-9.5])
    ax2.set_xlim([xlim1+1.5, xlim2-1.5])
    
    ax[0].set_xlabel('Theta ($^\circ$)')
    ax[0].set_ylabel('Rho (km)')
    
    fig,ax[1], l,w=pltRec(lines, xc,yc, fig=fig, ax=ax[1])

    
    if np.sign(ang) < 0:
            loc=4
            tlocx=.05
            tlocx2=.60
            
    else: 
            
            loc=3
            tlocx=.60
            tlocx2=.05
        
    
    ax4 = inset_axes(ax[1], width="30%", height="30%", loc=loc)

    plotlines(dikeset, 'grey', ax4, alpha=0.1)
    plotlines(lines, 'r', ax4)
    ax4.xaxis.tick_top()
    ax4.xaxis.set_label_position('top') 
    
    if ang> 0:
        ax4.yaxis.tick_right()
        ax4.yaxis.set_label_position('right') 
        

    FixCartesianLabels(ax4)
    FixCartesianLabels(ax[1])
    FixAxisAspect(ax[0],ax[1])
    
    if np.sign(ang)>0:
        tlocy=0.95
    else:
        tlocy=0.25
        

    t=ax[0].text(tlocx2,tlocy-0.10,'Angle Mean $(^\circ)$: '+str( round(ang, 2)), transform=ax[0].transAxes)

    t=ax[0].text(tlocx2,tlocy-0.05,'Rho Mean (km): '+str(round(lines['rho'].mean()/1000,1)), transform=ax[0].transAxes)

        
    t=ax[0].text(tlocx2,tlocy-0.20,'Angle Range $(^\circ)$: '+str( round(np.ptp(lines['theta']), 2)), transform=ax[0].transAxes)

    t=ax[0].text(tlocx2,tlocy-0.15,'Rho Range (km): '+str(round(np.ptp(lines['rho'])/1000,1)), transform=ax[0].transAxes)


    tlocy=.95
    ax[1].text(tlocx,tlocy,'Length: '+str(round(l,1)), transform=ax[1].transAxes)
    ax[1].text(tlocx,tlocy-0.05,'Width: '+str(round(w,1)), transform=ax[1].transAxes)
    ax[1].text(tlocx,tlocy-0.10,'Aspect: '+str(round(l/w,1)), transform=ax[1].transAxes)
    ax[1].text(tlocx,tlocy-0.15,'Size: '+str(len(lines)), transform=ax[1].transAxes)
    #plt.tight_layout()
    
    print('Angle Mean: '+str( round(ang, 2)))
    print('Rho Mean: '+str(round(lines['rho'].mean(),1)))
    print('Length: '+str(round(l,1)))
    print('Width: '+str(round(w,1)))
    
    


    return fig,ax

def checkoutClusterCart(dikeset, label, fig, ax):
    """
    Visualize and annotate cluster information in a Cartesian coordinate system.

    Parameters:
        dikeset (DataFrame): A DataFrame containing line data with columns 'xc', 'yc', 'Labels', 'theta', and 'rho'.
        label (int): The label of the cluster to be analyzed.
        fig (Figure): The Figure object for plotting.
        ax (Axes): The Axes object for plotting.

    Returns:
        ax (Axes): The updated Axes object after plotting and annotation.

    This function visualizes cluster information in a Cartesian coordinate system. It takes a DataFrame containing line data,
    a cluster label, and Figure and Axes objects for plotting. The function generates a plot of the cluster's line segments,
    annotates the plot with information such as the mean angle, length, width, and size of the cluster, and adjusts the plot
    based on the cluster's orientation.

    The function also prints the cluster's angle mean, rho mean, length, and width for reference.

    Example usage:
    fig, ax = plt.subplots()
    ax = checkoutClusterCart(line_data, 2, fig, ax)
    plt.show()
    """

    xc=dikeset['xc'].iloc[0]
    yc=dikeset['yc'].iloc[0]
    
    mask=(dikeset['Labels']==label)
    lines=dikeset[mask]
    size=len(lines)
    
    if abs(np.sum(np.sign(lines['theta'].values))) < size: 
        crossZero=True

        ang=np.mean(abs(lines['theta'].values))
        tol=6
        if np.isclose(ang,0, atol=4):
            ang=np.mean((lines['theta'].values))
    else:
        crossZero=False
        ang=np.average((lines['theta']))
        
        
    print(ang)
    xlim1=min(lines['theta'])-2
    xlim2=max(lines['theta'])+2
    ylim1=(min(lines['rho'])-10000)/1000
    ylim2=(max(lines['rho'])+10000)/1000
    
    if np.mean(lines['rho'].values) > 0 :
        left, bottom, width, height = [0.2, 0.2, 0.2, 0.2]
        loc='lower left'
    else:
        left, bottom, width, height = [0.2, 0.55, 0.2, 0.2]
        loc='upper right'


    fig,ax, l,w=pltRec(lines, xc,yc, fig=fig, ax=ax)

    
    if np.sign(ang) < 0:
            loc=4
            tlocx=.05
            tlocx2=.60
            
    else: 
            
            loc=3
            tlocx=.60
            tlocx2=.05
    
    FixCartesianLabels(ax)

    
    if np.sign(ang)>0:
        tlocy=0.95
    else:
        tlocy=0.25
        

    t=ax.text(tlocx2,tlocy,'Angle Mean $(^\circ)$: '+str( round(ang, 2)), transform=ax.transAxes)
    ax.text(tlocx2,tlocy-0.05,'Length (km): '+str(round(l/1000,1)), transform=ax.transAxes)
    ax.text(tlocx2,tlocy-0.10,'Width (m): '+str(round(w,1)), transform=ax.transAxes)
    ax.text(tlocx2,tlocy-0.15,'Size: '+str(len(lines)), transform=ax.transAxes)
    #plt.tight_layout()
    
    print('Angle Mean: '+str( round(ang, 2)))
    print('Rho Mean: '+str(round(lines['rho'].mean(),1)))
    print('Length: '+str(round(l,1)))
    print('Width: '+str(round(w,1)))
    return ax

def CheckoutBy(dikeset, lines, col, maximum=True, minimum=False):
    """
    Display cluster information and plots for lines with maximum or minimum values in a specified column.

    Parameters:
        dikeset (DataFrame): A DataFrame containing line data with columns 'xc', 'yc',
        'Labels', 'theta', and 'rho'.
        lines (DataFrame): A DataFrame containing line data with a 'Label' column and 
        the specified 'col' column.
        col (str): The column name in 'lines' to identify the maximum or minimum values.
        maximum (bool): If True, display information and plots for lines with the maximum 'col' value.
        minimum (bool): If True, display information and plots for lines with the minimum 'col' value.

    Returns:
        fig, ax (tuple): A tuple containing the Figure and Axes objects for the generated plots.

    Notes:
    This function allows you to analyze and visualize clusters of lines based 
    on the maximum or minimum values in a specified column. It can display 
    cluster information and plots for both maximum and minimum values independently.

    Parameters 'dikeset' and 'lines' should contain the necessary line data,
    and 'col' should specify the column in 'lines' to determine maximum or 
    minimum values. You can choose to display information for maximum values,
    minimum values, or both by setting the 'maximum' and 'minimum' flags accordingly.

    The function utilizes the 'checkoutCluster' function to generate cluster 
    information and plots and prints the label and the corresponding maximum or 
    minimum value of the specified column.

    Example usage:
    fig, ax = CheckoutBy(lines, dikeset, 'theta', maximum=True, minimum=True)
    plt.show()
    """


    if maximum:
        loc=np.nanargmax(lines[col].values)
        label=lines['Label'].iloc[loc]
        print("Label", int(label), "maximum", col)
        
        fig,ax=checkoutCluster(dikeset, label)
        
    if minimum:
        loc=np.nanargmin(lines[col].values)
        label=lines['Label'].iloc[loc]
        print("Label", int(label), "minimum", col)
        fig,ax=checkoutCluster(dikeset, label)
        
    return fig,ax

def RotateOverlap(lines):
    """
    Calculate the overlap ratio and maximum overlap count of lines after rotation.

    Parameters:
        lines (DataFrame): A DataFrame containing line data with 'theta', 'Xstart', and 'Xend' columns.

    Returns:
        overlap (float): The ratio of overlapping line segments after rotation.
        max_overlap_count (int): The maximum count of overlapping line segments along the rotated x-axis.

    This function takes a DataFrame of line data and performs the following steps:
    1. Calculates the mean angle ('theta') of the lines.
    2. Rotates the line data by the complementary angle (90 degrees minus the mean angle).
    3. Transforms the 'Xstart' values to ensure they are in the correct order.
    4. Computes the total length of the line segments.
    5. Determines the overlapping segments along the rotated x-axis and calculates the overlap ratio.

    The overlap ratio is defined as the ratio of the length of overlapping 
    segments to the total length of the segments. Higher ratios are more overlap
    
    A ratio of 1 would indicate it has two segments that are totally overlapping 
    while a ratio of 2 would indicate 3 segments totally overlapping.
    
    The maximum overlap count represents the maximum number of overlapping 
    segments along the rotated x-axis.

    Example usage:
    overlap, max_overlap_count = RotateOverlap(line_data)
    print("Overlap Ratio:", overlap)
    print("Maximum Overlap Count:", max_overlap_count)
    """


    theta=np.mean(lines['theta'].values)
    
    dfRotated=rotateData(lines, (90-theta))
    dfRotated=transformXstart(dfRotated)
    
    
    Xstart=dfRotated['Xstart'].to_numpy()
    Ystart=dfRotated['Ystart'].to_numpy()
    Xend=dfRotated['Xend'].to_numpy()
    Yend=dfRotated['Yend'].to_numpy()
    step=1
    totalL=np.sum( np.sqrt( (Xstart-Xend)**2 + (Ystart-Yend)**2))
    xs=[np.arange(min(x,y), max(x,y), step) for x,y in zip(np.floor(Xstart),np.ceil(Xend))]
    
    l=np.max([len(xi) for xi in xs])
    xs_sameLength=[np.append(xi,[np.nan]*(l-len(xi))) for xi in xs]

    arr=np.vstack(xs_sameLength) #better
    u,xcounts=np.unique(arr[~np.isnan(arr)], return_counts=True)
    
    overlapx=np.float64(np.sum(xcounts[xcounts>1]-1)*step)
    

    overlap=overlapx/(totalL-overlapx+0.00001)
    if overlapx > totalL:
        overlap=overlapx/totalL
        
    return overlap, np.max(xcounts)

def enEchelonAngleTwist(d, avgtheta):
    """
    Calculate the angle twist and statistical significance for en échelon features.

    Parameters:
        d (DataFrame): A DataFrame containing en échelon line segment data with 'Xstart', 'Xend', 'Ystart', 'Yend', 
                       'Xmid', and 'Ymid' columns.
        avgtheta (float): The average angle (in degrees) for alignment comparison.

    Returns:
        angle_twist (float): The angle twist between the average angle and the en échelon feature's orientation.
        p_value (float): The p-value indicating the statistical significance of the alignment.
        
        Xstart (float): The starting x-coordinate of the midpoints line.
        Xend (float): The ending x-coordinate of the midpoints line.
        Ystart (float): The corresponding y-coordinate at the starting x-coordinate.
        Yend (float): The corresponding y-coordinate at the ending x-coordinate.

    This function analyzes features represented as line segments and performs the following steps:
    1. Fits a linear regression model to the midpoints of the  line feature to estimate its orientation.
    3. Computes the p-value for the linear regression to assess the statistical significance of the alignment.
    4. If the p-value is above 0.05, the angle twist is set to 0 (indicating no significant alignment).
    5. If the p-value is below 0.05 (inidicating significant allignment of midpoints)
    it caluclates the angle difference between the line created by the midpoints
    and the average angle of the lines

    Example usage:
    angle_twist, p_value, Xstart, Xend, Ystart, Yend = enEchelonAngleTwist(en_echelon_data, 30.0)
    print("Angle Twist (degrees):", angle_twist)
    print("P-Value:", p_value)
    print("Starting Coordinates:", Xstart, Ystart)
    print("Ending Coordinates:", Xend, Yend)
    """


    warnings.filterwarnings("ignore")
    if len(d) <3: 
        return 0, 1, d["Xstart"], d["Xend"], d["Ystart"],d["Yend"]
    
    
    Xmid=d['Xmid'].to_numpy()
    Ymid=d['Ymid'].to_numpy()
    
    slope, intercept, r_value, p_value_bad, std_err = stats.linregress(Xmid, Ymid)
    thetaM=np.rad2deg(np.arctan(-1/(slope+ 0.0000000000000001)))
    # we don't actually want p from linregress... 
    Xstart=min(Xmid)
    Xend=max(Xmid)
    Ystart=slope*Xstart+intercept 
    Yend=slope*Xend+intercept 
    tdiff=CyclicAngleDist([avgtheta], [thetaM])
    #print(tdiff, thetaM, p_value)
    mod = sm.OLS(Ymid,Xmid)
    fii = mod.fit()
    
    p_values = fii.summary2().tables[1]['P>|t|']
    
    if p_values.values[0] > 0.05:
        tdff=0
    else: 
        tdiff=CyclicAngleDist([avgtheta], [thetaM])
        
    if np.isnan(tdiff):
        tdiff=0
    
    
    return tdiff, p_values.values[0], Xstart, Xend, Ystart, Yend

def examineClusterShort(clusters):
    """
    Analyze and summarize information about clustered line segments.

    Parameters:
        clusters (DataFrame): A DataFrame containing line data with columns 'Labels', 'theta', 'rho', and 'HashID'.

    Returns:
        clusters_data (DataFrame): A DataFrame containing summarized information for each cluster.

    This function analyzes and summarizes information about clustered line segments. 
    It calculates various statistics for each cluster,
    including the starting and ending coordinates, the average rho (distance from the origin), 
    the average theta (angle), and the size
    (number of lines) within each cluster. Additionally, it computes the midpoint of each 
    cluster and determines if it is a clustered or
    non-clustered line segment.

    The function returns a DataFrame ('clusters_data') containing the summarized 
    information for each cluster, including its label,
    coordinates, average rho, average theta, size, and a hash identifier.
    
    Notes:
        This is the shorter and faster version of examineCluster

    Example usage:
    cluster_summary = examineClusterShort(cluster_data)
    print(cluster_summary.head())
    """
    # Function code here

    clabel=np.unique(clusters['Labels'])
    nclusters=len(clabel)-1
    notclustered=sum([clusters['Labels']==-1][0])
    xc,yc=HT_center(clusters)
    clusters_data=pd.DataFrame()
    ids=np.arange(0, len(clusters),1)
    if "HashID" not in clusters.columns:
        clusters=giveHashID(clusters)

    #print("Found", nclusters)
    #print( sum( clusters > -1), "clustered out of", len(p))
    sizes=[]
    EEDikes=pd.DataFrame()
    
    cmask=np.full(len(clusters),False)
    for i in np.unique(clusters['Labels']): 
        clustered=True
        mask=clusters['Labels']==i
        lines=clusters[mask]
        cmask[mask]=True
        if (i == -1 or len(lines)<2):
            clustered=False
            #continue
            
        x,y=endpoints2(lines)

        x0=(max(x)-min(x))/2
        y0=(max(y)-min(y))/2
        size=np.sum(mask)
        avgrho=np.average(lines['rho'])
        
        if abs(np.sum(np.sign(lines['theta'].values))) < size: 
            crossZero=True
    
            avgtheta=np.mean(abs(lines['theta'].values))
            if np.isclose(avgtheta,0, atol=6):
                avgtheta=np.mean((lines['theta'].values))
        else:
            crossZero=False
            avgtheta=np.average((lines['theta']))
        w,l,r,Xe, Ye, Xmid, Ymid=fit_Rec(lines, xc, yc)

        x1, x2, y1, y2=clustered_lines(x,y,avgtheta, l, xmid=Xmid, ymid=Ymid)

        
        hashlines=hash((lines['HashID'].values.tostring()))
        clusters_data=clusters_data.append({ "Label": i, "Xstart": x1, "Ystart":y1, "Xend": x2,
                                            "Yend":y2, "X0": x0, "Y0": y0, "AvgRho":avgrho,
                                            "AvgTheta":avgtheta, "ClusterHash": hashlines, "Size":size}, ignore_index=True)
    return clusters_data
                                       
def examineClusters(clusters, enEchelonCutofff=7, ifEE=False, MaxNNSegDist=0.5, skipUnlinked=True, xc=None, yc=None):
    """
    Analyze and summarize information about clusters of line segments.

    Parameters:
        clusters (DataFrame): A DataFrame containing line data with columns 'Xstart', 'Ystart', 'Xend', 'Yend', 'seg_length',
                              'ID', 'rho', 'theta', 'Labels', and 'PerpOffsetDist'.
        enEchelonCutofff (int): A cutoff value for en échelon angle differences.
        ifEE (bool): A flag indicating whether en échelon analysis should be performed.
        MaxNNSegDist (float): The maximum normalized nearest neighbor segment distance for the TrustFilter.
        skipUnlinked (bool): A flag to skip unlinked segments.
        xc (float): The x-coordinate of the center (optional).
        yc (float): The y-coordinate of the center (optional).

    Returns:
        clusters_data (DataFrame): A DataFrame containing summarized information for each cluster.
        
            1. `Label`: Cluster label or identifier.
            2. `Xstart`: Starting x-coordinate of the clustered line.
            3. `Ystart`: Starting y-coordinate of the clustered line.
            4. `Xend`: Ending x-coordinate of the clustered line.
            5. `Yend`: Ending y-coordinate of the clustered line.
            6. `X0`: Midpoint of the x-coordinate range of the cluster.
            7. `Y0`: Midpoint of the y-coordinate range of the cluster.
            8. `AvgRho`: Average rho for the cluster.
            9. `AvgTheta`: Average angle (theta) for the cluster.
            10. `AvgSlope`: Average slope of the lines in the cluster.
            11. `AvgIntercept`: Average intercept of the lines in the cluster.
            12. `RhoRange`: Range of rho values within the cluster.
            13. `Aspect`: Aspect ratio, calculated as the length (l) divided by the width (w).
            14. `Xmid`: X-coordinate of the midpoint of the fitted rectangle
            15. `Ymid`: Y-coordinate of the midpoint of the fitted rectangle
            16. `PerpOffsetDist`: Average perpendicular offset distance for the lines in the cluster.
            17. `PerpOffsetDistRange`: Range of perpendicular offset distances within the cluster.
            18. `NormPerpOffsetDist`: Normalized perpendicular offset distance.
            19. `ThetaRange`: Range of theta values within the cluster.
            20. `StdRho`: Standard deviation of rho values within the cluster.
            21. `StdTheta`: Standard deviation of theta values within the cluster.
            22. `R_Width`: Width (w) of the cluster.
            23. `R_Length`: Length (l) of the cluster.
            24. `Size`: Number of lines in the cluster.
            25. `R_error`: Square root of the error (r) in the cluster's line fit.
            26. `Linked`: Indicates whether the lines in the cluster are considered linked or not.
            27. `SegmentLSum`: Sum of the lengths of line segments within the cluster.
            28. `ClusterHash`: A hash identifier for the cluster.
            29. `ClusterCrossesZero`: Indicates whether the cluster's angles cross zero.
            30. `EnEchelonAngleDiff`: Twist angle difference for features within the cluster.
            31. `Overlap`: Overlap of line segments within the cluster.
            32. `nOverlapingSegments`: Number of overlapping segments within the cluster.
            33. `EEPvalue`: P-value related to en échelon analysis.
            34. `MaxSegNNDist`: Maximum normalized nearest neighbor segment distance.
            35. `MedianSegNNDist`: Median normalized nearest neighbor segment distance.
            36. `MinSegNNDist`: Minimum normalized nearest neighbor segment distance.
            37. `TrustFilter`: A filter indicating trustworthiness based on the
            maximum normalized nearest neighbor segment distance.
            38. 'xc': X-coordinate of HT origin
            39. 'yc': Y-coordinate of HT origin
            40: 'Date_Changed': date string of generation or change time

        evaluation (DataFrame): A DataFrame containing summary statistics of the clusters.
        
        1. `nClusters`: The number of clusters in the `clusters_data` DataFrame.
        2. `nDikePackets`: The number of clusters with an overlap greater than 0.1 (presumably indicating some form of overlap between line segments).
        3. `AverageRhoRange`: The average range of rho values within the clusters.
        4. `MaxRhoRange`: The maximum range of rho values within the clusters.
        5. `StdRhoRange`: The standard deviation of the range of rho values within the clusters.
        6. `AverageThetaRange`: The average range of theta values within the clusters.
        7. `MaxThetaRange`: The maximum range of theta values within the clusters.
        8. `StdThetaRange`: The standard deviation of the range of theta values within the clusters.
        9. `AvgClusterSize`: The average size (number of lines) of the clusters.
        10. `ClusterSizeStd`: The standard deviation of the size of the clusters.
        11. `ClusterMax`: The maximum size (number of lines) among the clusters.
        12. `AverageL`: The average length (l) of the clusters.
        13. `MaxL`: The maximum length (l) among the clusters.
        14. `StdL`: The standard deviation of the length (l) of the clusters.
        15. `AverageW`: The average width (w) of the clusters.
        16. `MaxW`: The maximum width (w) among the clusters.
        17. `StdW`: The standard deviation of the width (w) of the clusters.
        18. `nTrustedDikes`: The number of clusters that pass a trust filter (presumably based on some criteria).
        19. `MaxEEAngleDiff`: The maximum en échelon angle difference among the clusters.
        20. `AverageEAngleDiff`: The average en échelon angle difference among the clusters.
        21. `Date`: The date when this summary information was generated.

        

    This function analyzes and summarizes information about clusters of line segments. 
    It calculates various statistics for each cluster,
    including coordinates, average rho (distance from the origin), average theta 
    (angle), size (number of lines) within each cluster, and
    other cluster-related information. It also includes a TrustFilter based o
    n the maximum normalized nearest neighbor segment distance.

    The function returns two DataFrames: 'clusters_data' contains summarized 
    information for each cluster, and 'evaluation' contains summary
    statistics of the clusters.

    Example usage:
    cluster_summary, cluster_evaluation = examineClusters(cluster_data,
                                                          enEchelonCutofff=10, 
                                                          ifEE=True, 
                                                          MaxNNSegDist=0.6, 
                                                          skipUnlinked=True)

    """

    clabel=np.unique(clusters['Labels'])
    if xc is None:
        xc,yc=HT_center(clusters)
        
        
    clusters_data=pd.DataFrame()
    ids=np.arange(0, len(clusters),1)
    if "HashID" not in clusters.columns:
        clusters=giveHashID(clusters)


    sizes=[]
    EEDikes=pd.DataFrame()
    
    cmask=np.full(len(clusters),False)
    perpoffsetMean=clusters['PerpOffsetDist'].mean()
    perpoffsetCutoff=clusters['PerpOffsetDist'].mean()+clusters['PerpOffsetDist'].std()
    pofdCutoff=perpoffsetCutoff/perpoffsetMean
    for i in np.unique(clusters['Labels']): 
        clustered=True
        mask=clusters['Labels']==i
        lines=clusters[mask]
        cmask[mask]=True
        if (i == -1 or len(lines)<2):
            clustered=False
            if skipUnlinked:
                continue
        
        x,y=endpoints2(lines)
        avgrho=np.average(lines['rho'])
        x0=(max(x)-min(x))/2
        y0=(max(y)-min(y))/2
        size=np.sum(mask)
        

        if abs(np.sum(np.sign(lines['theta'].values))) < size: 
            crossZero=True
    
            avgtheta=np.mean(abs(lines['theta'].values))
            tol=6
            if np.isclose(avgtheta,0, atol=4):
                avgtheta=np.mean((lines['theta'].values))
        else:
            crossZero=False
            avgtheta=np.average((lines['theta']))

        rrange=max(lines['rho'])-min(lines['rho'])
        trange=CyclicAngleDist([lines['theta'].min()], [lines['theta'].max()])

        
        
        stdrho=np.std(lines['rho'])
        stdt=np.std(lines['theta'])
        segmentL=lines['seg_length'].sum()
        

        w,l,r,Xe, Ye, Xmid, Ymid=fit_Rec(lines, xc, yc)

        x1, x2, y1, y2=clustered_lines(x,y,avgtheta, l, xmid=Xmid, ymid=Ymid)

        
        hashlines=hash((lines['HashID'].values.tostring()))

        if (l-w > 0): 
            EE=enEchelonAngleTwist(lines, avgtheta)
            tdiff=EE[0]
            EE_pvalue=EE[1]
        else: 
            tdiff=0
            EE_pvalue=1
        slope=-1/(np.tan(avgtheta)+0.00000000001)
        step=lines['seg_length'].min()+1
        
        npofd=lines['PerpOffsetDist'].mean()/perpoffsetMean
        b=avgrho/np.sin((np.deg2rad(avgtheta))+0.00000000001)+yc-xc*slope
        
        if size>1: 
            overlap, nOverlap=RotateOverlap(lines)
        else:
            overlap=0
            nOverlap=0
            
        if size>3: 
            X=(np.vstack((lines['Xmid'].values, lines['Ymid'].values)).T)
            NN= NearestNeighbors(n_neighbors=1, algorithm='auto', metric='euclidean').fit(X)
            distances, nbrs = NN.kneighbors()            
            MaxDist=np.max(distances)
            MinDist=np.min(distances)
            MedianDist=np.median(distances)
        else:
            MaxDist=np.nan
            MinDist=np.nan
            MedianDist=np.nan
            
        if MaxDist/l < 0.5 and npofd>1:
            Trust=True
        else:
            Trust=False
            
        clusters_data=clusters_data.append({ "Label": i, "Xstart": x1, "Ystart":y1, "Xend": x2,
                                            "Yend":y2, "X0": x0, "Y0": y0, "AvgRho":avgrho,
                                            "AvgTheta":avgtheta, "AvgSlope": slope, "AvgIntercept": b ,
                                            "RhoRange":rrange, "Aspect": l/w, 'Xmid': Xmid, 'Ymid': Ymid,
                                            "PerpOffsetDist": lines['PerpOffsetDist'].mean(),
                                            "PerpOffsetDistRange": lines['PerpOffsetDist'].max()-lines['PerpOffsetDist'].min(),
                                            "NormPerpOffsetDist":npofd,
                                            "ThetaRange": trange, "StdRho": stdrho, 
                                            "StdTheta": stdt, "R_Width": w, "R_Length": l, 
                                            "Size":size, "R_error":np.sqrt(r), "Linked":clustered,
                                            "SegmentLSum":segmentL, "ClusterHash":hashlines, 
                                            'ClusterCrossesZero': crossZero,
                                            "EnEchelonAngleDiff":tdiff, "Overlap": overlap,
                                            'nOverlapingSegments':nOverlap, 
                                            "EEPvalue":EE_pvalue, "MaxSegNNDist": MaxDist/l,
                                            "MedianSegNNDist": MedianDist/l, "MinSegNNDist": MinDist/l, "TrustFilter": Trust},
                                            ignore_index=True)
    
    clusters_data.astype({'Linked':'bool'})
   
    if skipUnlinked:
        nclusters= len(clusters_data)
    else:
        nclusters=np.sum(clusters_data['Size']>1)
        

    now = datetime.now()
    date = now.strftime("%d %b, %Y")
    clusters_data=clusters_data.assign(Date_Changed=date, xc=xc, yc=yc)
    evaluation=pd.DataFrame({ "nClusters": nclusters,
                              "nDikePackets": np.sum(clusters_data['Overlap']>0.1),
                              "AverageRhoRange": np.average(clusters_data["RhoRange"]),
                              "MaxRhoRange": np.max(clusters_data["RhoRange"]),
                              "StdRhoRange": np.std(clusters_data["RhoRange"]),                              
                              "AverageThetaRange": np.average(clusters_data["ThetaRange"]),
                              "MaxThetaRange": np.max(clusters_data["ThetaRange"]), 
                              "StdThetaRange": np.std(clusters_data["ThetaRange"]),
                              "AvgClusterSize": np.average(clusters_data["Size"]),
                              "ClusterSizeStd": np.std(clusters_data["Size"]),
                              "ClusterMax": clusters_data["Size"].max(), 
                              "AverageL": np.average(clusters_data["R_Length"]),
                              "MaxL": np.max(clusters_data['R_Length']),
                              "StdL": np.std(clusters_data["R_Length"]),
                              "AverageW": np.average(clusters_data["R_Width"]),
                              "MaxW": np.max(clusters_data['R_Width']),
                              "StdW": np.std(clusters_data['R_Width']),
                              "nTrustedDikes":np.sum(clusters_data['TrustFilter']>0),
                              "MaxEEAngleDiff":np.max(clusters_data["EnEchelonAngleDiff"]),
                              "AverageEAngleDiff":np.average(clusters_data["EnEchelonAngleDiff"]),
                              "Date":date},
                              index=[0])
    

    return clusters_data, evaluation

def evaluationOnClusters(clusters_data):
    """
    Calculate evaluation metrics based on cluster data.

    Args:
    clusters_data (DataFrame): A DataFrame containing cluster information.

    Returns:
    DataFrame: A summary DataFrame with evaluation metrics.

    This function computes various evaluation metrics based on the provided cluster data. 
    The metrics include information about the number of clusters, cluster sizes, 
    rho and theta range statistics, average lengths and widths of clusters, 
    en échelon angle differences, and more. 

    The resulting summary DataFrame provides insights into the distribution and 
    characteristics of the clusters, making it useful for further analysis and 
    interpretation of the data.

    Example:
    >>> evaluation = evaluationOnClusters(clusters_data)
    """


    now = datetime.now()
    date = now.strftime("%d %b, %Y")
    nclusters=np.sum(clusters_data['Size']>1)
    evaluation=pd.DataFrame({ "nClusters": nclusters,
                          "nDikePackets": np.sum(clusters_data['Overlap']>0.2),
                          "nTrustedDikes":np.sum(clusters_data['TrustFilter']>0),
                          "AverageRhoRange": np.average(clusters_data["RhoRange"]),
                          "MaxRhoRange": np.max(clusters_data["RhoRange"]),
                          "StdRhoRange": np.std(clusters_data["RhoRange"]),                              
                          "AverageThetaRange": np.average(clusters_data["ThetaRange"]),
                          "MaxThetaRange": np.max(clusters_data["ThetaRange"]), 
                          "StdThetaRange": np.std(clusters_data["ThetaRange"]),
                          "AvgClusterSize": np.average(clusters_data["Size"]),
                          "ClusterSizeStd": np.std(clusters_data["Size"]),
                          "ClusterMax": clusters_data["Size"].max(), 
                          "AverageL": np.average(clusters_data["R_Length"]),
                          "MaxL": np.max(clusters_data['R_Length']),
                          "StdL": np.std(clusters_data["R_Length"]),
                          "AverageW": np.average(clusters_data["R_Width"]),
                          "MaxW": np.max(clusters_data['R_Width']),
                          "StdW": np.std(clusters_data['R_Width']),
                          "MaxEEAngleDiff":np.max(clusters_data["EnEchelonAngleDiff"]),
                          "AverageEAngleDiff":np.average(clusters_data["EnEchelonAngleDiff"]),
                          "Date":date},
                          index=[0])
    return evaluation
    
def checkAllClusterChange(lines1, lines2):
    """
    Compare two sets of line clusters to check if they are the same.

    Args:
    lines1 (DataFrame): The first set of line clusters as a DataFrame.
    lines2 (DataFrame): The second set of line clusters as a DataFrame.

    Returns:
    bool: True if the clusters are the same, False otherwise.

    This function compares two sets of line clusters to determine whether they are 
    identical or different. It does this by sorting and hashing the ClusterHash 
    values of both sets and then comparing the resulting hash values. If the hash 
    values are the same, the clusters are considered the same; otherwise, they 
    are considered different.

    Example:
    >>> are_clusters_identical = checkAllClusterChange(cluster_data1, cluster_data2)
    """
    hash1 = np.sort(lines1['ClusterHash'])
    hash2 = np.sort(lines2['ClusterHash'])
    hash1.flags.writeable = False
    hash2.flags.writeable = False

    if hash(str(hash1)) == hash(str(hash2)):
        print("These clusters are the same")
        return True
    else: 
        print("These clusters are not the same")
        return False


def checkIndividualClusterChange(df1, df2):
    """
    Compare individual line clusters between two sets of data frames.

    Args:
    df1 (DataFrame): The first set of line clusters as a DataFrame.
    df2 (DataFrame): The second set of line clusters as a DataFrame.

    Returns:
    tuple: A tuple containing two NumPy arrays - eqLabels and diffLabels.
           - eqLabels: An array of labels that are found in both df1 and df2.
           - diffLabels: An array of labels that are unique to either df1 or df2.

    This function compares individual line clusters between two sets of data frames, 
    df1 and df2. It identifies which cluster labels are common (eqLabels) and which 
    are unique to each data frame (diffLabels). It also provides information about 
    the number of overlapping clusters and the total number of clusters in each data 
    frame.

    Example:
    >>> eqLabels, diffLabels = checkIndividualClusterChange(cluster_data1, cluster_data2)
    """
    l1 = df1['ClusterHash'].values
    l2 = df2['ClusterHash'].values
    HashList = np.array([l1, l2])
    
    # Pick the array with the longest length
    longest = max(HashList, key=lambda col: len(col))
    shortest = min(HashList, key=lambda col: len(col))

    same = np.in1d(longest, shortest)  # Outputs length of first input
    s = np.sum(same)
    eqLabels = np.arange(0, len(longest), 1)[same]
    diffLabels = np.arange(0, len(longest), 1)[~same]
  
    print(s, len(df1), len(df2))
    
    if (s == len(df1) and s == len(df2)):
        print("These clusters are ALL the same")
    else: 
        print("These clusters are not all the same")
        print("df1 has", str(len(df1)), 'while df2 has', str(len(df2)), "identified clusters")
        print("They have", str(s), "overlapping clusters")
        
    return eqLabels, diffLabels




def extendLines(lines, save=False, name='Longlines.csv'):

    """
    Extends lines in dataframe up to min(x)*.5 and max(x)*1.5
        

    Parameters
    ----------
    df : pandas dataframe
        DESCRIPTION.

    Returns
    -------
    lines: pandas.dataframe 
    

    """
    t,r=whichForm(lines)
    xc,yc=HT_center(lines)
    xmid=(lines['Xstart'].values+lines['Xend'].values)/2
    ymid=(lines['Ystart'].values+lines['Yend'].values)/2
    a = np.cos(np.deg2rad(lines[t].values))
    b = np.sin(np.deg2rad(lines[t].values))
    
    x0 = xmid
    y0 = ymid
    
    xmax=np.max([lines['Xstart'].max(), lines['Xend'].max()])*1.5
    xmin=np.min([lines['Xstart'].min(), lines['Xend'].min()])*0.5
    
    l=900000
    
    x1 = (x0 + l/2 * (-b))
    y1 = (y0 + l/2 * (a))
    x2 = (x0 - l/2 * (-b))
    y2 = (y0 - l/2 * (a))
    
    # m=(-1/np.tan(np.deg2rad(lines['AvgTheta'].values)))+0.000001
    # b=(lines['AvgRho'].values)/np.sin(np.deg2rad(lines['AvgTheta'].values))

    
    # print(xmax, xmin)
    # xstart=[xmax]*len(lines)
    # xend=[xmin]*len(lines)
    
    # ystart=m*xstart+b
    # yend=m*xend+b
    longlines=lines.copy()
    longlines['Xstart']=x1
    longlines['Ystart']=y1
    longlines['Xend']=x2
    longlines['Yend']=y2
    
    if save: 
        writeToQGIS(longlines, name)
    
    return longlines


        
def OutputRectangles(clusters):
    """
    Compute the coordinates of bounding rectangles for each cluster in a set of line clusters.

    Args:
    clusters (DataFrame): A DataFrame containing line clusters with attributes 'Labels', 'Xmid', 'Ymid'.

    Returns:
    tuple: A tuple containing two NumPy arrays - Xs and Ys.
           - Xs: An array of X-coordinates for the corners of bounding rectangles for each cluster.
           - Ys: An array of Y-coordinates for the corners of bounding rectangles for each cluster.

    This function computes the coordinates of bounding rectangles for each cluster in a set of 
    line clusters. It uses the 'Xmid' and 'Ymid' attributes of the clusters to determine the 
    center points and calculates the coordinates of the corners of rectangles that enclose 
    the clusters.

    Example:
    >>> Xs, Ys = OutputRectangles(cluster_data)
    """
    clabel = np.unique(clusters['Labels'])
    nclusters = len(clabel)
    xc, yc = HT_center(clusters)
    Xs = np.zeros((nclusters, 5))
    Ys = np.zeros((nclusters, 5))
    
    for i in np.unique(clusters['Labels']): 
        clustered = True
        mask = clusters['Labels'] == i
        lines = clusters[mask]
        
        if (i == -1 or len(lines) < 2):
            clustered = False
            
        w, l, r, Xe, Ye, Xmid, Ymid = fit_Rec(lines, xc, yc)
        Xs[i - 1, :] = Xe
        Ys[i - 1, :] = Ye
        
    return Xs, Ys

                    
                    