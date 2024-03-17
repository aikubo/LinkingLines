# LinkingLines Package
 # Written by aikubo
 # Version: 2.1.0
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul  7 14:45:57 2022

@author: akh
"""

from .PlotUtils import DotsLines, whichForm, plotlines, SetupAGUFig, labelSubplots
import numpy as np
from .PrePostProcess import transformXstart
import matplotlib.pyplot as plt


def dilation(df, binWidth=1.0, averageWidth=1.0, method='Expanded'):
    """
    Calculates dilation values for a dataset along East-West (EW) and North-South (NS) directions.

    Dilation is calculated by dividing the data into bins and computing the dilation for each bin based on the specified method. 
    There are three methods available for calculating dilation: 'Expanded', 'Average', and 'Total'. 
    - 'Expanded' considers partial and fully contained line segments in each bin, calculating fractional expansions. 
    - 'Average' calculates the average dilation of fully contained line segments within each bin. 
    - 'Total' sums the dilation values of fully contained line segments within each bin.

    Parameters
    ----------
    df : pandas.DataFrame
        The input DataFrame containing data to calculate dilation. Must include columns for line segment coordinates.
    binWidth : float, optional
        The width of bins used for dilation calculation along both EW and NS directions. Default is 1.0.
    averageWidth : float, optional
        The average width of line segments used in the dilation calculation. Default is 1.0.
    method : {'Expanded', 'Average', 'Total'}, optional
        The method used for calculating dilation:
        - 'Expanded': Considers partial and fully contained segments, calculating fractional expansions.
        - 'Average': Averages dilation values of fully contained segments within each bin.
        - 'Total': Sums dilation values of fully contained segments within each bin.
        Default is 'Expanded'.

    Returns
    -------
    EWDilation : numpy.ndarray
        An array of calculated dilation values along the East-West direction.
    NSDilation : numpy.ndarray
        An array of calculated dilation values along the North-South direction.
    binx : numpy.ndarray
        The bin edges along the East-West direction.
    biny : numpy.ndarray
        The bin edges along the North-South direction.


    Notes
    -----
    The choice of method ('Expanded', 'Average', 'Total') depends on the specific requirements of the geological analysis. 
    'Expanded' provides detailed calculations including partial segments, 'Average' offers smoother dilation values, 
    and 'Total' gives a cumulative measure of dilation within bins.
    """

    t,r=whichForm(df)
    df=transformXstart(df)

    xs=[ min(df['Xstart'].min(), df['Xend'].min() ), max( df['Xstart'].max(), df['Xend'].max() )]
    ys=[ min(df['Ystart'].min(), df['Yend'].min() ), max( df['Ystart'].max(), df['Yend'].max() )]

    if np.ptp(xs) < binWidth or np.ptp(ys) < binWidth:
        binWidth=np.min( [np.ptp(xs)/10, np.ptp(ys)/10])

    binx=np.arange(xs[0]-binWidth, xs[1]+binWidth, binWidth)
    biny=np.arange(ys[0]-binWidth, ys[1]+binWidth, binWidth)

    FractionEW=np.cos(np.deg2rad(df[t].values))
    # for theta=0, Fraction EW is 1, for theta=90 fractonEW is 0
    FractionNS=abs(np.sin(np.deg2rad(df[t].values)))
    # its the opposite of above

    Xstart=df['Xstart'].values
    Ystart=df['Ystart'].values
    Xend=df['Xend'].values
    Yend=df['Yend'].values

    y=np.array([Ystart,Yend]).T
    Ystart=np.min(y, axis=1)
    Yend=np.max(y, axis=1)+1


    x=np.array([Xstart,Xend]).T
    Xstart=np.min(x, axis=1)
    Xend=np.max(x, axis=1)+1
    # for now sort Yend and Ystart, can't keep it but for now it's okay



    EWDilation=np.zeros_like(biny)
    NSDilation=np.zeros_like(binx)

    TotalEWDilation=0
    TotalNSDilation=0


    for i in range(len(binx)-1):

        # Case 1: Line Passes completely thru bin
        maskx1= np.logical_and( (Xstart< binx[i]), (Xend>binx[i+1]))

        #Case 2: Starts in Bin
        maskx2=np.logical_and( (Xstart> binx[i]), (Xstart<binx[i+1]) )

        #Case3: Ends in Bin
        maskx3=np.logical_and( (Xend> binx[i]), (Xend<binx[i+1]) )

        maskx=np.logical_or(maskx1, np.logical_or(maskx2, maskx3))

        #print(np.sum(maskx), "in binx",binx[i], "-", binx[i+1] )
        l=np.where(binx==binx[i])
        #print(np.sum(maskx))
        if np.isclose(np.sum(maskx), 0.0):
            NSDilation[l]=0
        elif method == 'Average':
            NSDilation[l]= np.mean(FractionNS[maskx])*averageWidth*np.sum(maskx)
        elif method== 'Total':
            NSDilation[l]= np.sum(FractionNS[maskx])*averageWidth
        elif method =='Expanded':
             NSDilation[l]= (np.sum(FractionNS[maskx1]) +
                             np.sum( (abs(Xstart[maskx2]-binx[i])/binWidth)*FractionNS[maskx2]) +
                             np.sum( (abs(Xend[maskx3]-binx[i])/binWidth)*FractionNS[maskx3]))*averageWidth

    for j in range(len(biny)-1):
         # Case 1: Line Passes completely thru bin
        masky1= np.logical_and( (Ystart< biny[j]), (Yend>biny[j+1]))

        #Case 2: Starts in Bin
        masky2=np.logical_and( (Ystart> biny[j]), (Ystart<biny[j+1]) )

        #Case3: Ends in Bin
        masky3=np.logical_and( (Yend> biny[j]), (Yend<biny[j+1]) )

        masky=np.logical_or(masky1, np.logical_or(masky2, masky3))

        #print(np.sum(masky), "in biny",biny[j], "-", biny[j+1] )
        l=np.where(biny==biny[j])
        #print(np.sum(masky))
        if np.isclose(np.sum(masky), 0.0):
            EWDilation[l]=0
        elif method == 'Average':
            EWDilation[l]= np.mean(FractionEW[masky])*averageWidth*np.sum(masky)
        elif method== 'Total':
            EWDilation[l]= np.sum(FractionEW[masky])*averageWidth
        elif method =='Expanded':
            EWDilation[l]= (np.sum(FractionEW[masky1]) +
                            np.sum( (abs(Ystart[masky2]-biny[j])/binWidth)*FractionEW[masky2]) +
                            np.sum( (abs(Yend[masky3]-biny[j])/binWidth)*FractionEW[masky3]))*averageWidth



    return EWDilation, NSDilation, binx, biny


def TripleDilationPlot(df, lines, shape=['half', 'portrait'], kwargs=None):
    """
    Creates a triple-panel plot to visualize dilation results.

    This function generates a visualization consisting of three panels: a main plot displaying line segments
    in Cartesian coordinates, a histogram of North-South (NS) dilation values below the main plot, and a histogram
    of East-West (EW) dilation values to the right of the main plot. It utilizes the `dilation` function to calculate
    dilation values based on the provided line segment data.

    Parameters
    ----------
    df : pandas.DataFrame
        A DataFrame containing data points to be visualized in the main panel.
    lines : pandas.DataFrame
        A DataFrame containing line segment data used for dilation calculation.
    shape : list, optional
        A list specifying the final figure size and orientation. The first element determines the figure's size ('full', 
        'half', 'quarter'), and the second element specifies the orientation ('landscape', 'portrait'). 
        Defaults to ['half', 'portrait'].
    kwargs : dict, optional
        Additional keyword arguments to pass to the `dilation` function for calculating dilation. 
        This allows for customization of the dilation calculation, such as bin width and method. Defaults to None.

    Returns
    -------
    fig : matplotlib.figure.Figure
        The figure object containing the generated triple-panel plot.
    ax : list of matplotlib.axes._subplots.AxesSubplot
        A list containing the three axes objects for the main panel, NS dilation histogram, and EW dilation histogram.

    Notes
    -----
    - The main panel plot displays line segments as provided in the `lines` DataFrame.
    - Dilation histograms help visualize the distribution of dilation values across NS and EW directions, providing insight
      into the stretching or expansion patterns in the dataset.
    - Customization options for the dilation calculation are provided via the `kwargs` parameter, allowing the user to 
      adjust aspects such as bin width and calculation method according to their analysis needs.
    """

    fig = SetupAGUFig(shape[0], shape[1])



    xlim=[ np.min( [lines['Xstart'].min(), lines['Xend'].min()]), np.max( [lines['Xstart'].max(), lines['Xend'].max()])]
    ylim=[ np.min( [lines['Ystart'].min(), lines['Yend'].min()]), np.max( [lines['Ystart'].max(), lines['Yend'].max()])]
    aspect=np.diff(xlim)[0]/np.diff(ylim)[0]
    #aspect is w/h, xlim is w and ylim is h

    gskw = dict(width_ratios = [ .75, .25],
                height_ratios= [ .25,.75])

    gs = gridspec.GridSpec(2, 2, **gskw)
    ax_main = plt.subplot(gs[1, 0])



    ax_xDist = plt.subplot(gs[0, 0], sharex=ax_main, adjustable='box')
    ax_yDist = plt.subplot(gs[1, 1], sharey=ax_main, adjustable='box')




    EWDilation, NSDilation, binx,biny=dilation(df, **kwargs)
    EWDilation1, NSDilation1, binx1,biny1=dilation(lines, **kwargs)

    m=lines['TrustFilter']==1

    EWDilation2, NSDilation2, binx2,biny2=dilation(lines[m], **kwargs)

    ys=[np.min(biny), np.max(biny)]


    #ax_xDist.plot(binx[:-1], NSDilation[:-1])

    ax_xDist.fill_between(binx1,0, NSDilation1, alpha=0.6, color='r')
    ax_xDist.fill_between(binx,0, NSDilation, alpha=0.6, color='b')

    f1=ax_yDist.fill_between(EWDilation1,0,biny1, alpha=0.6, color='r', label='All Linked' )
    f2=ax_yDist.fill_between(EWDilation,0,biny, alpha=0.6, color='b', label='Segments' )


    ax_xDist.set(ylabel='NS Dilaton (m)')
    ax_yDist.set(xlabel='EW Dilaton (m)')

    ax_xDist.tick_params(axis='x',          # changes apply to the x-axis
                        which='both',      # both major and minor ticks are affected
                        bottom=False,      # ticks along the bottom edge are off
                        top=False,         # ticks along the top edge are off
                        labelbottom=False) # labels along the bottom edge are off)
    ax_yDist.tick_params(axis='y',
                         which='both',
                         left=False,
                         right=False,
                         labelleft=False)

    ax_yDist.set_ylim([ys[0], ys[1]])
    labelSubplots([ax_main, ax_xDist, ax_yDist])
    plotlines(lines, 'r', ax_main, alpha=0.6)
    plotlines(df, 'b', ax_main, alpha=0.6)
    yMean=np.average(biny2, weights=EWDilation2)
    xMean=np.average(binx2, weights=NSDilation2)
    print(yMean, xMean)
    m=ax_yDist.axhline(y=yMean, color='navy', linestyle=":", label='Mean')
    ax_xDist.axvline(x=xMean,color='navy', linestyle=":")


    ax_main.axhline(y=yMean, color='navy', linestyle=":")
    ax_main.axvline(x=xMean,color='navy', linestyle=":")


    ax_yDist.legend(loc="lower left")
    print("")
    print("EW")
    print("Max segment")
    print(np.max(EWDilation))
    print("Max Linked")
    print(np.max(EWDilation1))
    print("Max Filtered")
    print(np.max(EWDilation2))
    print("XRange (m)")
    print( xlim[1]-xlim[0])
    print('max strain')
    print("segment")
    print(np.max(EWDilation)/( xlim[1]-xlim[0])*100)
    print("Linked")
    print(np.max(EWDilation1)/( xlim[1]-xlim[0])*100)
    print("Filtered")
    print(np.max(EWDilation2)/( xlim[1]-xlim[0])*100)
    print("")
    print("NS")
    print("Max segment")
    print(np.max(NSDilation))
    print("Max Linked")
    print(np.max(NSDilation1))
    print("Max Filtered")
    print(np.max(NSDilation2))
    print("YRange (m)")
    print( ylim[1]-ylim[0])
    print('max strain')
    print("segment")
    print(np.max(NSDilation)/( ylim[1]-ylim[0])*100)
    print("Linked")
    print(np.max(NSDilation1)/( ylim[1]-ylim[0])*100)
    print("Filtered")
    print(np.max(NSDilation2)/( ylim[1]-ylim[0])*100)
    print("")

    return fig, [ax_main, ax_xDist, ax_yDist]
