#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul  7 14:45:57 2022

@author: akh
"""

from plotmod import DotsLines, whichForm
import numpy as np 
from PrePostProcess import transformXstart
import matplotlib.pyplot as plt


def dilation(df, binWidth=1, averageWidth=1, method='Expanded'):
    """
    Calculate dilation values for a given dataset along EW and NS directions.

    This function calculates dilation values along the East-West (EW) and North-South (NS) directions based on the input DataFrame. It divides the data into bins and computes the dilation for each bin. The method parameter allows you to choose between different dilation calculation methods.

    Parameters:
        df (pandas.DataFrame): The input DataFrame containing data.
        binWidth (float, optional): The width of bins used for dilation calculation. Default is 1.
        averageWidth (float, optional): The average width used in dilation calculation. Default is 1.
        method (str, optional): The dilation calculation method. Options are 'Average', 'Total', and 'Expanded'. Default is 'Expanded'.
        
 The differences between the three dilation calculation methods (Expanded, Average, and Total) lie in how they compute dilation values from the data. Dilation is a measure of how much a geological feature has been stretched or expanded, typically expressed as a ratio. Let's explain the differences between these methods:

1. **Expanded Dilation:**
   
   - **Method Description:** Expanded dilation calculates dilation by considering how much a line segment expands within a bin. It accounts for segments that partially cross the bin's boundaries and calculates dilation accordingly.

   - **Calculation:** For each bin, it sums the fractional expansions of line segments within the bin. Fractional expansion is the portion of a line segment that enters the bin. It considers both the portion of segments that start or end within the bin and the portions that cross into the bin.

   - **Advantages:** Expanded dilation provides a more detailed measure of dilation because it considers partial segments within bins. It can capture dilation from segments that only partially cross the bin boundaries.

   - **Use Cases:** Expanded dilation is useful when you want to capture the full contribution of line segments to dilation, including those that partially overlap with the bin.

2. **Average Dilation:**

   - **Method Description:** Average dilation calculates dilation by averaging the dilation values of all line segments entirely contained within each bin.

   - **Calculation:** For each bin, it calculates the dilation value for each segment entirely contained within the bin. It then averages these segment dilation values to obtain the bin's dilation value.

   - **Advantages:** Average dilation provides a simpler and smoother measure of dilation. It may reduce noise caused by the partial overlap of line segments with bin boundaries.

   - **Use Cases:** Average dilation is suitable when you want a smoother representation of dilation and are less concerned about capturing the effects of partial segment overlap.

3. **Total Dilation:**

   - **Method Description:** Total dilation calculates dilation by summing the dilation values of all line segments entirely contained within each bin.

   - **Calculation:** Similar to Average Dilation, it calculates the dilation value for each segment entirely contained within the bin. Instead of averaging, it sums these segment dilation values to obtain the bin's dilation value.

   - **Advantages:** Total dilation provides a cumulative measure of dilation within each bin. It considers the total dilation contributed by all segments within the bin.

   - **Use Cases:** Total dilation is useful when you want to understand the cumulative effect of line segments on dilation within each bin. It may be preferred when analyzing the combined impact of multiple segments.

In summary, the choice between Expanded, Average, or Total dilation depends on your specific geological analysis needs. Expanded dilation is more detailed and suitable for capturing partial segment contributions. Average dilation provides a smoother representation, while Total dilation gives a cumulative measure of dilation within each bin. The choice should be based on the level of detail and accuracy required for your analysis.
    Returns:
        numpy.ndarray: The EW dilation values.
        numpy.ndarray: The NS dilation values.
        numpy.ndarray: The bin edges along the EW direction.
        numpy.ndarray: The bin edges along the NS direction.

    Example:
        # Calculate dilation using the 'Expanded' method
        EW_Dilation, NS_Dilation, binx, biny = dilation(df, binWidth=1, averageWidth=1, method='Expanded')
    """
    # Function code goes here...

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
    Create a triple-panel plot showing dilation results.

    This function creates a triple-panel plot to visualize dilation results. It includes a main panel with line segments
    in Cartesian coordinates on the left, a histogram of NS dilation values on the bottom, and a histogram of EW dilation
    values on the right.

    Parameters:
        df (pandas.DataFrame): A DataFrame containing data points.
        lines (pandas.DataFrame): A DataFrame containing line segment data.
        shape (list, optional): A list specifying the final figure size and orientation.
            Defaults to ['half', 'portrait'].
        kwargs (dict, optional): Additional keyword arguments to pass to the `dilation` function for calculating dilation.
            Defaults to None.

    Returns:
        matplotlib.figure.Figure: The modified Figure object.
        list of matplotlib.axes._subplots.AxesSubplot: A list of three axes objects (main panel, NS dilation histogram, and EW dilation histogram).

    Example:
        # Create a triple-panel dilation plot
        fig, ax = TripleDilationPlot(data_df, lines_df, shape=['half', 'portrait'], kwargs={'rstep': 0.1, 'tstep': 10})
    """


    fig = SetupJGRFig(shape[0], shape[1])
    
    
   
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
    
    
    
    #ax_yDist.plot(EWDilation2,biny2, color='k', label='Filtered' )
    
    #ax_xDist.plot(binx2,NSDilation2, color='k' )
    # EWorder=np.argsort([np.sum(i) for i in EWl])
    # colors=['gray', 'green', 'yellow']
    # for i in EWorder:
    #     ax_yDist.fill_between(EWl[i],0,yl[i], alpha=0.6, color=colors[i] )
        
    # NSorder=np.argsort([np.sum(i) for i in NSl])

    # for i in NSorder:
    #     ax_xDist.fill_between(xl[i],0, NSl[i], alpha=0.6, color=colors[i])
        
    
    
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


