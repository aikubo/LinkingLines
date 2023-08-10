#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 22 15:24:34 2023

@author: A Kubo Hutchison (akubo@uoregon.edu)

In order to replicate and explore the results in  Kubo Hutchison et al 
"Multiscale Spatial Patterns in Giant Dike Swarms Identified through Objective Feature Extraction"
we have created this demo file.

This demo is designed to work with the code in  https://github.com/aikubo/Linking-and-Clustering-Dikes
and requires several packages which are included in that repository.

First, please install the correct Conda or Mamba environment. This demo file may not work 
without the correct packages. 

Use the command: 
    conda create --name dikes --file requirements.txt
    
Please check the paths to make sure you have the necessary files in the correct places.

This demo requires the file which contains the Spanish Peaks dike segments which
we have called 'DikeMountain_SpanishPeaks_3857.csv' . This dataset is included in the 
paper supplement. It is a .csv file with the dike locations in a WKT format.

We use the pandas DataFrame data structure extensively and one can familarize themselves 
with the syntax here: https://pandas.pydata.org/docs/ before running the demo.

The variable "dikeset" refers the DataFrame which contains the segment data
while the variable "lines" refers to the clustered lines. 
"""
#previously published packages
import pandas as pd 
import matplotlib.pyplot as plt 
import numpy as np 
import scipy.cluster.hierarchy as sch
import os
import seaborn as sns
from datetime import datetime


# Packages written for the paper
from htMOD import AKH_HT as HT
from htMOD import MidtoPerpDistance, HT_center
from clusterMod import HT_AGG_custom as AggHT
from plotmod import *
from PrePostProcess import * 
from examineMod import examineClusters, checkoutCluster, TopHTSection
from ResultsPlots import *



def DemoKuboHutchisonetAl(save=False):

    """Load and Preprocess dataset
    
    We load the dataset which was digitized using QGIS and then convert the 
    WKT strings to a dataframe with columns [Xstart,Ystart, Xend, Yend, seg_length, midpoint, xc,yc]
    
    We exclude any dikes which cannot be fit to a single straight line. 
    """
    dikeset=pd.read_csv('DikeMountain_SpanishPeaks_3857.csv')
    dikeset=DikesetReProcess(dikeset)

   
    
    """ Hough Transform
    The major input in the hough transform code is the center (xc,yc)
    which is assumed to be the mean of the x and y coordinates in the dataset.
    
    it outputs theta, rho which then are added to the dataset as columns
    """
    
    theta, rho, xc,  yc=HT(dikeset)
    print(xc,yc)
    dikeset= MidtoPerpDistance(dikeset, xc, yc)
    dikeset['theta']=theta
    dikeset['rho']=rho
    
    
    """Clustering
    
    The two main inputs in the clustering algorithm are the theta and rho cut offs
    this determines which segments will be clustered together.
    Details are in clusterMod.py 
    
    The algorithm takes the cutoffs and the dataset as inputs and the linkage, rotation,
    and distance metric can be changed. 
    We rotate the dataset so the average angle is set at 20 degreees this decreases 
    the instances of clustecrossing from -90 to 90. 
    The outputs are not rotated and retain their original orientation. 
    the dataframe now has the column label which is the cluster ID
    
    AggHT also output Z which is the cluster linkage matrix which can be helpful 
    if you need to rerun the clustering.
    """
    dtheta=2 
    drho=np.floor(dikeset['seg_length'].mean())
    
    dikeset, Z=AggHT(dikeset, dtheta, drho, linkage='complete', rotate=True, metric='Euclidean')
    lines,evaluation=examineClusters(dikeset) #this function takes all the data with the labels 
    # and combines them into lines 

    
    dikeset=dikeset.assign(Date_Changed=date) #add date chaned to dataframe
    lines=lines.assign(Date_Changed=date)
    lines=lines.assign(Rho_Threshold=drho)
    lines=lines.assign(Theta_Threshold=dtheta)

    if save:
        #save the outputs
        name=j+str(int(drho))+".csv"
        npname=j+str(int(drho))+"LINKAGE"+".npy"
        now = datetime.now() 
        date = now.strftime("%d %b, %Y")
        date2=now.strftime('%d_%m_%Y')
        lines.to_csv(name, index=False)
        dikeset.to_csv(i, index=False)
    
    """#Plotting Radial Centers"""
    fig2,Cart=plt.subplots()
    
    fig = SetupJGRFig((90,170), 'landscape') #set up figure size
    dist_width=0.10
    gskw = dict(height_ratios= [ dist_width,1-dist_width], width_ratios=[(1-dist_width)/2,(1-dist_width)/2,dist_width])
    
    gs = gridspec.GridSpec(2, 3, **gskw)
    ht1 = plt.subplot(gs[:,0])
    ht2= plt.subplot(gs[1,1])
    tdist=plt.subplot(gs[0,1])
    rdist=plt.subplot(gs[1,2])
    

    
    theta=np.linspace(-90,90)
   
    """Fitting Radial Centers
    This uses functions in fitRadialCenters.py to fit data to the equation 
    
    rho= X*cos(theta)+Y*sin(theta)
    where (X,Y) is the location of a radial center in cartesian coordinates
    
    Radial fit returns a pandas dataframe of the location wiht other metrics 
    such as how many segments intersect with the center.
    """
    xc,yc=HT_center(dikeset) #find the HT center of the data
    DotsHT(fig, ht1, lines,  ColorBy='Ymid', cmap='magma') #plot the theta and rho data
    DotsHT(fig, ht2, lines, alpha=0.05, ColorBy=None, axlabels=(True,False)) #plot the theta and rho data
    plotlines(lines, 'k', Cart, alpha=0.2) #plot the clustered lines 
    

    #purple
    mask2=(dikeset['Ymid']>4.52e6)  #filter dataset
    Purp=RadialFit(dikeset[mask2]) #fit based on filtered dataset
    PurpClose, Purp=NearCenters(lines, Purp, tol=2500) # pickout all clustered lines which intersect within 2.5km of center
    PurpCloseSeg, Purp=NearCenters(dikeset, Purp, tol=2500) # pickout all segments which intersect within 2.5km of center
    rhoPurp=CenterFunc(theta,Purp['Center'][0][0], Purp['Center'][0][1], xc,yc)/1000 #calculate line for plotting
    ht2.plot(theta,rhoPurp, color='purple', linewidth=2) #plot radial line over data
    DotsHT(fig, ht2, PurpCloseSeg, color='purple', ColorBy=None, axlabels=(True,False)) #plot the theta an rho data
    plotlines(PurpClose, 'purple', Cart) #plot the clustered lines in purple
    
    #green
    mask2=(dikeset['Ymid']<4.52e6) & (dikeset['Ymid']>4.48e6)
    green=RadialFit(dikeset[mask2])
    GreenClose, green=NearCenters(lines, green, tol=2500)
    GreenCloseSeg, green=NearCenters(dikeset, green, tol=2500)
    rhoGreen=CenterFunc(theta,green['Center'][0][0], green['Center'][0][1], xc,yc)/1000
    ht2.plot(theta,rhoGreen, color='green', linewidth=2)
    DotsHT(fig, ht2, GreenCloseSeg, color='green', ColorBy=None, axlabels=(True,False))
    plotlines(GreenClose, 'green', Cart)
    
    radLabels=np.concatenate( (PurpCloseSeg['HashID'].values, GreenCloseSeg['HashID'].values))
    lin=dikeset[(~np.in1d(dikeset['HashID'].values, radLabels)) & ( (dikeset['theta'].values<-55) )]
    print(len(lin))
    DotsHT(fig, ht2, lin, color='blue', ColorBy=None, axlabels=(True,False))
    
    lin=lin.assign(Structure='Linear')
    GreenCloseSeg=GreenCloseSeg.assign(Structure='Radial 1')
    PurpCloseSeg=PurpCloseSeg.assign(Structure='Radial 2')
    labelSeg=pd.concat((lin, GreenCloseSeg, PurpCloseSeg))
    
    
    """Histogram of rho and theta"""
    rbins=np.arange( dikeset['rho'].min(), dikeset['rho'].max(), 5000)
    
    rdist.hist(lin['rho'], color='blue', alpha=0.6, bins=rbins,  orientation='horizontal')
    rdist.hist(GreenCloseSeg['rho'], color='green', bins=rbins,  alpha=0.6, orientation='horizontal')
    rdist.hist(PurpCloseSeg['rho'], color='purple', bins=rbins,   alpha=0.6, orientation='horizontal')
    
    tbins=np.arange(-90,90,10)
    tdist.hist(lin['theta'], color='blue', alpha=0.6, bins=tbins) #, density=True)
    tdist.hist(GreenCloseSeg['theta'], color='green', alpha=0.6, bins=tbins) # , density=True)
    tdist.hist(PurpCloseSeg['theta'], color='purple', alpha=0.6, bins=tbins) #, density=True)
    
    tdist.tick_params(labelbottom = False, bottom = False)
    #tdist.yaxis.set_major_formatter(PercentFormatter(xmax=1))
    rdist.tick_params(labelleft = False, left = False)
    #rdist.xaxis.set_major_formatter(PercentFormatter(xmax=1))
    ht2.tick_params(labelleft = False, left = False)
    
    """Labelling dataframes"""
    radLabelsl=np.concatenate( (PurpClose['Label'].values, GreenClose['Label'].values))
    linl=lines[((~np.in1d(lines['Label'].values, radLabelsl)) & (lines['AvgTheta'].values<-55)) ]
    linl=linl.assign(Structure='Linear')
    GreenClose=GreenClose.assign(Structure='Radial 1')
    PurpClose=PurpClose.assign(Structure='Radial 2')
    
    GreenCloseSeg=GreenCloseSeg.assign(Structure='Radial 1')
    PurpCloseSeg=PurpCloseSeg.assign(Structure='Radial 2')
    
    labelLines=pd.concat((linl, GreenClose, PurpClose))
    labelSubplots([ht1, ht2, tdist, rdist])
    plt.tight_layout()
    
    Centers=pd.concat((green, Purp)) #make one dataframe called centers
    
    CloseLines=pd.concat((GreenClose, PurpClose, linl)) #make one dataframe called closelines
    CloseSegments=pd.concat((GreenCloseSeg, PurpCloseSeg))
    if save:
        writeCenterWKT(Centers, 'SpanishPeaksRadialFitsCenters.csv')
        writeToQGIS(CloseLines, 'SPRadialFitsLines.csv')
        writeToQGIS(CloseSegments, 'SPRadialFitsSegments.csv')
    
        fig.savefig('SPRadialFitsLinesFig.pdf', dpi=600)

    