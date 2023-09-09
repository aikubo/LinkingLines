#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 20 16:42:50 2022

@author: akh
"""

from PrePostProcess import * 
import pandas as pd 
import numpy as np 
from htMOD import AKH_HT as HT 
from htMOD import HT_center, MidtoPerpDistance, moveHTcenter
from plotmod import *
import matplotlib.pyplot as plt 
from clusterMod import HT_AGG_custom
import scipy.cluster.hierarchy as sch
from examineMod import examineClusters, checkoutCluster, TopHTSection, checkoutClusterCart
import os
import labellines
from matplotlib.patches import Rectangle
from matplotlib.gridspec import GridSpec
import matplotlib.colors as colors
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
from fitRadialCenters import RadialFit, CenterFunc, NearCenters, writeCenterWKT, RadialAzimuthal, RipleyRadial
from matplotlib.ticker import PercentFormatter
from mpl_toolkits.mplot3d import Axes3D
from synthetic import fromHT, makeLinear2, makeRadialSwarmdf


import seaborn as sns

# def WriteSIfiles():
#     datasets=['/home/akh/myprojects/Linking-and-Clustering-Dikes/dikedata/spanish peaks/SpanishPeaks_Complete_euc_3_2013.csv',
#               'dikedata/deccandata/AllDeccanLinked_euc_18_10_2022.csv',
#               '/home/akh/myprojects/Linking-and-Clustering-Dikes/dikedata/crb/AllCRBLinked_euc_18_10_2022.csv']
    
#             'Xstart', 'Ystart', 'Xend', 'Yend', 'X0', 'Y0', 'AvgRho',
#            'AvgTheta', 'AvgSlope', 'AvgIntercept', 'RhoRange', 'Aspect', 'Xmid',
#            'Ymid', 'ThetaRange', 'StdRho',
#            'StdTheta', 'Size', 'R_error', 'Linked',
#            'SegmentLSum', 'ClusterHash', 'ClusterCrossesZero',
#            'EnEchelonAngleDiff', 'Overlap', 'nOverlapingSegments', 'EEPvalue',
#            'MaxSegNNDist', 'MedianSegNNDist', 'MinSegNNDist', 'TrustFilter',
#            'Date_Changed', 'Rho_Threshold', 'Theta_Threshold',
#            'Dike Cluster Width (m)', 'Dike Cluster Length (km)', 'Average Rho (m)',
#            'Average Theta ($^\circ$)', 'SwarmID', 'HashID',
#            'yc', 'xc'

def HTFigure1():
    SMALL_SIZE = 8
    MEDIUM_SIZE = 8
    BIGGER_SIZE = 12
    dikeset=pd.read_csv('/home/akh/myprojects/Linking-and-Clustering-Dikes/dikedata/AllDatasets_18_10_2022.csv')
    plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
    plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
    plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
    plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
    plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
    plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
    plt.rc('figure', titlesize=MEDIUM_SIZE)  # fontsize of the figure title
    plt.rcParams['legend.title_fontsize'] = SMALL_SIZE
    



    fig,ax=plt.subplots(3,1)
    fig.set_size_inches( 100/25.4, 225/25.4,)
    #temp=pd.read_csv('/home/akh/myprojects/Linking-and-Clustering-Dikes/dikedata/deccandata/AllDeccan_PreProcessed.csv')
    #DotsHT(fig, ax[0],temp, ColorBy='SwarmID', palette="hls", rhoScale=True)
    

    for i, p in enumerate([ 'Deccan','CRBG', 'SpanishPeaks']):
        temp=dikeset[dikeset['Region']==p]
        if p=='Deccan':
            temp=pd.read_csv('/home/akh/myprojects/Linking-and-Clustering-Dikes/dikedata/deccandata/AllDeccan_PreProcessed.csv')
        theta,rho, x, y=HT(temp)
        temp=temp.assign(theta=theta, rho=rho, abstheta=abs(theta))

        ax[i].set_xlabel('Theta ($^\circ$)')
        ax[i].set_ylabel('Rho (km)')
        
        
        sns.scatterplot(temp, x='theta', y=rho/1000, hue='abstheta', style='SwarmID', legend=False, palette='viridis', ax=ax[i], alpha=0.4, edgecolor="black")
        #DotsHT(fig, ax[j],temp, ColorBy='SwarmID', palette="hls", rhoScale=True)
        #ax[i].set_ylabel(None)
        #ax[i].tick_params(axis="y",direction="in", pad=-22)
        
    plt.tight_layout()
    ax[0].set_ylabel('Rho (km)')
    labelSubplots(ax,fontsize=12, labels=['g', 'h', 'i'])
    fig.savefig('/home/akh/myprojects/Linking-and-Clustering-Dikes/dikedata/All_HT3_coloredbyThetaLarge.pdf', dpi=600)
    fig.savefig('/home/akh/myprojects/Linking-and-Clustering-Dikes/dikedata/All_HT3_coloredbyThetaLarge.png', dpi=600)
    
    fig,ax=plt.subplots(3,1)
    fig.set_size_inches( 110/25.4, 225/25.4,)
    #temp=pd.read_csv('/home/akh/myprojects/Linking-and-Clustering-Dikes/dikedata/deccandata/AllDeccan_PreProcessed.csv')
    #DotsHT(fig, ax[0],temp, ColorBy='SwarmID', palette="hls", rhoScale=True)

    for i, p in enumerate([ 'Deccan','CRBG', 'SpanishPeaks']):
        temp=dikeset[dikeset['Region']==p]
        if p=='Deccan':
            temp=pd.read_csv('/home/akh/myprojects/Linking-and-Clustering-Dikes/dikedata/deccandata/AllDeccan_PreProcessed.csv')
        theta,rho, x, y=HT(temp)
        temp=temp.assign(theta=theta, rho=rho, abstheta=abs(theta))
        ax[i].set_xlabel('Theta ($^\circ$)')
        ax[i].set_ylabel('Rho (km)')
        
        sns.scatterplot(temp, x='theta', y=rho/1000, hue='abstheta', style='SwarmID', legend='brief', palette='viridis', ax=ax[i], alpha=0.4, edgecolor="black")
        #DotsHT(fig, ax[j],temp, ColorBy='SwarmID', palette="hls", rhoScale=True)
        #ax[i].set_ylabel(None)
        ax[i].tick_params(axis="y",direction="in")
        
    plt.tight_layout()
    ax[0].set_ylabel('Rho (km)')
    fig.savefig('/home/akh/myprojects/Linking-and-Clustering-Dikes/dikedata/All_HT3_coloredbyTheta_wlegendLarge_portrait.pdf', dpi=600)

    
    
def PubBoxnWhisker(df, p='icefire', hue='Status', l='Length (km)'):
    
    SMALL_SIZE = 8.5
    MEDIUM_SIZE = 10
    BIGGER_SIZE = 12
    
    plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
    plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
    plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
    plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
    plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
    plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
    plt.rc('figure', titlesize=MEDIUM_SIZE)  # fontsize of the figure title
    plt.rcParams['legend.title_fontsize'] = SMALL_SIZE
    
    fig = plt.figure(figsize=(230/25.4, 190/25.4), constrained_layout=True)

    gs = GridSpec(2, 4, figure=fig)
    #lengths
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[0, 1:])

    ax1.set_yscale('log')
    ax2.set_xscale('log')

    g=sns.boxplot(x='Region', y=l,  hue=hue, data=df, ax=ax1, palette=p)


    g.legend_.remove()

    g=sns.boxplot(x=l, y='SwarmID',  hue=hue, data=df, ax=ax2, palette=p)
    # sns.stripplot(x='R_Length', y='SwarmID',data=dfAll, ax=ax[1],
    #               size=.8, color=".05", linewidth=0)

    ax2.xaxis.grid(True)
    ax2.set(ylabel="")
    
    #widths
    ax3 = fig.add_subplot(gs[1, 0])
    ax4 = fig.add_subplot(gs[1, 1:])

    ax3.set_yscale('log')
    ax4.set_xscale('log')
    
    df2=df.loc[(df['Status']=='Filtered Cluster')]# & (np.in1d(df['SwarmID'].values, ['Deccan:Central','Deccan:Narmada-Tapi','Deccan:Saurashtra','CRBG:CJDS']))]
    
    g=sns.boxplot(x='Region', y='Dike Cluster Width (m)', data=df2, ax=ax3, palette=p)

    #g.legend_.remove()

    g=sns.boxplot(x='Dike Cluster Width (m)', y='SwarmID', data=df2, ax=ax4, palette=p)
    #sns.move_legend(g, "lower left")
    # sns.stripplot(x='R_Length', y='SwarmID',data=dfAll, ax=ax[1],
    #               size=.8, color=".05", linewidth=0)
    ax4.xaxis.grid(True)
    ax4.set(ylabel="")
    
    
    sns.despine(offset=5)
    
    
    labelSubplots([ax1, ax2, ax3, ax4], fontsize=12)
    
    return fig, [ax1, ax2, ax3, ax4]
    

def PubWidthsVsLengths(dfAllTrusted, palette="icefire", hue='Region', style='Region'):
        
    SMALL_SIZE = 8.5
    MEDIUM_SIZE = 10
    BIGGER_SIZE = 12
    
    plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
    plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
    plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
    plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
    plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
    plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
    plt.rc('figure', titlesize=MEDIUM_SIZE)  # fontsize of the figure title
    plt.rcParams['legend.title_fontsize'] = SMALL_SIZE
    
    fig, [ax_main, ax_xDist, ax_yDist]=plotScatterHist(dfAllTrusted, "Dike Cluster Width (m)", 
                                                       "Dike Cluster Length (km)",
                                                       hue=hue, log_scale=(True, True),
                                                       palette=palette,
                                                       style=style,
                                                       size='quarter')
    rat1=1 #10^3
    ax_main, l1=plotRatioLine(ax_main, dfAllTrusted['Dike Cluster Width (m)'].values, rat1, line_kw=dict(color='b', ls=":", label="$10^{3}$"))
    rat2=.1 #10^2
    ax_main, l2=plotRatioLine(ax_main, dfAllTrusted['Dike Cluster Width (m)'].values, rat2, line_kw=dict(color='b', ls=":", label="$10^{4}$"))
    rat3=.01 #10^1
    ax_main, l3=plotRatioLine(ax_main, dfAllTrusted['Dike Cluster Width (m)'].values, rat3, line_kw=dict(color='b', ls=":", label="$10^{5}$"))
    ax_main.set_ylim((1,1500))
    ax_main.set_xlim((5,30000))
    a=[l1, l2, l3]
    
    labellines.labelLine(
        l1[0],
        30,
        label=r"$10^{3}$",
        ha="left",
        va="bottom",
        align=False,
        backgroundcolor="none",
        fontsize=8.5
    )
    
    labellines.labelLine(
        l2[0],
        45,
        label=r"$10^{2}$",
        ha="left",
        va="bottom",
        align=False,
        backgroundcolor="none",
        fontsize=8.5
    )
    
    labellines.labelLine(
        l3[0],
        300,
        label=r"$10^{1}$",
        ha="left",
        va="bottom",
        align=False,
        backgroundcolor="none",
        fontsize=8.5
    )
    labelSubplots([ax_main, ax_xDist, ax_yDist], fontsize=12)
    
    return fig, [ax_main, ax_xDist, ax_yDist]

def PubTwistOverlap(df, partial=True, return_EE=False, tol=.10):
        
    SMALL_SIZE = 8.5
    MEDIUM_SIZE = 10
    BIGGER_SIZE = 12
    
    plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
    plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
    plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
    plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
    plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
    plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
    plt.rc('figure', titlesize=MEDIUM_SIZE)  # fontsize of the figure title
    plt.rcParams['legend.title_fontsize'] = SMALL_SIZE
    
    print(len(df))
    df=df.assign(Overlap_Normalized= (df['Overlap'].values*df['SegmentLSum'].values/df['Size'].values+0.001), 
                 TwistAngle=df['EnEchelonAngleDiff'].values+0.001,
                 b=df['SegmentLSum'].values/df['Size'].values/2)
    #only plot ones with twist >0 and overlap >0 
    if partial:
        df=df[ (df['Overlap_Normalized']>0.1) & (df['TwistAngle']>1)]
    
    print(len(df))
    norm=colors.LogNorm(vmin=df['b'].min(), vmax=df['b'].max())
    fig, [ax_main, ax_xDist, ax_yDist]= plotScatterHist(df, 'TwistAngle',
                                                        'Overlap_Normalized', 
                                                        hue="b",
                                                        log_scale=(True,True),
                                                        palette='Blues',
                                                        hue_norm=norm)
    ax_main.set_ylabel('Average Overlap per Segment (m)')
    ax_main.set_xlabel('Twist Angle ($^\circ$)')
    
    bs=np.array([200, 400, 1000])
    cs=[1]#np.array([0.8, 1])
    twist=np.linspace(1,90)
    markers=['r:', 'k:']
    close=pd.DataFrame()
    for b in bs:
        for i,c in enumerate(cs):
            o=b-b*c*np.cos(np.deg2rad(twist))
            l=ax_main.plot(twist,o*2,  markers[i], label="b="+str(b))
            labellines.labelLine(
                l[0],
                4,
                label="$b=$"+str(b)+" m",
                ha="left",
                va="bottom",
                backgroundcolor="none",
                fontsize=8.5
            )
            o2=b*(1-np.cos(np.deg2rad(df['TwistAngle'].values)))
            o_calc=df['b'].values*(1-np.cos(np.deg2rad(df['TwistAngle'].values)))
            o_close= (abs(o2-o_calc)/o2 < tol)
            print(df.loc[o_close])
            close=pd.concat((close,df.loc[o_close]), ignore_index=True)
            
    labelSubplots([ax_main, ax_xDist, ax_yDist], fontsize=12, labels=['d','e','f'])
    
    
    #ax.set_xscale('log')

    # calculated_overlap=dfAll['SegmentLSum'].values/dfAll['Size']-dfAll['StdRho'].values/2
    # fig, ax=plt.subplots()
    # sns.scatterplot(data=dfAll, x=dfAll['Overlap'].values*dfAll['SegmentLSum'].values, y=calculated_overlap, hue='TwistAngle', ax=ax)
    # ax.plot(calculated_overlap,calculated_overlap, 'r-.', label='1:1')
    # ax.set_xscale('log')
    # ax.set_yscale('log')
    # ax.set_ylim((min(dfAll['Overlap'].values*dfAll['SegmentLSum'].values), max(dfAll['Overlap'].values*dfAll['SegmentLSum'].values)))
    # ax.set_xlim((min(dfAll['Overlap'].values*dfAll['SegmentLSum'].values), max(dfAll['Overlap'].values*dfAll['SegmentLSum'].values)))
    # ax.legend()
    # ax.set_xlabel('Observed Overlap')
    # ax.set_ylabel('Calculated Overlap (Pollard 1987)')
    
    if return_EE:
        return fig, [ax_main, ax_xDist, ax_yDist], close
    else:
        return fig, [ax_main, ax_xDist, ax_yDist]

def PubWL_TwistOverlap(dfAllTrusted):
    fig1,ax=PubWidthsVsLengths(dfAllTrusted , hue='SwarmID', style='Region')
    fig1.savefig('/home/akh/myprojects/Linking-and-Clustering-Dikes/dikedata/WidthVsLength_trusted2.pdf', dpi=600)

    fig2,ax=PubTwistOverlap(dfAllTrusted)
    fig2.savefig("/home/akh/myprojects/Linking-and-Clustering-Dikes/dikedata/TwistOverlap.pdf", dpi=600)
    
    combinePlots(fig1, fig2, "/home/akh/myprojects/Linking-and-Clustering-Dikes/dikedata/WidthsLengthsTwistOverlap.pdf")

def ExampleClusters():
    df1=pd.read_csv('dikedata/crb/CJDS_FebStraightened.csv')
    df2=pd.read_csv('dikedata/deccandata/Central_preprocesed.csv')
    df3=pd.read_csv('dikedata/deccandata/NarmadaTapi_preprocesed.csv')
    
    SMALL_SIZE = 8
    MEDIUM_SIZE = 8
    BIGGER_SIZE = 12
    
    plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
    plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
    plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
    plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
    plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
    plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
    plt.rc('figure', titlesize=MEDIUM_SIZE)  # fontsize of the figure title
    plt.rcParams['legend.title_fontsize'] = SMALL_SIZE
    
    
    fig,ax=plt.subplots(1,3)
    fig.set_size_inches((190/25.4, 80/25.4))
    
    label1=428
    ax[0]=checkoutClusterCart(df1, label1, fig, ax[0])
    FixAxisAspect(ax[1],ax[0])
    ax[0].title.set_text('CRBG:CJDS - '+str(label1))
    
    label2=1896
    ax[1]=checkoutClusterCart(df2, label2, fig, ax[1])
    FixAxisAspect(ax[0],ax[1])
    ax[1].title.set_text('Deccan:Central - '+str(label2))
    
    label3=2608
    ax[2]=checkoutClusterCart(df3, label3, fig, ax[2])
    FixAxisAspect(ax[0],ax[2])
    ax[2].title.set_text('Deccan: N-T - '+str(label3))
    labelSubplots(ax, fontsize=12)
    
    plt.tight_layout()
    
    fig.savefig("/home/akh/myprojects/Linking-and-Clustering-Dikes/dikedata/ThreeExampleClusters.pdf", dpi=600)

def DilationPlots():
    dfs=['/home/akh/myprojects/Linking-and-Clustering-Dikes/dikedata/crb/allCRB_dikes_PreProcessed.csv', 
         '/home/akh/myprojects/Linking-and-Clustering-Dikes/dikedata/deccandata/AllDeccan_PreProcessed.csv']
         
    lines=['dikedata/crb/AllCRBLinked_euc_18_10_2022.csv',
           'dikedata/deccandata/AllDeccanLinked_euc_18_10_2022.csv']
    widths=[8,10]
    names=['CRB', 'Deccan']
    shapes=[ ['half', 'portrait'], ['half', 'landscape']]    
    for path1, path2, w, name,shape in zip(dfs,lines,widths,names, shapes):
        dikeset=pd.read_csv(path1)
        lines=pd.read_csv(path2)
        fig,axes=TripleDilationPlot(dikeset, lines, shape=shape, kwargs=dict(binWidth=1700, averageWidth=w))
        with open(name+"_limits.txt", 'w') as f:
            f.write('Y limits')
            f.write(str(axes[0].get_ylim()))
            f.write('X limits')
            f.write(str(axes[0].get_xlim()))
            
        print(axes[0].get_ylim())
        print(axes[0].get_xlim())
        fig.savefig('/home/akh/myprojects/Linking-and-Clustering-Dikes/dikedata/'+name+"TripleDilatonTest.pdf", dpi=600)
        writeToQGIS(lines,'/home/akh/myprojects/Linking-and-Clustering-Dikes/dikedata/'+name+"Qgis_Nov22.csv")
    
def PubTopLines():
    dfs=['/home/akh/myprojects/Linking-and-Clustering-Dikes/dikedata/crb/allCRB_dikes_PreProcessed.csv', 
         '/home/akh/myprojects/Linking-and-Clustering-Dikes/dikedata/deccandata/AllDeccan_PreProcessed.csv']
         
    lines=['dikedata/crb/AllCRBLinked_euc_18_10_2022.csv',
           'dikedata/deccandata/AllDeccanLinked_euc_18_10_2022.csv']
    
    SMALL_SIZE = 8
    MEDIUM_SIZE = 8
    BIGGER_SIZE = 12
    path='/home/akh/myprojects/Linking-and-Clustering-Dikes/dikedata/toplinesFigure/'
    plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
    plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
    plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
    plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
    plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
    plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
    plt.rc('figure', titlesize=MEDIUM_SIZE)  # fontsize of the figure title
    plt.rcParams['legend.title_fontsize'] = SMALL_SIZE
    plt.rcParams['text.usetex']=False
    names=['CRB', 'Deccan']
    tstep=7
    rstep=[]
    for path1, path2, name in zip(dfs,lines,names):
        print(name)
        dikeset=pd.read_csv(path1)
        lines=pd.read_csv(path2)
        rstep= np.ptp(lines['AvgRho'].values)*0.025 /1000# max(lines['Rho_Threshold'].unique())*10/1000
        #
        toplines, fig,axes, reds=TopHTSection(lines, dikeset, rstep, tstep, n=3)
        linesname=path+name+"TopLines_"+str(tstep)+str(int(rstep))+".csv"
        levels=np.zeros(len(lines))
        
        l=np.in1d(lines['Label'].values,toplines['Label'].values)
        #levels[l]==toplines['Level'].values
        lines=lines.assign(TopLines=l, TopHTLevels=levels)
        #writeToQGIS(toplines, linesname)
        plt.tight_layout()
        # with open(path+name+"_limits_toplines.txt", 'w') as f:
        #     f.write('Y limits')
        #     f.write(str(axes[1].get_ylim()))
        #     f.write("\n")
        #     f.write('X limits')
    
        #     f.write(str(axes[1].get_xlim()))
        #     f.write("\n")
        #     for i in range(5):
        #         f.write('level:'+str(i)+"  "+str(reds(i)))
        #         f.write("\n")
                
        #     f.write('level:'+str(-1)+"  g")
        #     f.write("\n")
        #fig.savefig('/home/akh/myprojects/Linking-and-Clustering-Dikes/dikedata/toplinesFigure/'+name+"toplines_rstep"+str(int(rstep))+".pdf", dpi=600)
        
def PubRadialSwarms():
    #fig = SetupJGRFig((230,90), 'landscape')
    gskw = dict(height_ratios= [ .30,.70])
    #fig2=SetupJGRFig('full', 'landscape')
    
    gs = gridspec.GridSpec(2, 4, **gskw)
    # htCRBG1 = plt.subplot(gs[0,0])
    # htDeccan1= plt.subplot(gs[0,2])
    
    # htCRBG = plt.subplot(gs[0,1])
    # htDeccan= plt.subplot(gs[0,3])
    fig, [htCRBG1, htCRBG, htDeccan1, htDeccan]=plt.subplots(1,4)
    fig.set_size_inches( (230/25.4, 95/25.4))
    fig2,CRBG=plt.subplots()#plt.subplot(gs[1,0])
    fig3,Deccan=plt.subplots()#plt.subplot(gs[1,1:])
    
    theta=np.linspace(-90,90)
    #CRBG 
    dikeset=pd.read_csv('/home/akh/myprojects/Linking-and-Clustering-Dikes/dikedata/crb/allCRB_dikes_PreProcessed.csv')
    lines=pd.read_csv('dikedata/crb/AllCRBLinked_euc_18_10_2022.csv')
    xc,yc=HT_center(dikeset)
    DotsHT(fig, htCRBG1, lines,  ColorBy='Ymid', cmap='magma')
    DotsHT(fig, htCRBG, lines, alpha=0.05, ColorBy=None, axlabels=(True,False))
    plotlines(lines, 'k', CRBG, alpha=0.2)
    

    #all-black

    All=RadialFit(lines)
    AllClose, All=NearCenters(lines, All, tol=100000)
    rhoAll=CenterFunc(theta,All['Center'][0][0], All['Center'][0][1], xc,yc)/1000
    htCRBG.plot(theta,rhoAll, 'k-.', linewidth=2)
    
    #htCRBG.axhline(y=dikeset['rho'].std()/1000*2)
    #htCRBG.axhline(y=dikeset['rho'].std()/-1000*2)
    htCRBG.tick_params(labelleft=False)   
    
    # Red 
    mask=lines['TrustFilter']==1
    mlines=lines[mask]
    mask2=dikeset['Ymid']<5e6
    Red=RadialFit(dikeset[mask2])
    RedClose,Red=NearCenters(lines,Red)
    plotlines(RedClose, 'r', CRBG)
    
    rhoRed=CenterFunc(theta,Red['Center'][0][0], Red['Center'][0][1], xc,yc)/1000
    htCRBG.plot(theta,rhoRed, 'r', linewidth=2)
    DotsHT(fig, htCRBG, RedClose, color='r', ColorBy=None, axlabels=(True,False))
    
    #Green
    mask2=dikeset['Ymid']>5.0e6
    Green=RadialFit(dikeset[mask2])
    GreenClose, Green=NearCenters(lines, Green)
    rhoGreen=CenterFunc(theta,Green['Center'][0][0], Green['Center'][0][1], xc,yc)/1000
    htCRBG.plot(theta,rhoGreen, 'g', linewidth=2)
    DotsHT(fig, htCRBG, GreenClose, color='g', ColorBy=None, axlabels=(True,False))
    
    plotlines(GreenClose, 'g', CRBG)
    
    CRBG.plot( Green['Center'][0][0], Green['Center'][0][1], '*g', markersize=10, markeredgecolor='k')
    CRBG.plot( Red['Center'][0][0], Red['Center'][0][1], '*r', markersize=10, markeredgecolor='k')
    CRBG.plot(  All['Center'][0][0], All['Center'][0][1], '*k', markersize=10, markeredgecolor='white')
    

    CRBG.plot(xc,yc, "wX", markeredgecolor="black", markersize=10)
    # circle=plt.Circle((xc,yc), 2*dikeset['rho'].std(), color='grey', fill=False, linestyle='--')
    # CRBG.add_patch(circle)
    CRBGCenters=pd.concat((All, Red, Green))
    CRBGCloseLines=pd.concat((AllClose, RedClose, GreenClose))
    writeToQGIS(CRBGCloseLines, 'dikedata/RadialFits/CRBGRadialFitsLines.csv')
    writeCenterWKT(CRBGCenters, 'dikedata/RadialFits/CRBGRadialFitsCenters.csv')
    
    
    #Deccan
    
    lines=pd.read_csv('dikedata/deccandata/AllDeccanLinked_euc_18_10_2022.csv')
    dikeset=pd.read_csv('/home/akh/myprojects/Linking-and-Clustering-Dikes/dikedata/deccandata/AllDeccan_PreProcessed.csv')
    xc,yc=HT_center(dikeset)
    DotsHT(fig, htDeccan1, lines, ColorBy='Ymid', cmap='magma')
    DotsHT(fig, htDeccan, lines, alpha=0.05, ColorBy=None, axlabels=(True,False))
    #htDeccan.axhline(y=dikeset['rho'].std()/1000*2)
    #htDeccan.axhline(y=dikeset['rho'].std()/-1000*2)
    htDeccan.tick_params(labelleft=False)  
    plotlines(lines, 'k', Deccan, alpha=0.2)

    #all-black

    All=RadialFit(lines)
    AllClose, All=NearCenters(lines, All, tol=100000)
    rhoAll=CenterFunc(theta,All['Center'][0][0], All['Center'][0][1], xc,yc)/1000
    htDeccan.plot(theta,rhoAll, 'k-.', linewidth=2)

    #Red
    mask=lines['TrustFilter']==1
    mlines=lines[mask]
    mask2=mlines['Ymid']<2.3e6
    Red=RadialFit(mlines[mask2])
    RedClose,Red=NearCenters(lines,Red)
    plotlines(RedClose, 'r', Deccan)
    
    rhoRed=CenterFunc(theta,Red['Center'][0][0], Red['Center'][0][1], xc,yc)/1000
    htDeccan.plot(theta,rhoRed, 'r', linewidth=2)
    DotsHT(fig, htDeccan, RedClose, color='r', ColorBy=None, axlabels=(True,False))
    
    #Green
    mask2=(mlines['Ymid']>2.38e6) & (mlines['Ymid']<2.42e6) 
    Green=RadialFit(mlines[mask2])
    GreenClose, Green=NearCenters(lines, Green)
    rhoGreen=CenterFunc(theta,Green['Center'][0][0], Green['Center'][0][1], xc,yc)/1000
    htDeccan.plot(theta,rhoGreen, 'g', linewidth=2)
    DotsHT(fig, htDeccan, GreenClose, color='g', ColorBy=None, axlabels=(True,False))
    
    plotlines(GreenClose, 'g', Deccan)
    

    #purple
    mask2=(dikeset['Ymid']>2.64e6) 
    Purp=RadialFit(dikeset[mask2])
    PurpClose, Purp=NearCenters(lines, Purp)
    #PurpCloseSeg, Purp=NearCenters(dikeset, Purp)
    rhoPurp=CenterFunc(theta,Purp['Center'][0][0], Purp['Center'][0][1], xc,yc)/1000
    htDeccan.plot(theta,rhoPurp, color='purple', linewidth=2)
#    DotsHT(fig, htDeccan, PurpCloseSeg, color='purple', ColorBy=None, axlabels=(True,False))
    plotlines(PurpClose, 'purple', Deccan)
    
    labelSubplots([htCRBG1, htCRBG, htDeccan1, htDeccan])
    plt.tight_layout()
    

    #yellow
    mask2=(dikeset['Ymid']<1.9e6) 
    Yellow=RadialFit(dikeset[mask2])
    YellowClose, Yellow=NearCenters(lines, Yellow)
    #YellowCloseSeg, Yellow=NearCenters(dikeset, Yellow)
    #YellowClose=YellowClose.assign("Color")
    rhoYellow=CenterFunc(theta,Yellow['Center'][0][0], Yellow['Center'][0][1], xc,yc)/1000
    htDeccan.plot(theta,rhoYellow, color='y', linewidth=2)
    #DotsHT(fig, htDeccan, YellowCloseSeg, color='yellow', ColorBy=None, axlabels=(True,False))
    #plotlines(YellowClose, 'yellow', Deccan)
    
    Deccan.plot(  Yellow['Center'][0][0], Yellow['Center'][0][1], '*', color='yellow', markersize=10, markeredgecolor='k')

    Deccan.plot(  Purp['Center'][0][0], Purp['Center'][0][1], '*', color='purple', markersize=10, markeredgecolor='k')
    Deccan.plot(  Green['Center'][0][0], Green['Center'][0][1], '*g', markersize=10, markeredgecolor='k')
    Deccan.plot(  Red['Center'][0][0], Red['Center'][0][1], '*r', markersize=10, markeredgecolor='k')
    Deccan.plot(  All['Center'][0][0], All['Center'][0][1], '*', color='black', markersize=10, markeredgecolor='k')
    Deccan.plot(xc,yc, "wX", markeredgecolor="black", markersize=10)
    
    # circle=plt.Circle((xc,yc), 2*dikeset['rho'].std(), color='grey', fill=False, linestyle='--')
    # Deccan.add_patch(circle)
    labelSubplots([htCRBG1, htCRBG, htDeccan1, htDeccan])

    
    

    DeccanCenters=pd.concat((All, Purp, Red, Green, Yellow))
    DeccanCloseLines=pd.concat((AllClose, PurpClose, RedClose, GreenClose, YellowClose))
    writeToQGIS(DeccanCloseLines, 'dikedata/RadialFits/DeccanRadialFitsLines.csv')
    writeCenterWKT(DeccanCenters, 'dikedata/RadialFits/DeccanRadialFitsCenters.csv')
    fig.tight_layout()
    #fig.savefig('dikedata/RadialFits/PubRadialFits.pdf', dpi=600)
    
    #fig.savefig('dikedata/RadialFits/PubRadialFits.png', dpi=900)
    
def PubSPeaksRadial():
    
    # fig, [ht1, ht2]=plt.subplots(1,2)
    # fig.set_size_inches( (115/25.4, 95/25.4))
    fig2,Cart=plt.subplots()
    
    fig = SetupJGRFig((90,170), 'landscape')
    dist_width=0.10
    gskw = dict(height_ratios= [ dist_width,1-dist_width], width_ratios=[(1-dist_width)/2,(1-dist_width)/2,dist_width])
    #fig2=SetupJGRFig('full', 'landscape')
    
    gs = gridspec.GridSpec(2, 3, **gskw)
    ht1 = plt.subplot(gs[:,0])
    ht2= plt.subplot(gs[1,1])
    tdist=plt.subplot(gs[0,1])
    rdist=plt.subplot(gs[1,2])
    
    # htCRBG = plt.subplot(gs[0,1])
    # htDeccan= plt.subplot(gs[0,3])
    
    theta=np.linspace(-90,90)
   
    dikeset=pd.read_csv('/home/akh/myprojects/Linking-and-Clustering-Dikes/dikedata/spanish peaks/SpanishPeaks_3857_preprocessed.csv')
    lines=pd.read_csv('/home/akh/myprojects/Linking-and-Clustering-Dikes/dikedata/spanish peaks/SpanishPeaks_Complete_euc_3_2013.csv')
    xc,yc=HT_center(dikeset)
    DotsHT(fig, ht1, lines,  ColorBy='Ymid', cmap='magma')
    DotsHT(fig, ht2, lines, alpha=0.05, ColorBy=None, axlabels=(True,False))
    plotlines(lines, 'k', Cart, alpha=0.2)
    

    #purple
    mask2=(dikeset['Ymid']>4.52e6) 
    Purp=RadialFit(dikeset[mask2])
    PurpClose, Purp=NearCenters(lines, Purp, tol=2500)
    PurpCloseSeg, Purp=NearCenters(dikeset, Purp, tol=2500)
    rhoPurp=CenterFunc(theta,Purp['Center'][0][0], Purp['Center'][0][1], xc,yc)/1000
    ht2.plot(theta,rhoPurp, color='purple', linewidth=2)
    DotsHT(fig, ht2, PurpCloseSeg, color='purple', ColorBy=None, axlabels=(True,False))
    plotlines(PurpClose, 'purple', Cart)
    
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
    
    Centers=pd.concat((green, Purp))
    
    CloseLines=pd.concat((GreenClose, PurpClose, linl))
    CloseSegments=pd.concat((GreenCloseSeg, PurpCloseSeg))
    writeCenterWKT(Centers, 'dikedata/RadialFits/SpanishPeaksRadialFitsCenters.csv')
    writeToQGIS(CloseLines, 'dikedata/RadialFits/SPRadialFitsLines.csv')
    writeToQGIS(CloseSegments, 'dikedata/RadialFits/SPRadialFitsSegments.csv')

    #fig.savefig('dikedata/RadialFits/SPRadialFitsLinesFig.pdf', dpi=600)
    #fig.savefig('dikedata/RadialFits/SPRadialFitsLinesFig.png', dpi=600)
    
def RadialAnalysis():
    linestyles=["-", "--", "-.", ":", (0, (1, 10))]
    sp=pd.read_csv('/home/akh/myprojects/Linking-and-Clustering-Dikes/dikedata/RadialFits/SpanishPeaksRadialFitsCenters.csv')
    dec=pd.read_csv('/home/akh/myprojects/Linking-and-Clustering-Dikes/dikedata/RadialFits/DeccanRadialFitsCenters.csv')
    crbg=pd.read_csv('/home/akh/myprojects/Linking-and-Clustering-Dikes/dikedata/RadialFits/CRBGRadialFitsCenters.csv')
    sp=sp.assign(ls=linestyles[:len(sp)], label=["SP"+str(i) for i in np.arange(len(sp))])
    sp=sp.assign(Colors=[ 'green', 'purple'], name="Spanish Peaks")
    dec=dec.assign(ls=linestyles[:len(dec)], label=["Dec"+str(i) for i in np.arange(len(dec))])
    dec=dec.assign(Colors=['black', 'purple', 'red', 'green', 'yellow'] , name="Deccan")
    
    crbg=crbg.assign(ls=linestyles[:len(crbg)], label=["CRBG"+str(i) for i in np.arange(len(crbg))])
    crbg=crbg.assign(Colors=['black', 'red', 'green'] , name="CRBG")
    
    fig = SetupJGRFig((90,190), 'landscape')
    
    ax=np.array([fig.add_subplot(121), fig.add_subplot(122)])
    Rtols=np.array([0.0001,0.005, 0.01, 0.1, 0.2, 0.5, 0.75, 1])
    convert1= lambda x: np.fromstring(x[1:-1], sep=" ")
    convert2= lambda x: np.fromstring(x[1:-1], sep=",")
    n=36
    l=n/360
    
    K_Pure=np.ones(35)*l/n
    #ax[1].plot(np.arange(0,350,10), K_Pure, 'k-.')
    
    for i, ls in zip([sp,dec,crbg], ["-", "--", (0, (1, 10))]):
        
        for ind,j in i.iterrows():
            c=np.floor(convert1(j['Center']))
            label=j['name']+": ("+str(int(c[0]))+"E, "+str(int(c[1]))+"N)"
            ax[0].plot(Rtols, convert2(j['ExpandingR']), color=j['Colors'], ls=ls, label=label)# j['label'])
            k=convert1(j['AllK'])-np.arange(0,350,10)
            ax[1].plot(np.arange(0,350,10), k, color=j['Colors'], ls=ls, label=label) # j['label'])
    ax[0].legend()
    ax[0].yaxis.set_major_formatter(PercentFormatter(xmax=1))
    
    l=10000
    df2=makeLinear2(l, 75, 2, 0, 1000, ndikes=150, label=1)
    df2=df2.append(makeLinear2(l, 30, 2, 0, 700, ndikes=150, label=2))
    df2=df2.append(makeLinear2(l, -30, 2, 0, 500, ndikes=150, label=3))
    for i in ['Xstart', 'Ystart', 'Xend', 'Yend']:
        df2[i]=df2[i].values+10000
    #df2['Label']=df2['Label'].astype(str)
    df2=DikesetReProcess(df2)
    xdist=df2['Xmid'].values
    ydist=df2['Ymid'].values
    
    rAngle=np.rad2deg(np.arctan2(xdist,ydist))+180
    
    ripley1,_=RipleyRadial(rAngle)
    
    l=5000
    df3=makeRadialSwarmdf(l, center=[0,0], ndikes=5)
    df3=DikesetReProcess(df3)
    xdist=df3['Xmid'].values
    ydist=df3['Ymid'].values
    
    rAngle=np.rad2deg(np.arctan2(xdist,ydist))+180
    ripley2,_=RipleyRadial(rAngle)
    ax[1].plot(np.arange(0,350,10), ripley1-np.arange(0,350,10), color='grey', linewidth=5, label=label)
    ax[1].plot(np.arange(0,350,10), ripley2-np.arange(0,350,10), color='grey', linewidth=3, label=label)
    
    
    ax[0].set_xscale('log')
    ax[0].set_ylabel('Percentage of Swarm Intersecting within Center')
    ax[0].set_xlabel('Percentage of Swarm Total Radius')
    ax[0].xaxis.set_major_formatter(PercentFormatter(xmax=1, decimals=2))
    
    ax[1].set_ylabel('Ripley K (k-s)')
    ax[1].set_xlabel('Angular Distance ($^\circ$)')
    plt.tight_layout()
    fig.savefig('/home/akh/myprojects/Linking-and-Clustering-Dikes/Publication Figures/SupRipleyK.pdf')
    
def comparisonFig():
    
    averagewidth=[10,8,2]
    name=['Deccan', 'CJDS', 'West Peak']
    EruptedVolume=np.array([ 1.3e6, 210000, 100])
    volScaled=(EruptedVolume/np.min(EruptedVolume))**(1/3)
    
    figV = SetupJGRFig((60,90), 'portrait')

    vol = figV.add_subplot(111, projection='3d')
    
    
    offset=20
    x, y, z = np.indices((46, 46, int(np.sum(volScaled)+2*offset+10)))
    cube1=(x > 23-volScaled[2]/2) & (y > 23-volScaled[2]/2) & (z > volScaled[0]+2*offset+volScaled[1]) & (x < 23+volScaled[2]/2) & (y < 23+volScaled[2]/2) & (z < volScaled[0]+2*offset+volScaled[1]+volScaled[2])
    cube2=(x > 23-volScaled[1]/2) & (y > 23-volScaled[1]/2) & (z > volScaled[0]+offset) & (x < 23+volScaled[1]/2) & (y < 23+volScaled[1]/2) & (z < volScaled[0]+offset+volScaled[1])
    cube3=(x > 23-volScaled[0]/2) & (y > 23-volScaled[0]/2) & (z > 1) & (x < 23+volScaled[0]/2) & (y < 23+volScaled[0]/2) & (z < volScaled[0])
    
    voxelarray= cube1 | cube2 | cube3
    colors = np.empty(voxelarray.shape, dtype=object)
    colors[cube1] = 'red'
    colors[cube2] = 'blue'
    colors[cube3] = 'green'
    
    # and plot everything
    vol.voxels(voxelarray, facecolors=colors, alpha=0.4)
    vol.set_box_aspect([1,1,1])
    vol.xaxis.set_ticklabels([])
    vol.yaxis.set_ticklabels([])
    vol.zaxis.set_ticklabels([])

    fig = SetupJGRFig((90,90), 'landscape')
    area = fig.add_subplot(131)
    length = fig.add_subplot(132)
    
    
    fig2 = SetupJGRFig((65,90), 'portrait')
    axSeg=fig2.add_subplot(111)
    
    axArea=axSeg.twinx()
    #axSeg=np.array([fig2.add_subplot(311), fig2.add_subplot(312), fig2.add_subplot(313) ])

    d=['dikedata/deccandata/AllDeccan_PreProcessed.csv',
        'dikedata/crb/allCRB_dikes_PreProcessed.csv',
        'dikedata/RadialFits/SPRadialFitsSegments.csv'
        ]
    
    c=['green', 'blue', 'red']
    
        
    dfAll=pd.read_csv('/home/akh/myprojects/Linking-and-Clustering-Dikes/dikedata/AllDatasetsLinked_18_10_2022.csv')
    dfAll=dfAll.assign(R_Width=dfAll['R_Width']+8)
    dfAll['Dike Cluster Length (km)']=dfAll['R_Length']/1000
    dfAll['Dike Cluster Width (m)']=dfAll['R_Width']+2
    
    SpRadial=pd.read_csv('/home/akh/myprojects/Linking-and-Clustering-Dikes/dikedata/RadialFits/SPRadialFitsLines.csv')
    CRBGlines=pd.read_csv('dikedata/crb/AllCRBLinked_euc_18_10_2022.csv')
    declines=pd.read_csv('dikedata/deccandata/AllDeccanLinked_euc_18_10_2022.csv')
    nlines=[len(declines), len(CRBGlines),np.sum(SpRadial['Structure']=='Radial 1')]
    
    dfAll=dfAll.assign(Status=[ 'Filtered Cluster' if i==1 else 'Cluster' for i in dfAll['TrustFilter']])
    dfAllTrusted=dfAll[ dfAll['TrustFilter']==1]
    dfAll=dfAll[dfAll['Linked']==1]

    rw=[np.median(dfAllTrusted['R_Width'].loc[dfAllTrusted['Region']==j]) for j in ['Deccan', 'CRBG', 'SpanishPeaks']]
    rl=[np.median(dfAllTrusted['R_Length'].loc[dfAllTrusted['Region']==j]) for j in ['Deccan', 'CRBG', 'SpanishPeaks']]
    offset=0
    offset2=0
    offset3=0
    areas=[]
    aspects=[]
    dikes=[ 29000,4000,np.sum(SpRadial[ SpRadial['Structure']=='Radial 1']['Size'])]
    scaleCRBG= lambda x: x/x[1]
    scaleSP= lambda x:x/x[2]
    
    dikes=scaleSP(dikes)*2
    ndikes=[]
    i=0
    for dikeset,color, label,w,l, dikes, volume in zip(d,c, name, rw, rl, nlines, EruptedVolume): 
        
        df=pd.read_csv(dikeset)
        if label=='West Peak':
            df=df[df['Structure']=='Radial 1']
        #df=DikesetReProcess(df)
        ndikes.append(len(df))
        xc,yc=HT_center(df)
        xs=[ min(df['Xstart'].min(), df['Xend'].min()), max(df['Xstart'].max(), df['Xend'].max())]
        ys=[ min(df['Ystart'].min(), df['Yend'].min()), max(df['Ystart'].max(), df['Yend'].max())]
        width=np.min([np.ptp(xs),np.ptp(ys)])
        
        height=np.max([np.ptp(xs),np.ptp(ys)])
        areas.append(width*height)
        
        Xedges=np.array([0, 0, width, width, 0])
        Yedges=np.array([height, 0, 0, height, height])
        area.plot(Xedges, Yedges, color, linewidth=2, alpha=0.5)
        rect=Rectangle( (0, 0), width, height ,facecolor=color, alpha=0.3, label=label)
        area.add_patch(rect)
    
        seglength=np.median(df['seg_length'])
        
        # Xedges=np.array([offset, offset, offset+w, offset+w, offset])
        # Yedges=np.array([l, 0, 0, l, l])
        
        # theta=np.random.randint(1,90, int(nsegs))
        # a=np.array([1,-1])
        # c=np.random.choice(a, int(nsegs))
        # theta=theta*c
        
        # rho=np.random.randint(5,20, int(nsegs))
        # segDf=fromHT(theta,rho, length=10, scale=10, xc=0, yc=offset3)
        
        #plotlines(segDf, color, seg, equal=False)
        length.plot( Xedges, Yedges, color, linewidth=2, alpha=0.5, label=label)
        
        
        # rect1=Rectangle( (0, 0), 50, len(df)  ,facecolor=color, alpha=0.3, label=label)
        # rect2=Rectangle( (100, 0), 50, dikes,facecolor=color, alpha=0.3, label=label)
        # axSeg[i].add_patch(rect1)
        # axSeg[i].add_patch(rect2)
        # axSeg[i].autoscale_view()
        # axSeg[i].set_ylim((0,25000))
        
        axSeg.scatter( dikes, volume, color=color, marker="p")
        axSeg.scatter( len(df), volume, color=color, marker="*")
        axArea.scatter( dikes,width*height, color=color, marker="^")
        axArea.scatter( len(df),width*height, color=color, marker="P")
        
        #seg.plot( [offset2, offset2], [0,seglength], color=color )
        
        offset=offset+w+100
        offset2=offset2+10
        offset3=offset3+50
        aspects.append(l/w)
        
        i=i+1
    

    areas=np.array(areas)
    aspects=np.array(aspects)
    ndikes=np.array(ndikes)
    nlines=np.array(nlines)
        
    axSeg, l1= plotRatioLine(axSeg, ndikes, 50, line_kw=dict(color='k', ls=":", label="50 $km^3$ per Segment"))
    axSeg, l2=plotRatioLine(axSeg, nlines, 200, line_kw=dict(color='k', ls=":", label="200 $km^3$ per Cluster"))

    axArea, l3=plotRatioLine(axArea, nlines, np.mean(areas/nlines), line_kw=dict(color='orange', ls=":", label="$10^7 km^2 per Cluster$"))
    axArea, l4=plotRatioLine(axArea, ndikes,np.mean(areas/ndikes), line_kw=dict(color='orange', ls=":", label="$10^7 km^2 per Segment$"))
    
    # labellines.labelLine(
    #     l1[0],
    #     6000,
    #     label=r"50 $km^3$ per Segment",
    #     ha="left",
    #     va="center",

    #     backgroundcolor="none",
    #     fontsize=8
    # )
    
    # labellines.labelLine(
    #     l2[0],
    #     1500,
    #     label=r"200 $km^3$ per Cluster",
    #     ha="left",
    #     va="center",

    #     backgroundcolor="none",
    #     fontsize=8
    # )
    
    # labellines.labelLine(
    #     l3[0],
    #     500,
    #     label=r"$1.5x10^8 km^2 per Cluster$",
    #     ha="left",
    #     va="center",

    #     backgroundcolor="none",
    #     fontsize=8
    # )
    
    # labellines.labelLine(
    #     l4[0],
    #     1000,
    #     label=r"$4x10^7 km^2 per Segment$",
    #     ha="left",
    #     va="center",

    #     backgroundcolor="none",
    #     fontsize=8
    # )
    
    
    #axArea.set_xscale('log')
    axSeg.set_xscale('log')
    axSeg.set_yscale('log')
    axArea.set_yscale('log')

   # areas[2]=2607815002.00

    t=['Cluster Length', 'Cluster Width', 'Area', 'Aspect', 'Erupted Volume', 'Segments', 'Lines']
    print("Scaled by CRBG")
    for i,j in zip(t, [rl, rw, areas, aspects, EruptedVolume, ndikes, nlines]):
        print(i)
        print(scaleCRBG(j))
        
    print("Scaled by Spanish Peaks")
    for i,j in zip(t, [rl, rw, areas, aspects, EruptedVolume, ndikes, nlines]):
        print(i)
        print(scaleSP(j))

    area.xaxis.set_ticklabels([])
    area.yaxis.set_ticklabels([])
    
    length.xaxis.set_ticklabels([])
    length.yaxis.set_ticklabels([])
    
    seg.xaxis.set_ticklabels([])
    seg.yaxis.set_ticklabels([])
    area.legend()
    plt.tight_layout()
    
    figV.savefig('/home/akh/myprojects/Linking-and-Clustering-Dikes/Publication Figures/comparisonPart1.pdf', dpi=600)
    fig.savefig('/home/akh/myprojects/Linking-and-Clustering-Dikes/Publication Figures/comparisonPart2.pdf', dpi=600)
    fig2.savefig('/home/akh/myprojects/Linking-and-Clustering-Dikes/Publication Figures/comparisonPart3.pdf', dpi=600)
    
def RadialIdentification():
    
    # how far from the HT origin is a Radial Swarm Identifiable
    SMALL_SIZE = 8
    MEDIUM_SIZE = 8
    BIGGER_SIZE = 12
    
    plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
    plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
    plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
    plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
    plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
    plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
    plt.rc('figure', titlesize=MEDIUM_SIZE)  # fontsize of the figure title
    plt.rcParams['legend.title_fontsize'] = SMALL_SIZE
    
    angles=np.random.randint(0,180, 100)-90
    #angles[abs(angles)>90]=angles[abs(angles)>90]%90*np.sign(angles[abs(angles)>90])
    rhos=np.random.normal(0,1000, 100)
    df=fromHT(angles, rhos)
    fig, ax=plt.subplots(2,1)
    w=95/25.4 ## cm/(cm/in)
    l=115/25.4
    fig.set_size_inches(l,w)
    
    dist=np.array([0,1,2.5,10,30,50,100,200,300,1000])*1000
    r1=[]
    r10=[]
    r1000=[]
    
    for d in dist:
        d=d/np.sqrt(2)
        df=fromHT(angles, rhos, xc=d, yc=d, scale=100000, CartRange=10000000)
        
        Center=RadialFit(df, xc=0, yc=0)

        r1.append(Center['RSq'])
    
        
        df=fromHT(angles, rhos*10, xc=d, yc=d, scale=100000, CartRange=10000000)
        Center=RadialFit(df, xc=0, yc=0)

        r10.append(Center['RSq'])
        
        df=fromHT(angles, rhos*100, xc=d, yc=d, scale=100000, CartRange=10000000)
        Center=RadialFit(df, xc=0, yc=0)
        r1000.append(Center['RSq'])


    ax[0].plot(dist/1000, r1, 'r^-', label="$\mu$=1 km")
    ax[0].plot(dist/1000, r10, 'bp-', label="$\mu$=10 km")
    ax[0].plot(dist/1000, r1000, 'g*-', label="$\mu$=100 km")
    ax[0].set_xlabel('Distance from Origin (km)')
    
    ax[1].plot(dist/np.std(rhos), r1, 'r^-', label="$\mu$=1 km")
    ax[1].plot(dist/np.std(rhos*10), r10, 'bp-', label="$\mu$=10 km")
    ax[1].plot(dist/np.std(rhos*100), r1000, 'g*-', label="$\mu$=100 km")
    ax[0].legend()
    ax[1].set_xlabel('Distance from Origin / $\mu$')
    ax[1].set_xlim(0,15)
    
    for a in ax: 
        a.set_ylabel('$R_{sq}$')
    labelSubplots(ax)
    plt.tight_layout()
    
    fig.savefig('Publication Figures/RadialIdeniticationHTOrigin.pdf', dpi=600)
