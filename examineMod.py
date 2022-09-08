#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  1 12:54:02 2021

@author: akh
"""
import numpy as np 
import pandas as pd 
from fitRectangle import *
from htMOD import HT_center
from plotmod import HThist, plotlines, DotsHT, clustered_lines, pltRec, FixCartesianLabels
from sklearn.neighbors import NearestNeighbors
from sklearn.preprocessing import scale
from PrePostProcess import *
from htMOD import rotateData2
from clusterMod import *

from PrePostProcess import completePreProcess, whichForm, midPoint
import scipy.cluster.hierarchy as sch
from scipy.spatial.distance import squareform
from scipy.cluster.hierarchy import dendrogram
from scipy.spatial.distance import pdist, squareform

import matplotlib.ticker as mticker
import matplotlib.pyplot as plt 
from scipy import stats
import warnings
from matplotlib import cm
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import os 
from datetime import datetime

def checkoutCluster(dikeset, label):
    
    xc=dikeset['xc'].iloc[0]
    yc=dikeset['yc'].iloc[0]
    
    fig, ax=plt.subplots(1,2) 
    fig.set_size_inches( 12,6)
    fig.suptitle("Label "+str(int(label)))
    mask=(dikeset['Labels']==label)
    lines=dikeset[mask]
    
    xlim1=min(lines['theta'])-2
    xlim2=max(lines['theta'])+2
    ylim1=min(lines['rho'])-10000
    ylim2=max(lines['rho'])+10000
    
    ax[0].scatter(dikeset['theta'], dikeset['rho'], c='grey', alpha=0.2) 
    ax[0].scatter(lines['theta'], lines['rho'], c='r')
    ax[0].plot(lines['theta'].mean(), lines['rho'].mean(), 'g*')
    
    ax[0].plot([xlim1, xlim1, xlim2, xlim2, xlim1], [ylim1,ylim2, ylim2, ylim1, ylim1], 'k')
    
    if np.mean(lines['rho'].values) > 0 :
        left, bottom, width, height = [0.2, 0.2, 0.2, 0.2]
    else:
        left, bottom, width, height = [0.2, 0.55, 0.2, 0.2]
    
    ax2 = fig.add_axes([left, bottom, width, height])
    ax2.scatter(lines['theta'], lines['rho'], c='r')
    ax2.plot(lines['theta'].mean(), lines['rho'].mean(), 'g*')
    ax2.scatter(dikeset['theta'], dikeset['rho'], c='grey', alpha=0.1) 
    ax2.set_ylim([ylim1+9500, ylim2-9500])
    ax2.set_xlim([xlim1+1.5, xlim2-1.5])
    
    ax[0].set_xlabel('Theta ($^\circ$)')
    ax[0].set_ylabel('Rho (m)')
    
    fig,ax[1], l,w=pltRec(lines, xc,yc, fig=fig, ax=ax[1])

    
    if np.sign(lines['theta'].mean()) < 0:
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
    
    if lines['theta'].mean()> 0:
        ax4.yaxis.tick_right()
        ax4.yaxis.set_label_position('right') 
        

    FixCartesianLabels(ax4)
    FixCartesianLabels(ax[1])
    
    " Need to fix the aspect ratios now"
    
    figW0, figH0 = ax[0].get_figure().get_size_inches()
    # Axis size on figure
    _, _, w0, h0 = ax[0].get_position().bounds
    # Ratio of display units
    disp_ratio0 = (figH0 * h0) / (figW0 * w0)
    
    figW1, figH1 = ax[1].get_figure().get_size_inches()
    # Axis size on figure
    _, _, w1, h1 = ax[1].get_position().bounds
    # Ratio of display units
    disp_ratio1 = (figH1 * h1) / (figW1 * w1)
    
    if h0 > h1:
        ratio=h0/h1
        y_min, y_max = ax[1].get_ylim()
        deltay=(y_max-y_min)*ratio/2
        ax[1].set_ylim(y_min-deltay, y_max+deltay)
        print('reset Y')
        
    if w0>w1:
        ratio=w0/w1
        y_min, y_max = ax[1].get_xlim()
        deltay=(y_max-y_min)*ratio/2
        ax[1].set_xlim(y_min-deltay, y_max+deltay)
        print('reset X')
    
    if np.sign(lines['rho'].mean())>0:
        tlocy=0.95
    else:
        tlocy=0.25
        

    t=ax[0].text(tlocx2,tlocy-0.10,'Angle Mean: '+str( round(lines['theta'].mean(), 2)), transform=ax[0].transAxes)

    t=ax[0].text(tlocx2,tlocy-0.05,'Rho Mean: '+str(round(lines['rho'].mean(),1)), transform=ax[0].transAxes)

        
    t=ax[0].text(tlocx2,tlocy-0.20,'Angle Range: '+str( round(np.ptp(lines['theta']), 2)), transform=ax[0].transAxes)

    t=ax[0].text(tlocx2,tlocy-0.15,'Rho Range: '+str(round(np.ptp(lines['rho']),1)), transform=ax[0].transAxes)


    tlocy=.95
    ax[1].text(tlocx,tlocy,'Length: '+str(round(l,1)), transform=ax[1].transAxes)
    ax[1].text(tlocx,tlocy-0.05,'Width: '+str(round(w,1)), transform=ax[1].transAxes)
    ax[1].text(tlocx,tlocy-0.10,'Aspect: '+str(round(l/w,1)), transform=ax[1].transAxes)
    ax[1].text(tlocx,tlocy-0.15,'Size: '+str(len(lines)), transform=ax[1].transAxes)
    plt.tight_layout()
    
    print('Angle Mean: '+str( round(lines['theta'].mean(), 2)))
    print('Rho Mean: '+str(round(lines['rho'].mean(),1)))
    print('Length: '+str(round(l,1)))
    print('Width: '+str(round(w,1)))
    
    


    return fig,ax

# def overlap(min1, max1, min2, max2):
#     return max(0, min(max1, max2) - max(min1, min2))

def CheckoutBy(dikeset, lines, col, maximum=True, minimum=False):
    
    loc=np.nanargmax(lines[col].values)
    label=lines['Label'].iloc[loc]
    print("Label", int(label))
    
    fig,ax=checkoutCluster(dikeset, label)
    
    if minimum:
        loc=np.where(lines[col].values==np.min(lines[col].values))
        label=lines['Label'].loc[loc]
        fig,ax=checkoutCluster(dikeset, label)
        
    return fig,ax


def RotateOverlap(lines):
    theta=np.mean(lines['theta'].values)
    
    dfRotated=rotateData2(lines, (90-theta))
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
        avgtheta=np.average((lines['theta']))
        w,l,r,Xe, Ye, Xmid, Ymid=fit_Rec(lines, xc, yc)

        x1, x2, y1, y2=clustered_lines(x,y,avgtheta, l, xmid=Xmid, ymid=Ymid)

        
        hashlines=hash((lines['HashID'].values.tostring()))
        clusters_data=clusters_data.append({ "Label": i, "Xstart": x1, "Ystart":y1, "Xend": x2,
                                            "Yend":y2, "X0": x0, "Y0": y0, "AvgRho":avgrho,
                                            "AvgTheta":avgtheta, "ClusterHash": hashlines, "Size":size}, ignore_index=True)
    return clusters_data

                                            
def examineClusters(clusters, enEchelonCutofff=7, ifEE=False, MaxNNSegDist=0.5, skipUnlinked=False):
    """ Must have  ['Xstart', 'Ystart', 'Xend', 'Yend','seg_length', 'ID', 'rho', 'theta', 'Labels'] """
    #fig,ax=plt.subplots(1,3)
    clabel=np.unique(clusters['Labels'])
    
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
            if skipUnlinked:
                continue
            
        x,y=endpoints2(lines)

        x0=(max(x)-min(x))/2
        y0=(max(y)-min(y))/2
        size=np.sum(mask)
        rrange=max(lines['rho'])-min(lines['rho'])
        trange=(max(lines['theta'])-min(lines['theta']))%180
        avgrho=np.average(lines['rho'])
        avgtheta=np.average((lines['theta']))
        
        # check if it includes different signed angles
        if abs(np.sum(np.sign(lines['theta'].values))) < size: 
            crossZero=True
        else:
            crossZero=False
            
        
        stdrho=np.std(lines['rho'])
        stdt=np.std(lines['theta'])
        segmentL=lines['seg_length'].sum()
        
        
        # rotate data so all are at 20 deg then fit rectangle 
        # then calculate squares Error 
        #rotation_angle=-1*avgtheta #+20
        #rotatedLines=rotateData2(lines, rotation_angle)
        #print(avgtheta, rotation_angle)
        w,l,r,Xe, Ye, Xmid, Ymid=fit_Rec(lines, xc, yc)

        x1, x2, y1, y2=clustered_lines(x,y,avgtheta, l, xmid=Xmid, ymid=Ymid)

        
        hashlines=hash((lines['HashID'].values.tostring()))
        
        if l-w > 0: 
            EE=enEchelonAngleTwist(lines, avgtheta)
            tdiff=EE[0]
            EE_pvalue=EE[1]
        else: 
            tdiff=0
            EE_pvalue=1
        slope=-1/(np.tan(avgtheta)+0.00000000001)
        step=lines['seg_length'].min()+1
        

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
            
        if MaxDist/l < 0.5:
            Trust=True
        else:
            Trust=False
            
        clusters_data=clusters_data.append({ "Label": i, "Xstart": x1, "Ystart":y1, "Xend": x2,
                                            "Yend":y2, "X0": x0, "Y0": y0, "AvgRho":avgrho,
                                            "AvgTheta":avgtheta, "AvgSlope": slope, "AvgIntercept": b ,
                                            "RhoRange":rrange, "Aspect": l/w, 'Xmid': Xmid, 'Ymid': Ymid,
                                            "PerpOffsetDist": lines['PerpOffsetDist'].mean(),
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
        
    notclustered=len(clusters)-nclusters
    now = datetime.now()
    date = now.strftime("%d %b, %Y")
    clusters_data=clusters_data.assign(Date_Changed=date)
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
    
def ClusteredAll(dikeset,lines,cmask):
    notClustered=dikeset.iloc[~cmask]


def checkAllClusterChange(lines1, lines2):
    hash1=np.sort(lines1['ClusterHash'])
    hash2=np.sort(lines2['ClusterHash'])
    hash1.flags.writeable = False
    hash2.flags.writeable = False
    if hash(str(hash1)) == hash(str(hash2)):
        print("These clusters are the same")
        return True
    else: 
        print("These clusters are not the same")
        return False
    

def checkIndividualClusterChange(df1, df2):

    l1=df1['ClusterHash'].values
    l2=df2['ClusterHash'].values
    HashList=np.array([l1,l2])
    
    #pick the array with the longest length
    longest=max(HashList, key=lambda col: len(col))
    shortest=min(HashList, key=lambda col: len(col))

    same=np.in1d( longest,shortest) #outputs length of first input
    s=np.sum(same)
    eqLabels=np.arange(0,len(longest), 1)[same]
    diffLabels=np.arange(0,len(longest), 1)[~same]
  
    print(s, len(df1), len(df2))
    
    if (s==len(df1) and s==len(df2)):
        print("These clusters are ALL the same")
    else: 
        print("These clusters are not all the same")
        print("df1 has", str(len(df1)), 'while df2 has', str(len(df2)), "identified clusters")
        print("They have", str(s), "overlapping clusters")
        
    
    return eqLabels, diffLabels 
def checkClusterChange2(Segs, Line1, Line2):
    
    c1=np.unique(Line1['Label'])
    c2=np.unique(Line2['Label'])
    
    for i in c1: 
        lines=Segs[Segs['Label']==i]
        for j in c2: 
            
            x=lines['HashID'].values
            y=Segs[Segs['TrueLabel']==j]['HashID'].values
            result = np.where( x==y, x, 0)
            if np.sum(result) > 0:
                print(i,j, result, np.sum(result))
    

def RotateAndCluster(dikeset,dtheta, drho,**kwargs):
    dikesetOrig, ZOrig=HT_AGG_custom(dikeset, dtheta, drho)
    
    dikeset45=rotateData2(dikeset, 45)
    dikeset45=DikesetReProcess(dikeset45, HTredo=True)
    dikeset45, Z45=HT_AGG_custom(dikeset45, dtheta, drho)
    linesOrig=examineClusterShort(dikesetOrig)
    lines45=examineClusterShort(dikeset45)
    
    changedect=checkAllClusterChange(linesOrig, lines45)
    if changedect:
        return dikesetOrig
    if ~changedect:
        eqLabels, diffLabels=checkIndividualClusterChange(linesOrig, lines45)
        dikesetFinal=dikesetOrig.iloc[eqLabels[:len(dikesetOrig)]]

        return dikesetFinal

def errorAnalysis(lines, dikeset, plot=False):
    
    if plot: 
        fig,ax=plt.subplots(2,4)
        ax[0][0].set_ylabel('R_length')
        ax[0][0].scatter(lines['R_error'], lines['R_Length'])
        ax[0][1].scatter(lines['R_error'], lines['R_Width'])
        ax[0][1].set_ylabel('R_Width')
        
        ax[0][2].scatter(lines['R_error'], lines['Size'])
        ax[0][2].set_ylabel('Size')
        
        
        ax[0][3].hist(lines['R_error'], bins=50)
        ax[0][3].set_ylabel('Counts')
        
        ax[1][0].scatter(lines['R_error'], lines['AvgTheta'])
        ax[1][0].set_ylabel('AvgTheta')
        
        ax[1][1].scatter(lines['R_error'], lines['AvgRho'])
        ax[1][1].set_ylabel('AvgRho')
        
        ax[1][2].scatter(lines['R_error'], lines['ThetaRange'])
        ax[1][2].set_ylabel('ThetaRange')
        
        ax[1][3].scatter(lines['R_error'], lines['RhoRange'])
        ax[1][3].set_ylabel('RhoRange')
        
        for i in range(4):
            ax[1][i].set_xlabel("SS Error (m)")
    print("Error evaluation")
    print("average error:", lines['R_error'].mean())
    print("# clusters over 1 mil error:", max(lines['Label'])-len(lines))
    print("N% clustered", (np.sum(lines['Size'])/len(dikeset))*100)
        
def TopHTSection(lines, dikeset, rstep, tstep,n=1):
    fig,ax=plt.subplots(1,2)
    
    fig,ax[0],img=HThist(dikeset, rstep, tstep, ax=ax[0], fig=fig)
    h=img[0]
    
    xedges=img[1]
    yedges=img[2]
    
    
    [i,j]=np.unravel_index(h.argmax(), h.shape)
    
    
    xe=[xedges[i],xedges[i+1]]
    ye=[yedges[j],yedges[j+1]]
    print(h.max(), "clustered lines at", xe, "degrees", ye, "rho (m)" )
    #     masklat= (dikeset['latitude'] > lat1) & (dikeset['latitude'] < lat2)
    # masklong=(dikeset['longitude'] > lon1) & (dikeset['longitude'] < lon2)
    # masklatlong= (masklat==1) & (masklong==1)
    
    maskTheta=(lines['AvgTheta'] > xe[0]) & (lines['AvgTheta'] < xe[1])
    maskRho=(lines['AvgRho'] > ye[0]) & (lines['AvgRho'] < ye[1])
    mask=(maskTheta ==True) & (maskRho==True)
    ax[0].plot( [xe[0], xe[0], xe[1], xe[1], xe[0]], [ye[0], ye[1], ye[1], ye[0], ye[0]],'r')
    print(np.sum(lines['Size'].loc[mask]), "should equal", h.max())
    print(np.sum(lines['Size'].loc[mask])==h.max())
    toplines=lines[mask]
    plotlines(lines, 'k', ax[1], alpha=0.1)
    plotlines(toplines, 'r', ax[1])
    print(toplines['Size'].sum(), "dike segements")
    
    if n >1: 
        nextop=h.max()-1
        im,jm=np.where(h==nextop)
        for i,j in zip(im,jm):
            xe=[xedges[i],xedges[i+1]]
            ye=[yedges[j],yedges[j+1]]
           
            maskTheta=(lines['AvgTheta'] > xe[0]) & (lines['AvgTheta'] < xe[1])
            maskRho=(lines['AvgRho'] > ye[0]) & (lines['AvgRho'] < ye[1])
            mask=(maskTheta ==True) & (maskRho==True)
            toplines2=lines[mask]
            plotlines(toplines, 'r', ax[1])
            ax[0].plot( [xe[0], xe[0], xe[1], xe[1], xe[0]], [ye[0], ye[1], ye[1], ye[0], ye[0]],'r')
            toplines=toplines.append(toplines2)
    plt.tight_layout()
    
    return toplines, fig, ax




    


def extendLines(lines, save=False, name='Longlines.csv'):
    t,r=whichForm(lines)
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

def Run3Times(dikeset, dtheta, drho, plot=False, **kwargs):
    xc1,yc1=HT_center(dikeset)
    # run the first time 
    dikeset, Z1=HT_AGG_custom(dikeset, dtheta, drho, **kwargs)
    lines1, IC1=examineClusters(dikeset)
    
    # rotate 45 deg 
    dikeset2=rotateData2(dikeset, 45)
    theta2, rho2, xc2, yc2=AKH_HT(dikeset2)
    dikeset2['theta']=theta2
    dikeset2['rho']=rho2
    dikeset2, Z2=HT_AGG_custom(dikeset2, dtheta, drho, **kwargs)
    dikeset2=rotateData2(dikeset2,-45)
    lines2, IC2=examineClusters(dikeset2)
    
    #move to lower left
    dx=np.min([dikeset['Xstart'].min(), dikeset['Xend'].min()])
    dy=np.min([dikeset['Ystart'].min(), dikeset['Yend'].min()])
    theta3, rho3, xc3, yc3=AKH_HT(dikeset, xc=dx,yc=dy)
    dikeset3=dikeset.copy()
    dikeset3['theta']=theta3
    dikeset3['rho']=rho3
    
    dikeset3, Z3=HT_AGG_custom(dikeset3, dtheta, drho, **kwargs)
    
    lines3, IC3=examineClusters(dikeset3)
    
    #check the changes 
    eq12, diff12=checkClusterChange(lines1, lines2)
    eq13, diff13=checkClusterChange(lines1, lines3)
    eq32, diff32=checkClusterChange(lines3, lines2)
    
    Flines1=FilterLines(lines1)
    Flines2=FilterLines(lines2)
    Flines3=FilterLines(lines3)
    
    print("Comparing the filtered lines")
    eq12, diff12=checkClusterChange(Flines1, Flines2)
    eq13, diff13=checkClusterChange(Flines1, Flines3)
    eq32, diff32=checkClusterChange(Flines3, Flines2)
    
    if plot: 
        
        fig,ax=plt.subplots(1,2)
        plotlines(lines1, 'k', ax[0])
        plotlines(lines2, 'b', ax[0])
        plotlines(lines3, 'r', ax[0])
        
        ax[1].scatter(dikeset['theta'], dikeset['rho'], c=dikeset['Labels'].values, alpha=0.6, cmap=cm.Greys, marker='p')
        ax[1].scatter(dikeset2['theta'], dikeset2['rho'], c=dikeset2['Labels'].values, alpha=0.6, cmap=cm.Blues, marker='*')
        ax[1].scatter(dikeset3['theta'], dikeset3['rho'], c=dikeset3['Labels'].values, alpha=0.6, cmap=cm.Reds, marker='^')
        
        fig,ax=plt.subplots(3,2)
        
        plotlines(lines1.iloc[diff12], 'k', ax[0,0])
        plotlines(lines2.iloc[diff12], 'b', ax[0,1])
        
        plotlines(lines1.iloc[diff13], 'k', ax[1,0])
        plotlines(lines3.iloc[diff13], 'r', ax[1,1])
        
        plotlines(lines2.iloc[diff32[diff32<len(lines2)]], 'b', ax[2,0])
        plotlines(lines3.iloc[diff32], 'r', ax[2,1])
        
        
        
    return lines1, lines2, lines3
    

def persitance(df):
    t,r=whichForm(df)

    theta=df[t].values
    rho=df[r].values

    X2D = (np.vstack( (theta, rho-np.mean(rho))) ).T
    
    #use the scaled version of the distance metric 
    dtheta=2 
    drho=df['seg_length'].mean()
    metric= lambda x,y: CyclicEuclideanScaled(x,y,dtheta,drho)
    
    dist = squareform(pdist(X2D, metric))
    condensedD = squareform(dist)
    
    Y = sch.linkage(condensedD, method='complete')
    Z1 = sch.dendrogram(Y, orientation='left', no_plot=True)
    
    dcoord=np.array(Z1['dcoord'])
    icoord=np.array(Z1['icoord'])
    c=Z1['color_list']
    idx=Z1['leaves']
    
    
    #scaling for persistance
    #a1=(np.max(icoord)+np.max(dcoord))/2
    #a0=(np.min(icoord)+np.min(dcoord))/2
    
    #dcoord=(dcoord-a0)/(a1-a0)
    #icoord=(icoord-a0)/(a1-a0)


    x=np.max(dcoord)
    fig,ax=plt.subplots(2)
    ax[0].plot([1,1], [x,x], 'k-', linewidth=10)
    p=np.append(dcoord[:,1]-dcoord[:,0], dcoord[:,2]-dcoord[:,3])
    birth=np.array([ dcoord[:,0], dcoord[:,3]])+1
    death=np.array([ dcoord[:,1], dcoord[:,2]])+1
    
    ax[0].plot(birth, death, "*")
    
    ax[0].set_yscale('log')
    ax[0].set_xscale('log')
    ax[0].plot([1,x], [1,x], 'k-', linewidth=4)
    
    ax[1].hist(np.log(p+1), bins=20, color='r')
    
    ax[1].set_ylabel('Counts')
    ax[1].set_xlabel('Persistance')
    ax[0].set_xlabel('Birth')
    ax[0].set_ylabel('Death')
    # for ys, color in zip(dcoord, c):
    #     #ax[0].plot(xs, ys, color)
        
    #     birth=np.array([ys[0]+1, ys[3]+1])
    #     death=np.array([ys[1]+1, ys[2]+1])
    #     ax[0].plot(birth,death, "*", color=color)
    #     p=np.append(birth-death)
    ax[0].set_yscale('log')
    return fig, ax, Z1

def testValidity(lines, dtheta, drho):
    
    if any(lines['ThetaRange'] > dtheta):
        print("Failed Theta Validity Check 1")
    else: 
        print("Passed Theta Validity Check 1")
        
    if any(lines['RhoRange'] > drho):
        print("Failed Theta Validity Check 1")
    else: 
        print("Passed Theta Validity Check 1")
        
def OutputRectangles(clusters):

    clabel=np.unique(clusters['Labels'])
    nclusters=len(clabel)
    notclustered=sum([clusters['Labels']==-1][0])
    xc,yc=HT_center(clusters)
    clusters_data=pd.DataFrame()
    Xs=np.zeros((nclusters,5))
    Ys=np.zeros((nclusters,5))
    
    for i in np.unique(clusters['Labels']): 
        clustered=True
        mask=clusters['Labels']==i
        lines=clusters[mask]
        
        if (i == -1 or len(lines)<2):
            clustered=False
            #continue
            
        w,l,r,Xe, Ye, Xmid, Ymid=fit_Rec(lines, xc, yc)
        Xs[i-1,:]=Xe
        Ys[i-1,:]=Ye
        
    return Xs, Ys

def SensitivityAnalysis(dikeset, path=None):
    
    segL=np.median(dikeset['seg_length'].values)
    isExist = os.path.exists(path)
    IC_total=pd.DataFrame
    for dtheta in [1,2,5]:
        for drho in [0.25*segL, segL, 2*segL]:
            df_filename=path+'dtheta_'+str(dtheta)+'drho'+str(drho)+'.csv'
            lines_filename=path+'LINES_dtheta_'+str(dtheta)+'drho'+str(drho)+'.csv'
            Z_filename=path+'LinkageMat'+'dtheta_'+str(dtheta)+'drho'+str(drho)
            IC_filename=path+'Analysis'+'dtheta_'+str(dtheta)+'drho'+str(drho)
            if not isExist: 
                    dikeset2, Z=HT_AGG_custom(dikeset, dtheta, drho, linkage='complete')
                    lines,IC=examineClusters(dikeset)
                    dikeset2.to_csv(df_filename, index=False)
                    np.save(Z_filename,Z)
                    lines.to_csv(lines_filename, index=False)
                    IC.to_csv(IC_filename, index=False)
            else:
                #dikeset2=pd.read_csv(df_filename)
                #lines=pd.read_csv(lines_filename)
                IC=pd.read_csv(IC_filename)
                IC_total=IC_total.append(IC)
                
                    
                    