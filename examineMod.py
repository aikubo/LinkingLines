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
import statsmodels.api as sm
import matplotlib.ticker as mticker
import matplotlib.pyplot as plt 
from scipy import stats
import warnings
from matplotlib import cm
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

def checkoutCluster(dikeset, label):
    
    xc=dikeset['xc'].iloc[0]
    yc=dikeset['yc'].iloc[0]
    
    fig, ax=plt.subplots(1,2) 
    fig.set_size_inches( 10,5)
    fig.suptitle("Label:"+str(label))
    mask=(dikeset['Labels']==label)
    lines=dikeset[mask]
    
    xlim1=min(lines['theta'])-2
    xlim2=max(lines['theta'])+2
    ylim1=min(lines['rho'])-10000
    ylim2=max(lines['rho'])+10000
    
    ax[0].scatter(dikeset['theta'], dikeset['rho'], c='grey', alpha=0.1) 
    ax[0].scatter(lines['theta'], lines['rho'], c='r')
    ax[0].plot(lines['theta'].mean(), lines['rho'].mean(), 'g*')
    
    ax[0].plot([xlim1, xlim1, xlim2, xlim2, xlim1], [ylim1,ylim2, ylim2, ylim1, ylim1], 'k')
    left, bottom, width, height = [0.2, 0.2, 0.2, 0.2]
    ax2 = fig.add_axes([left, bottom, width, height])
    ax2.scatter(lines['theta'], lines['rho'], c='r')
    ax2.plot(lines['theta'].mean(), lines['rho'].mean(), 'g*')
    ax2.scatter(dikeset['theta'], dikeset['rho'], c='grey', alpha=0.1) 
    ax2.set_ylim([ylim1+8000, ylim2-8000])
    ax2.set_xlim([xlim1+1.5, xlim2-1.5])
    
    ax[0].set_xlabel('Theta ($^\circ$)')
    ax[0].set_ylabel('Rho (m)')
    
    fig,ax[1], l,w=pltRec(lines, xc,yc, fig=fig, ax=ax[1])
    
  
    left, bottom, width, height = [0.2, 0.7, 0.2, 0.2]
    
    if lines['theta'].mean()*lines['rho'].mean() > 0:
        loc=4
        tlocx=.05
        
    else: 
        
        loc=3
        tlocx=.60
    
    ax4 = inset_axes(ax[1], width="30%", height="30%", loc=loc)

    plotlines(dikeset, 'grey', ax4, alpha=0.1)
    plotlines(lines, 'r', ax4)
    ax4.xaxis.tick_top()
    ax4.xaxis.set_label_position('top') 
    
    if lines['theta'].mean()*lines['rho'].mean() < 0:
        ax4.yaxis.tick_right()
        ax4.yaxis.set_label_position('right') 
        

    FixCartesianLabels(ax4)
    FixCartesianLabels(ax[1])
    
    tlocy=.95
    ax[0].text(tlocx,tlocy-0.10,'Angle Mean: '+str( round(lines['theta'].mean(), 2)), transform=ax[0].transAxes)
    ax[0].text(tlocx,tlocy-0.05,'Rho Mean: '+str(round(lines['rho'].mean(),1)), transform=ax[0].transAxes)
    ax[0].text(tlocx,tlocy,'Size:'+str(len(lines['rho'])), transform=ax[1].transAxes)
    
    ax[0].text(tlocx,tlocy-0.20,'Angle Range: '+str( round(np.ptp(lines['theta']), 2)), transform=ax[0].transAxes)
    ax[0].text(tlocx,tlocy-0.15,'Rho Range: '+str(round(np.ptp(lines['rho']),1)), transform=ax[0].transAxes)
    
    ax[1].text(tlocx,tlocy,'Length: '+str(round(l,1)), transform=ax[1].transAxes)
    ax[1].text(tlocx,tlocy-0.05,'Width: '+str(round(w,1)), transform=ax[1].transAxes)
    ax[1].text(tlocx,tlocy-0.10,'Aspect: '+str(round(l/w,1)), transform=ax[1].transAxes)
    ax[1].text(tlocx,tlocy-0.15,'Size: '+str(len(lines)), transform=ax[1].transAxes)

    return fig,ax

# def overlap(min1, max1, min2, max2):
#     return max(0, min(max1, max2) - max(min1, min2))

def overlapSegments(lines):
    """
    Calculate overlap between lines 
    

    """
    Xstart=lines['Xstart'].to_numpy()
    Ystart=lines['Xstart'].to_numpy()
    Xend=lines['Xend'].to_numpy()
    Yend=lines['Yend'].to_numpy()
    step=lines['seg_length'].min()+1
    totalL=np.sum( np.sqrt( (Xstart-Xend)**2 + (Ystart-Yend)**2))
    xs=np.floor(np.concatenate([np.arange(min(x,y), max(x,y), step) for x,y in zip(Xstart,Xend)]))
    ys=np.floor(np.concatenate([np.arange(min(x,y), max(x,y), step) for x,y in zip(Ystart,Yend)]))
    u,ycounts=np.unique(ys, return_counts=True)
    u,xcounts=np.unique(xs, return_counts=True)
    
    overlapx=np.sum(xcounts[xcounts>1])*step
    overlapy=np.sum(ycounts[ycounts>1])*step
    
    overlap=np.sqrt(overlapx**2 +overlapy**2)
    #for i in len(Xstart):
        
    
    return overlap/totalL

def overlapSegments2(Xstart, Ystart, Xend, Yend, slope,step):
    """
    Calculate overlap between lines 
    

    """
    
    xs=[np.arange(min(x,y), max(x,y), step/2) for x,y in zip(Xstart,Xend)]
    #totalL=int(np.sum( np.sqrt( (Xstart-Xend)**2 + (Ystart-Yend)**2)))
    
    arr=np.zeros( (len(xs),)+xs[0].shape)
    for i,v in enumerate(xs):
        arr[i]=np.floor(v)
    u,xcounts=np.unique(arr, return_counts=True)
    
    overlapx=np.sum(xcounts[xcounts>1])*step
    overlapy=slope*overlapx
    
    overlap=np.sqrt(overlapx**2 +overlapy**2)
    #for i in len(Xstart):
        
    
    return overlap

def overlapSegments4(Xstart, Xend, slope,step, b):
    """
    Calculate overlap between lines 
    

    """
    step=1
    xs=[np.arange(min(x,y), max(x,y), step) for x,y in zip(np.floor(Xstart),np.ceil(Xend))]
    
    l=np.max([len(xi) for xi in xs])
    xs_sameLength=[np.append(xi,[np.nan]*(l-len(xi))) for xi in xs]

    arr=np.vstack(xs_sameLength) #better
    u,xcounts=np.unique(arr[~np.isnan(arr)], return_counts=True)
    
    overlapx=np.sum(xcounts>1)*step
    overlapy=slope*overlapx
    
    overlap=np.sqrt(overlapx**2 +overlapy**2)
        
    
    return overlap




def enEchelon(d, avgtheta):
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



def examineClusters(clusters, enEchelonCutofff=7, ifEE=False):
    """ Must have  ['Xstart', 'Ystart', 'Xend', 'Yend','seg_length', 'ID', 'rho', 'theta', 'Labels'] """
    #fig,ax=plt.subplots(1,3)
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
        
        if l/w > 1: 
            EE=enEchelon(lines, avgtheta)
            tdiff=EE[0]
            EE_pvalue=EE[1]
        else: 
            tdiff=0
            EE_pvalue=1
        slope=-1/(np.tan(avgtheta)+0.00000000001)
        step=lines['seg_length'].min()+1
        

        b=avgrho/np.sin((np.deg2rad(avgtheta))+0.00000000001)+yc-xc*slope
        
        if size>1: 
            
            if segmentL>l:
                overlap=segmentL-l
            else:
                overlap=overlapSegments4(x[0:size], x[size:size*2], slope,step, b)
        else:
            overlap=0
            
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
            
        clusters_data=clusters_data.append({ "Label": i, "Xstart": x1, "Ystart":y1, "Xend": x2,
                                            "Yend":y2, "X0": x0, "Y0": y0, "AvgRho":avgrho,
                                            "AvgTheta":avgtheta, "AvgSlope": slope, "AvgIntercept": b ,
                                            "RhoRange":rrange, 
                                            "PerpOffsetDist": lines['PerpOffsetDist'].mean(),
                                            "ThetaRange": trange, "StdRho": stdrho, 
                                            "StdTheta": stdt, "R_Width": w, "R_Length": l, 
                                            "Size":size, "R_error":np.sqrt(r), "Linked":clustered,
                                            "SegmentLSum":segmentL, "ClusterHash":hashlines, 
                                            'ClusterCrossesZero': crossZero,
                                            "EnEchelonAngleDiff":tdiff, "Overlap": overlap/segmentL,
                                            "EEPvalue":EE_pvalue, "MaxSegNNDist": MaxDist/l,
                                            "MedianSegNNDist": MedianDist/l, "MinSegNNDist": MinDist/l}, ignore_index=True)
    
    
        # if tdiff > enEchelonCutofff: 
        #     EEDikes=EEDikes.append({ "Label": i, "Xstart": EXstart, "Ystart":EYstart, "Xend": EXend,
        #                                     "Yend":EYend, "CXstart": x1, "CYstart":y1, "CXend": x2,
        #                                     "CYend":y2, "X0": x0, "Y0": y0, "AvgRho":avgrho,
        #                                     "AvgTheta":avgtheta, "RhoRange":rrange, 
        #                                     "PerpOffsetDist": lines['PerpOffsetDist'].mean(),
        #                                     "ThetaRange": trange, "StdRho": stdrho, 
        #                                     "StdTheta": stdt, "R_Width": w, "R_Length": l, 
        #                                     "Size":size, "R_error":np.sqrt(r), "Linked":clustered,
        #                                     "SegmentLSum":segmentL, "Hash":hashlines, 
        #                                     'ClusterCrossesZero': crossZero,
        #                                     "EnEchelonAngleDiff":tdiff, "Overlap":overlap/segmentL}, ignore_index=True)
    if 'Formation' in clusters.columns: 
        print("adding formation")
        
        for j in np.unique(clusters['Formation']):
            clusters_data[j]=''
        for i in np.unique(clusters['Labels']): 
            clustered=True
            mask=clusters['Labels']==i
            lines=clusters.loc[mask]
            forms, counts=np.unique(lines['Formation'], return_counts=True)
           
            for p in forms: 
                t=np.where(p==forms)[0][0]
                clusters_data[p][i]=counts[t]
            
    #Find KNN distance of the endpoints 
    # Check that the KNN is not just one endpoint
    # Add to lines output
    # X,Y=endpoints2(clusters_data)
    # #X=clusters_data['Xstart'].values
    # #Y=clusters_data['Ystart'].values
    # X=(np.vstack((X, Y)).T)
    # nbrs = NearestNeighbors(n_neighbors=3, algorithm='ball_tree').fit(X)
    # distances, indices = nbrs.kneighbors(X)
    # s=0
    # distAvg=[]
    # for i in range(int(len(indices)/2)): 
    #     if distances[i,1]==clusters_data['R_Length'].iloc[i]:
    #         s=s+1
    #         print("matched to endpoint", s)
    #         distAvg.append(min(distances[i,2],distances[i*2,2]))
    #         continue
            
    #     distAvg.append(min(distances[i,1],distances[i*2,1]))

    #clusters_data["KNN2"]=distAvg
    
    # X,Y=midpoint(clusters_data)
    # #X=clusters_data['Xstart'].values
    # #Y=clusters_data['Ystart'].values
    # X=(np.vstack((X, Y, clusters_data['AvgTheta'])).T)
    # X=np.array(clusters_data['AvgTheta']).reshape(-1,1)
    # #=scale(X)
    # #rint(X.shape)
    # nbrs = NearestNeighbors(n_neighbors=3, algorithm='ball_tree').fit(X)
    # distances, indices = nbrs.kneighbors(X)
    # #s=0
    # #distAvg=[]
    # for i in range(int(len(indices)/2)): 
    #     if distances[i,1]==clusters_data['R_Length'].iloc[i]:
    #         s=s+1
    #         print("matched to endpoint", s)
    #         distAvg.append(min(distances[i,2],distances[i*2,2]))
    #         continue
            
    #     distAvg.append(min(distances[i,1],distances[i*2,1]))
    # print(len(clusters_data), len(distances[:,1]))

    # clusters_data["KNN2"]=distances[:,1]
    

    # X=clusters_data['Xstart'].values
    # Y=clusters_data['Ystart'].values
    # X=(np.vstack((X, Y)).T)
    # nbrs = NearestNeighbors(n_neighbors=3, algorithm='ball_tree').fit(X)
    # distances, indices = nbrs.kneighbors(X)
    # s=0
    # distAvg=[]
    # for i in range(int(len(indices))): 
    #     if distances[i,1]==clusters_data['R_Length'].iloc[i]:
    #         s=s+1
    #         print("matched to endpoint", s)
    #         distAvg.append(distances[i,2])
    #         continue
            
    #     distAvg.append(distances[i,1])
        
    # clusters_data["KNN"]=distAvg

    clusters_data.astype({'Linked':'bool'})
    
    evaluation=pd.DataFrame({ "Ic": (len(clusters)-notclustered),
                              "nClusters": nclusters,
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
                              "NEEDikes":len(EEDikes),
                              "MaxEEAngleDiff":np.max(clusters_data["EnEchelonAngleDiff"]),
                              "AverageEAngleDiff":np.average(clusters_data["EnEchelonAngleDiff"])},
                              index=[0])
    
    #print(evaluation)
    #plt.tight_layout()
    if ifEE:
        return clusters_data, evaluation, EEDikes
    else:
        return clusters_data, evaluation

def ClusteredAll(dikeset,lines,cmask):
    notClustered=dikeset.iloc[~cmask]

def FilterLines(lines, size=4, MaxSegNNDist=0.5):
    mask= (lines['Size']>3)*(lines['MaxSegNNDist']<0.45)
    return lines[mask]


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
    print(np.sum(mask), "should equal", h.max())
    print(np.sum(mask)==h.max())
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
    
    return toplines, fig, ax

def plotlabel(df, linked, label, hashlabel=False, EE=False):
    fig,ax=plt.subplots()
    print(label)
    if hashlabel:
        l=linked['Label'].loc[linked['Hash']==label].astype(int).values[0]
        label=l 
        print(label)

    mask=df['Labels']==label
    xc,yc=HT_center(df)
    lines=df[mask]
    plotlines(lines,"r", ax)
    
    # x,y=endpoints2(lines)
    # w,l=fit_Rec(lines, xc, yc)
    # r=squaresError(lines,xc,yc)
    # avgtheta=np.average(lines['theta'])

    # x1, x2, y1, y2=clustered_lines(x,y,avgtheta)
    linkedlines=linked[linked['Label']==label]
    plotlines(linkedlines,"k", ax)
    
    
    ax.set_title("Label "+str(label)+ "| theta:"+str(linkedlines['AvgTheta'].values)+" rho:"+str(linkedlines['AvgRho'].values))
    ax.text( .80, .89, "W:"+str(linkedlines['R_Width'].values),transform=ax.transAxes )
    ax.text( .80, .84, "L:"+str(linkedlines['R_Length'].values),transform=ax.transAxes  )
    ax.text( .80, .80, "errors:"+str(linkedlines['R_error'].values),transform=ax.transAxes  )
    print("errors:"+str(linkedlines['R_error']))
    
    if EE: 
        Xmid=lines['Xmid'].to_numpy()
        Ymid=lines['Ymid'].to_numpy()
        tdiff, EXstart, EXend, EYstart, EYend=enEchelon(lines, linkedlines['AvgTheta'].values[0])
        ax.plot( [EXstart, EXend], [EYstart, EYend], 'g*-')
        ax.plot(Xmid, Ymid, 'bp', markersize=2)
        
        
    ax.axis('equal')

def checkLineClusterChange(lines1, lines2):
    hash1=np.sort(lines1['Hash'])
    hash2=np.sort(lines2['Hash'])
    hash1.flags.writeable = False
    hash2.flags.writeable = False
    if hash(str(hash1)) == hash(str(hash2)):
        print("These clusters are the same")
    else: 
        print("These clusters are not the same")
    
    return 

def checkClusterChange(df1, df2):

    
    l1=df1['ClusterHash'].values
    l2=df2['ClusterHash'].values
    same=np.in1d(l1,l2)
    s=np.sum(same)
    eqLabels=np.arange(0,len(df1), 1)[same]
    diffLabels=np.arange(0,len(df1), 1)[~same]
    # for i in range(len(df1)):
    #     for j in range(len(df2)):
    #         if df1['HashID'].iloc[i]==df2['HashID'].iloc[j]:
    #             s=s+1
    #             #ndikes=df1['Size'].iloc[i]+ndikes
    #             eqLabels.append(df1['HashID'].iloc[i])
    #             #errChange.append(df1['R_error'].iloc[i]-df2['R_error'].iloc[j])
    #             #print(ndikes)
    #             #eqLabels=np.append(eqLabels, [df1['Label'].iloc[i], df2['Label'].iloc[j]])
    print(s, len(df1), len(df2))
    
    if (s==len(df1) and s==len(df2)):
        print("These clusters are ALL the same")
    else: 
        print("These clusters are not all the same")
        print("df1 has", str(len(df1)), 'while df2 has', str(len(df2)), "identified clusters")
        print("They have", str(s), "overlapping clusters")
        
    
    return eqLabels, diffLabels #, errChange
    
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
