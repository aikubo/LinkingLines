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
from plotmod import HThist, plotlines, DotsHT
from sklearn.neighbors import NearestNeighbors
from sklearn.preprocessing import scale
from PrePostProcess import *
from htMOD import rotateData2
from clusterMod import CyclicEuclideanScaled, CyclicAngleDist

from PrePostProcess import completePreProcess, whichForm, midPoint
import scipy.cluster.hierarchy as sch
from scipy.spatial.distance import squareform
from scipy.cluster.hierarchy import dendrogram
from scipy.spatial.distance import pdist, squareform

import matplotlib.pyplot as plt 
from scipy import stats

def clustered_lines(xs, ys, theta):
    xstart=max(xs)
    ystart=max(ys)
    
    xend=min(xs)
    yend=min(ys)
    
    xmid=(xstart+xend)/2
    ymid=(ystart+yend)/2
    l=np.sqrt( (xstart-xend)**2 + (ystart-yend)**2)
    a = np.cos(np.deg2rad(theta))
    b = np.sin(np.deg2rad(theta))
    
    x0 = xmid
    y0 = ymid
    x1 = int(x0 + l/2 * (-b))
    y1 = int(y0 + l/2 * (a))
    x2 = int(x0 - l/2 * (-b))
    y2 = int(y0 - l/2 * (a))
    
    return x1, x2, y1, y2

def checkoutCluster(dikeset, label):
    fig, ax=plt.subplots(1,2) 
    mask=(dikeset['Labels']==label)
    lines=dikeset[mask]
    ax[0].scatter(dikeset['theta'], dikeset['rho'], c='k', alpha=0.1)
    ax[0].scatter(lines['theta'], lines['rho'], c='r')
    plotlines(lines, 'r', ax[1])
    
    
    return fig,ax

def enEchelon(lines, avgtheta):
    Xmid=lines['Xmid'].to_numpy()
    Ymid=lines['Ymid'].to_numpy()
    
    slope, intercept, r_value, p_value, std_err = stats.linregress(Xmid, Ymid)
    thetaM=np.rad2deg(np.arctan(-1/(slope+ 0.0000000000000001)))
    Xstart=min(Xmid)
    Xend=max(Xmid)
    Ystart=slope*Xstart+intercept 
    Yend=slope*Xend+intercept 
    tdiff=CyclicAngleDist([avgtheta], [thetaM])
    #print(tdiff, thetaM, p_value)
    if p_value > 0.05:
        tdff=0
    
    
    return tdiff, Xstart, Xend, Ystart, Yend


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
        size=sum(mask)
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
        rotation_angle=-1*avgtheta+20
        rotatedLines=rotateData2(lines, rotation_angle)
        #print(avgtheta, rotation_angle)
        w,l=fit_Rec(rotatedLines, xc, yc)
        r=squaresError(rotatedLines,xc,yc)
        #f r> 1000000: 
        #   continue 
        x1, x2, y1, y2=clustered_lines(x,y,avgtheta)
        #points=(np.vstack((x, y)).T)
        #sx,sy, ang=ellipse(points)
        Xe,Ye=RecEdges(lines, xc, yc)
        hashlines=hash((lines['HashID'].values.tostring()))
        tdiff, EXstart, EXend, EYstart, EYend=enEchelon(lines, avgtheta)
        clusters_data=clusters_data.append({ "Label": i, "Xstart": x1, "Ystart":y1, "Xend": x2,
                                            "Yend":y2, "X0": x0, "Y0": y0, "AvgRho":avgrho,
                                            "AvgTheta":avgtheta, "RhoRange":rrange, 
                                            "PerpOffsetDist": lines['PerpOffsetDist'].mean(),
                                            "ThetaRange": trange, "StdRho": stdrho, 
                                            "StdTheta": stdt, "R_Width": w, "R_Length": l, 
                                            "Size":size, "R_error":np.sqrt(r), "Linked":clustered,
                                            "SegmentLSum":segmentL, "Hash":hashlines, 
                                            'ClusterCrossesZero': crossZero,
                                            "EnEchelonAngleDist":tdiff}, ignore_index=True)
    
    
        if tdiff > enEchelonCutofff: 
            EEDikes=EEDikes.append({ "Label": i, "Xstart": EXstart, "Ystart":EYstart, "Xend": EXend,
                                            "Yend":EYend, "CXstart": x1, "CYstart":y1, "CXend": x2,
                                            "CYend":y2, "X0": x0, "Y0": y0, "AvgRho":avgrho,
                                            "AvgTheta":avgtheta, "RhoRange":rrange, 
                                            "PerpOffsetDist": lines['PerpOffsetDist'].mean(),
                                            "ThetaRange": trange, "StdRho": stdrho, 
                                            "StdTheta": stdt, "R_Width": w, "R_Length": l, 
                                            "Size":size, "R_error":np.sqrt(r), "Linked":clustered,
                                            "SegmentLSum":segmentL, "Hash":hashlines, 
                                            'ClusterCrossesZero': crossZero,
                                            "EnEchelonAngleDist":tdiff}, ignore_index=True)
    if 'Formation' in clusters.columns: 
        print("adding formation")
        
        for j in np.unique(clusters['Formation']):
            clusters_data[j]=''
        for i in np.unique(clusters['Labels']): 
            clustered=True
            mask=clusters['Labels']==i
            lines=clusters.loc[mask]
            forms, counts=np.unique(lines['Formation'], return_counts=True)
            print(forms, counts)
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
                              "AverageThetaRange": np.average(clusters_data["ThetaRange"]),
                              "StdRhoRange": np.std(clusters_data["RhoRange"]),
                              "StdThetaRange": np.std(clusters_data["ThetaRange"]),
                              "AvgClusterSize": np.average(clusters_data["Size"]),
                              "ClusterSizeStd": np.std(clusters_data["Size"]),
                              "ClusterMax": clusters_data["Size"].max(), 
                              "AverageL": np.average(clusters_data["R_Length"]),
                              "AverageW": np.average(clusters_data["R_Width"])},
                              "NEEDikes":len(EEDikes)
                              index=[0])
    
    #print(evaluation)
    #plt.tight_layout()
    if ifEE:
        return clusters_data, evalution, EEDikes
    else:
        return clusters_data, evaluation

def ClusteredAll(dikeset,lines,cmask):
    notClustered=dikeset.iloc[~cmask]

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
        
def TopHTSection(lines, rstep, tstep):
    fig,ax=plt.subplots(1,2)
    
    ax[0],img=HThist(lines['AvgRho'], lines['AvgTheta'], rstep=rstep, tstep=tstep, ax=ax[0])
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

    print(np.sum(mask), "should equal", h.max())
    print(np.sum(mask)==h.max())
    toplines=lines[mask]
    plotlines(lines, 'k', ax[1], alpha=0.1)
    plotlines(toplines, 'r', ax[1])
    print(toplines['Size'].sum(), "dike segements")
    
    return toplines

def plotlabel(df, linked, label, hashlabel=False):
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
    plotlines(linkedlines,"g", ax)
    
    
    ax.set_title("Label "+str(label)+ "| theta:"+str(linkedlines['AvgTheta'].values)+" rho:"+str(linkedlines['AvgRho'].values))
    ax.text( .89, .89, "W:"+str(linkedlines['R_Width'].values),transform=ax.transAxes )
    ax.text( .84, .84, "L:"+str(linkedlines['R_Length'].values),transform=ax.transAxes  )
    ax.text( .80, .80, "errors:"+str(linkedlines['R_error'].values),transform=ax.transAxes  )
    print("errors:"+str(linkedlines['R_error']))
    
def checkClusterChange(df1, df2):
    s=0
    eqLabels=[]
    oneTotwo=[]
    errChange=[]
    ndikes=0
    for i in range(len(df1)):
        for j in range(len(df2)):
            if df1['Hash'].iloc[i]==df2['Hash'].iloc[j]:
                s=s+1
                ndikes=df1['Size'].iloc[i]+ndikes
                eqLabels.append(df1['Hash'].iloc[i])
                #errChange.append(df1['R_error'].iloc[i]-df2['R_error'].iloc[j])
                print(ndikes)
                #eqLabels=np.append(eqLabels, [df1['Label'].iloc[i], df2['Label'].iloc[j]])
    print(s, len(df1), len(df2), ndikes)
    
    if (s==len(df1) and s==len(df2)):
        print("These clusters are the same")
    else: 
        print("These clusters are not the same")
        print("df1 has", str(len(df1)), 'while df2 has', str(len(df2)), "identified clusters")
        print("They have", str(s), "overlapping clusters")
        
    
    return eqLabels #, errChange
    
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