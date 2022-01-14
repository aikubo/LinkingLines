#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  1 13:12:50 2021

@author: akh
"""

from sklearn.preprocessing import scale, RobustScaler, PowerTransformer
from sklearn.cluster import AgglomerativeClustering as AGG
from sklearn.cluster import DBSCAN
from sklearn.preprocessing import StandardScaler
from sklearn.neighbors import NearestNeighbors
import numpy as np 
import pandas as pd
from scipy.cluster.hierarchy import dendrogram
from htMOD import *
from examineMod import *
from sklearn.decomposition import PCA
from sklearn.metrics.pairwise import pairwise_distances
from PrePostProcess import whichForm

def CylDistanceAngleModified(u, v):
    '''
    Calculates cylindrical distance between two points in 3D with angle modification

    Input:
        u: [theta, rho, z]
        v: [theta, rho, z]
    Output:
        dist: distance between u and v
    '''
    
    dist1 = u[0]-v[0]
    dist2 = dist1%180
    deltaT = min(dist1, dist2)
    return np.sqrt(u[1]**2 + v[1]**2 + 2*u[1]*v[1]*np.cos(np.deg2rad(deltaT)) + (u[2]-v[2])**2)

def PolarDistanceAngleModified(u, v):
    '''
    Calculates polar distance between two points in 2D with angle modification

    Input:  
        u: [theta, rho]
        v: [theta, rho]
    Output:
        dist: distance between u and v

    '''
    
    dist1 = abs(u[0]-v[0])
    dist2 = dist1%180
    deltaT = min(dist1, dist2)
    return np.sqrt(u[1]**2 + v[1]**2 + 2*u[1]*v[1]*np.cos(np.deg2rad(deltaT)))

def DistortionDistance(u,v):
    '''
    Calculates distance based on distortion calculation in 2D

    Input:
        u: [theta, rho]
        v: [theta, rho]
    Output:
        dist: distance between u and v

        For u==v, dist=0
        for u[0]==v[0], dist=u[1]-v[1]
        for u[1]==v[1] and dtheta-> 0, dist ->0
        )
    '''
    
    
    dtheta=abs(u[0]-v[0])%180
    d1=v[1]-u[1]/np.cos(np.deg2rad(dtheta))
    d2=v[1]/np.cos(np.deg2rad(dtheta))-u[1]
    
    return abs((1+dtheta)*max(d1,d2))

def DistortionDistance2(u,v):
    '''
    Calculates distance based on distortion calculation in 2D

    Input:
        u: [theta, rho]
        v: [theta, rho]
    Output:
        dist: distance between u and v

        For u==v, dist=0
        for u[0]==v[0], dist=u[1]-v[1]
        for u[1]==v[1] and dtheta-> 0, dist ->0
        )
    '''
    dtheta=abs(u[0]-v[0])%180
    d1=v[1]-u[1]/np.cos(np.deg2rad(dtheta))
    d2=v[1]/np.cos(np.deg2rad(dtheta))-u[1]
    
    return abs((1+dtheta)*(max(d1,d2)/min(d1,d2)))


def plotDendro(dist1, labels, title):

    #https://stackoverflow.com/questions/2982929/plotting-results-of-hierarchical-clustering-ontop-of-a-matrix-of-data-in-python    
    D=dist1
    condensedD = squareform(D)
    
    # Compute and plot first dendrogram.
    fig = pylab.figure(figsize=(8,8))
    
    ax1 = fig.add_axes([0.09,0.1,0.15,0.6])
    ax1.set_title(title)
    Y = sch.linkage(condensedD, method='complete')
    Z1 = sch.dendrogram(Y, labels=labels, orientation='left')
    #ax1.set_xticks()
    #ax1.set_yticks([])
    
    # Compute and plot second dendrogram.
    ax2 = fig.add_axes([0.3,0.76,0.6,0.2])
    Y = sch.linkage(condensedD, method='complete')
    Z2 = sch.dendrogram(Y,  labels=labels)
    #ax2.set_xticks(labels[Z1['leaves']])
    #ax2.set_yticks([])
    
    # Plot distance matrix.
    #add axis(left, bottom, width, height)
    axmatrix = fig.add_axes([0.3,0.1,0.6,0.6])
    idx1 = Z1['leaves']
    idx2 = Z2['leaves']
    D = D[idx1,:]
    D = D[:,idx2]
    im = axmatrix.matshow(D, aspect='auto', origin='lower', cmap=plt.cm.YlGnBu)
    axmatrix.set_xticks([])
    axmatrix.set_yticks([])
    
   # Plot colorbar.
    axcolor = fig.add_axes([0.91,0.1,0.02,0.6])
    pylab.colorbar(im, cax=axcolor)
    fig.show()
    #fig.savefig('dendrogram.png')
    return Z1
def PCA_custom(rho,theta,midist):
    pca=PCA()
    X=(np.vstack((theta, rho, midist)).T)
    pca.fit(X)
    print(pca.explained_variance_ratio_)
    print(pca.singular_values_)
    
    return pca

def scalesimple(rho,theta):
    X=(np.vstack((theta, rho)).T)
    X=scale(X)
    
    return X 

def scaleRobust(rho,theta):
    X=(np.vstack((theta, rho)).T)
    transformer=RobustScaler().fit(X)
    
    
    return transformer.transform(X) 

def analyzeScaling(theta, rho, midist):
    '''
    Analyze scaling of data
    
    Input:
        theta: theta values
        rho: rho values
        midist: midist values
        
    Output:

    '''
    fig,ax=plt.subplots(5,3)

    X=(np.vstack((theta, rho, midist)).T)
    transformer1=RobustScaler().fit(X)
    robust=transformer1.transform(X)
    scaler = StandardScaler().fit(X)

    normal=scaler.transform(X)
    
    pt = PowerTransformer(method='yeo-johnson', standardize=False).fit(X)
    power1=pt.transform(X)
    
    
    ax[0,0].hist(theta, bins='fd', color='b')
    ax[0,1].hist(rho, bins='fd', color='y')
    ax[0,2].hist(midist, bins='fd', color='r')


    ax[0,0].set_title("Unprocessed Theta")
    ax[0,1].set_title("Unprocessed Rho")
    ax[0,2].set_title("Unprocessed Offset")
    
    for i,c in zip( range(3), ['b', 'y', 'r']):
        ax[1,i].hist(robust[:,i], color=c)
        ax[2,i].hist(normal[:,i], color=c)
        ax[3,i].hist(power1[:,i], color=c)

    ax[1,0].set_title("Robust Theta")
    ax[1,1].set_title("Robust Rho")
    ax[1,2].set_title("Robust Offset")
    
    ax[2,0].set_title("Standard Theta")
    ax[2,1].set_title("Standard Rho")
    ax[2,2].set_title("Standard Offset")
    
    ax[3,0].set_title("Yeo-Johnson Theta")
    ax[3,1].set_title("Yeo-Johnson Rho")
    ax[3,2].set_title("Yeo-Johnson Offset")
    
    plt.tight_layout()

def intersection(df):
    inters=np.empty([len(df), len(df),2])
    for i in range(len(df)):
        x1=df['Xstart'].iloc[i]
        x2=df['Xend'].iloc[i]
        y1=df['Ystart'].iloc[i]
        y2=df['Yend'].iloc[i]
        for j in range(len(df)):
            if j == i :
                continue
            x3=df['Xstart'].iloc[j]
            x4=df['Xend'].iloc[j]
            y3=df['Ystart'].iloc[j]
            y4=df['Yend'].iloc[j]
            
            d=(x1-x2)*(y3-y4)-(y1-y2)*(x3-x4)+.0000001
            px=(x1*y2-y1*x2)*(x3-x4)-(x1-x2)*(x3*y4-y3*x4)
            py=(x1*y2-y1*x2)*(y3-y4)-(y1-y2)*(x3*y4-y3*x4)
            inters[i,j,0]=px/d
            inters[i,j,1]=py/d
    return inters

def HT_DB(dikeset, eps, samp):
    X=scalesimple(dikeset['rho'], dikeset['theta'])
    clustering = DBSCAN(eps=eps, min_samples=samp).fit(X)
    dikeset['Labels']=clustering
    return dikeset

def HT_AGG(dikeset,d):
    X=scalesimple(dikeset['rho'], dikeset['theta'])
    clusters = AGG(n_clusters=None, distance_threshold=d, linkage='complete').fit(X)
    
    return clusters

def HT_AGG_custom(dikeset,threshold, metric, dimensions=2, linkage='average'):
    '''
    Agglomerative clustering with custom metric on Hough transform data

    Input:
        dikeset: dataframe with Hough transform data
        metric: metric to use for clustering
        dimensions (optional): number of dimensions to use for clustering
        linkage (optional) : 
    Output: 
        dikeset: dataframe with cluster labels
        clusters: fitted AgglomeratieClustering object 

    '''
    t,r=whichForm(dikeset)

    #scale m unit values
    dikeset['ScaledRho']=dikeset[r].values - dikeset[r].mean()
    dikeset['ScaledPerpOffsetDist']=dikeset['PerpOffsetDist'].values-dikeset['PerpOffsetDist'].mean()
    # create X vector with theta and rho and midist
    if dimensions == 2:
        X=(np.vstack((dikeset[t], dikeset['ScaledRho'])).T)
    elif dimensions == 3:
        X=(np.vstack((dikeset[t], dikeset['ScaledRho'], dikeset['ScaledPerpOffsetDist'])).T)


    # Create a new AGG instance with the specified metric and dimensions
    M= pairwise_distances(X, metric=metric)

    # fit the data to the AGG model
    ag = AGG(n_clusters=None, affinity='precomputed', distance_threshold=threshold, linkage=linkage)
    clusters=ag.fit(M)
    dikeset['Labels']=clusters.labels_
    
    return dikeset, clusters

def AGGfull(dikeset):
    theta, rho, xc, yc= AKH_HT(dikeset)
    dikeset['rho']=rho
    dikeset['theta']=theta
    
    trange=2 
    rrange=5000 
    
    stdT=np.std(theta)
    stdRho=np.std(rho)
    d2=trange/stdT
    
    clustering=HT_AGG(dikeset,d2)
    dikeset['Labels']=clustering.labels_
    lines,IC=examineClusters(dikeset)
    
    #generate a kmz
    colorsSegments=labelcolors(dikeset['Labels'])
    colorsDikes=labelcolors(lines['Label'])
    
    #fig,ax=plotbyAngle(dikeset, lines, 20)
    
    ## Length and Width plots
    f,a=plt.subplots(2)
    a[0].hist(lines['R_Length'], bins=np.arange(0,100000,5000))
    a[0].set_xlabel('Dike Cluster Length')
    
    a[1].hist(lines['R_Width'],bins=np.arange(0,2000,100))
    a[1].set_xlabel('Dike Cluster Width')
    plt.tight_layout()
    
    f2,a2=plt.subplots(3)
    a2[0].scatter(lines['R_Length'], lines['R_Width'])
    a2[0].set_xlabel('Length')
    a2[0].set_ylabel('Width')
    
    a2[2].scatter(lines['Size'], lines['R_Width'])
    a2[2].set_xlabel('Size')
    a2[2].set_ylabel('Width')
    
    a2[1].scatter(lines['Size'], lines['R_Length'])
    a2[1].set_xlabel('Size')
    a2[1].set_ylabel('Length')
    plt.tight_layout()
    return dikeset, lines

def fullTree(model, **kwargs): 
    # Create linkage matrix and then plot the dendrogram
    
    # create the counts of samples under each node
    counts = np.zeros(model.children_.shape[0])
    n_samples = len(model.labels_)
    for i, merge in enumerate(model.children_):
        current_count = 0
        for child_idx in merge:
            if child_idx < n_samples:
                current_count += 1  # leaf node
            else:
                current_count += counts[child_idx - n_samples]
        counts[i] = current_count
    
    linkage_matrix = np.column_stack([model.children_, model.distances_,
                                      counts]).astype(float)
    
    # Plot the corresponding dendrogram
    dendrogram(linkage_matrix, **kwargs)

