#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  1 13:12:50 2021

@author: akh
"""

from sklearn.preprocessing import scale
from sklearn.cluster import AgglomerativeClustering as AGG
from sklearn.cluster import DBSCAN
from sklearn.neighbors import NearestNeighbors
import numpy as np 
import pandas as pd
from scipy.cluster.hierarchy import dendrogram
from htMOD import *
from examineMod import *

def scalesimple(rho,theta):
    X=(np.vstack((theta, rho)).T)
    X=scale(X)
    
    return X 

def HT_DB(dikeset, eps, samp):
    X=scalesimple(dikeset['rho'], dikeset['theta'])
    clustering = DBSCAN(eps=eps, min_samples=samp).fit(X)
    dikeset['Labels']=clustering
    return dikeset

def HT_AGG(dikeset,d):
    X=scalesimple(dikeset['rho'], dikeset['theta'])
    clusters = AGG(n_clusters=None, distance_threshold=d).fit(X)
    
    return clusters

def HT_AGG(dikeset,d):
    X=scalesimple(dikeset['rho'], dikeset['theta'])
    clusters = AGG(n_clusters=None, distance_threshold=d).fit(X)
    
    return clusters

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