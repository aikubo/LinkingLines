# LinkingLines Package
 # Written by aikubo
 # Version: 2.1.0
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  1 13:12:50 2021

@author: akh
"""

import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import dendrogram
from .HT import rotateData, HoughTransform
#from examineMod import *
from .PrePostProcess import whichForm
from scipy.spatial.distance import pdist, squareform
import scipy.cluster.hierarchy as sch
import matplotlib.pyplot as plt


def AggCluster(dikeset, dtheta, drho, dimensions=2, linkage='complete', rotate=False, metric='Euclidean'):
    """
    Agglomerative clustering with custom metric on Hough transform data.

    Parameters
    ----------
        dikeset : pandas DataFrame
            DataFrame with Hough transform data.
        dtheta : float
            Scaling factor for theta.
        drho : float
            Scaling factor for rho.
        dimensions : int, optional
            Number of dimensions to use for clustering (default is 2).
        linkage : str, optional
            Linkage method for hierarchical clustering (default is 'complete').
        rotate : bool, optional
            Whether to rotate the dataset (default is False).
        metric : str, optional
            Metric to use for clustering (default is 'Euclidean').

    Returns
    -------
        dikeset : pandas DataFrame
            DataFrame with cluster labels.
        Z : numpy ndarray
            The hierarchical clustering linkage matrix.
    """
    # if 'theta' is not dikeset.columns:
    # do the Hough transform
    if 'theta' not in dikeset.columns:
        print("Hough transform not found, performing Hough transform")
        dikeset, xc, yc=HoughTransform(dikeset) 

    t,r=whichForm(dikeset)
    angle=np.median(abs(dikeset[t]))-20

    if rotate:
        print("rotating dataset by", angle)
        dikeset=rotateData(dikeset,angle)

    #scale m unit values
    dikeset['ScaledRho']=dikeset[r].values - dikeset[r].mean()

    # create X vector with theta and rho and midist
    X=(np.vstack((dikeset[t], dikeset[r])).T)

    threshold=1



    X=X/[dtheta, drho]


    M= pdist(X, metric)

    if linkage=='complete':
        Z=sch.complete(M)
    elif linkage=='average':
        Z= sch.average(M)
    elif linkage=='single':
        Z=sch.single(M)


    labels=sch.fcluster(Z, t=threshold, criterion='distance')
    #rootnode, nodelist=sch.to_tree(Z)
    dikeset['Labels']=labels

    #unrotate
    if rotate:
        dikeset=rotateData(dikeset,-1*angle)


    return dikeset, Z #, rootnode, nodelist



def fullTree(model, **kwargs):
    """
    Generate and plot a full dendrogram for hierarchical clustering results.

    Parameters
    ----------
        model : sklearn.cluster.AgglomerativeClustering
            Fitted AgglomerativeClustering model.
        **kwargs : dict
            Additional keyword arguments to be passed to the dendrogram function.

    Returns
    -------
        None

    Note:
        This function generates and plots a full dendrogram showing the hierarchical clustering of data based on the provided AgglomerativeClustering model.
        It calculates counts of samples under each node and creates a linkage matrix for the dendrogram.

    Example:
        fullTree(clustering_model, color_threshold=0.5)
    """
    # Rest of the function code...

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

