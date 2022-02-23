#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 31 16:17:22 2022

@author: akh
"""
from matplotlib import pyplot as plt
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.spatial.distance import pdist, squareform
import numpy as np


np.random.seed(4711)  
X = np.random.uniform(-10,10,size=(10,2))
labels=[i for i in range(len(X))]
fig,ax=plt.subplots()
ax.scatter(X[:,0], X[:,1])

dist=squareform(pdist(X))

Y = linkage(squareform(dist), method='complete')
Z1 =dendrogram(Y, labels=labels, orientation='left')

