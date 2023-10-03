from synthetic import *
import numpy as np
from plotmod import plotlines, DotsHT
import pandas as pd
from htMOD import AKH_HT as HT
import matplotlib.pyplot as plt
from sklearn.cluster import DBSCAN
import matplotlib.lines as mlines
from sklearn.preprocessing import scale
from scipy.spatial.distance import pdist, squareform
from htMOD import MidtoPerpDistance
from examineMod import *
from matplotlib import cm
import scipy
import pylab
import scipy.cluster.hierarchy as sch
from scipy.spatial.distance import squareform
from sklearn.cluster import AgglomerativeClustering as AGG
from scipy.cluster.hierarchy import dendrogram
from PrePostProcess import *


def plotDendro(dist1, labels, title):

    # https://stackoverflow.com/questions/2982929/plotting-results-of-hierarchical-clustering-ontop-of-a-matrix-of-data-in-python
    D = dist1
    condensedD = squareform(D)

    # Compute and plot first dendrogram.
    fig = pylab.figure(figsize=(8, 8))

    ax1 = fig.add_axes([0.09, 0.1, 0.15, 0.6])
    ax1.set_title(title)
    Y = sch.linkage(condensedD, method='complete')
    Z1 = sch.dendrogram(Y, labels=labels, orientation='left', no_plot=True)
    # ax1.set_xticks()
    # ax1.set_yticks([])

    # Compute and plot second dendrogram.
    ax2 = fig.add_axes([0.3, 0.76, 0.6, 0.2])
    Y = sch.linkage(condensedD, method='complete')
    Z2 = sch.dendrogram(Y,  labels=labels, no_plot=True)
    # ax2.set_xticks(labels[Z1['leaves']])
    # ax2.set_yticks([])

    # Plot distance matrix.
    # add axis(left, bottom, width, height)
    axmatrix = fig.add_axes([0.3, 0.1, 0.6, 0.6])
    idx1 = Z1['leaves']
    idx2 = Z2['leaves']
    D = D[idx1, :]
    D = D[:, idx2]
   # im = axmatrix.matshow(D, aspect='auto', origin='lower', cmap=plt.cm.YlGnBu)
    axmatrix.set_xticks([])
    axmatrix.set_yticks([])

   # Plot colorbar.
    #axcolor = fig.add_axes([0.91, 0.1, 0.02, 0.6])
    #pylab.colorbar(im, cax=axcolor)
    #fig.show()
    # fig.savefig('dendrogram.png')
    return Z1

from clusterMod import CyclicEuclideanScaled

from PrePostProcess import completePreProcess, whichForm

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
    
    # for ys, color in zip(dcoord, c):
    #     #ax[0].plot(xs, ys, color)
        
    #     birth=np.array([ys[0]+1, ys[3]+1])
    #     death=np.array([ys[1]+1, ys[2]+1])
    #     ax[0].plot(birth,death, "*", color=color)
    #     p=np.append(birth-death)
    
    return fig, ax, Z1
    

# length, angle, angleSTD, rho, rhoSTD, ndikes=100, CartRange=300000
df1 = makeLinear2(5000, 20, 5, 10000, 1000, 50, CartRange=100000)
df2 = makeLinear2(5000, 21, 2, 20000, 1000, 50, CartRange=100000)
df3 = makeLinear2(5000, 19, 5, -20000, 1000, 50, CartRange=100000)

df = addSwarms([df1, df2, df3])
dfFrag = fragmentDikes(df)
dikes= completePreProcess(dfFrag)


theta, rho, xc, yc= HT(dikes, xc=0, yc=0)
dikes['theta']=theta
dikes['rho']=rho
dikes= MidtoPerpDistance(dikes, xc, yc)
midist=dikes['PerpOffsetDist'].values

X3D = (np.vstack( (theta, rho, midist) )).T
X2D = (np.vstack( (theta, rho-np.mean(rho))) ).T

fig,ax=plt.subplots(2); plotlines(dikes, 'k', ax[0], ColorBy='label')
DotsHT(fig, ax[1], dikes, ColorBy="label")
ax[0].axis('equal')

dist2=squareform(pdist(X2D, CyclicEuclidean))

from clusterMod import HT_AGG_custom

dikes2,clusters=HT_AGG_custom(dikes,5000, 2, 500)

lines, IC=examineClusters(dikes2)
fig, ax, Z=persitance(dikes)



df1 = makeLinear2(5000, 20, 5, 10000, 1000, 50, CartRange=100000)
df2 = makeLinear2(5000, 50, 2, 20000, 1000, 50, CartRange=100000)
df3 = makeLinear2(5000, 1, 5, -20000, 1000, 50, CartRange=100000)

df = addSwarms([df1, df2, df3])
dfFrag = fragmentDikes(df)
dikes= completePreProcess(dfFrag)


theta, rho, xc, yc= HT(dikes, xc=0, yc=0)
dikes['theta']=theta
dikes['rho']=rho
dikes= MidtoPerpDistance(dikes, xc, yc)
midist=dikes['PerpOffsetDist'].values

X3D = (np.vstack( (theta, rho, midist) )).T
X2D = (np.vstack( (theta, rho-np.mean(rho))) ).T

fig,ax=plt.subplots(2); plotlines(dikes, 'k', ax[0], ColorBy='label')
DotsHT(fig, ax[1], dikes, ColorBy="label")
ax[0].axis('equal')

dist2=squareform(pdist(X2D, CyclicEuclidean))

from clusterMod import HT_AGG_custom

# dikes2,clusters=HT_AGG_custom(dikes,5000, CyclicEuclidean)

# lines, IC=examineClusters(dikes2)
fig, ax, Z=persitance(dikes)



df = RandomDikes(100)
df['label']=np.ones(len(df))
dfFrag = fragmentDikes(df)
dikes= completePreProcess(dfFrag)


theta, rho, xc, yc= HT(dikes, xc=0, yc=0)
dikes['theta']=theta
dikes['rho']=rho
dikes= MidtoPerpDistance(dikes, xc, yc)
midist=dikes['PerpOffsetDist'].values

X3D = (np.vstack( (theta, rho, midist) )).T
X2D = (np.vstack( (theta, rho-np.mean(rho))) ).T

fig,ax=plt.subplots(2); plotlines(dikes, 'k', ax[0], ColorBy='label')
DotsHT(fig, ax[1], dikes, ColorBy="label")
ax[0].axis('equal')

dist2=squareform(pdist(X2D, CyclicEuclidean))

from clusterMod import HT_AGG_custom

# dikes2,clusters=HT_AGG_custom(dikes,5000, CyclicEuclidean)

# lines, IC=examineClusters(dikes2)
fig, ax, Z=persitance(dikes)

