#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 26 19:23:13 2022

@author: akh
"""
from PrePostProcess import * 
import pandas as pd 
import numpy as np 
from htMOD import AKH_HT as HT 
from htMOD import HT_center, MidtoPerpDistance
from plotmod import *
import matplotlib.pyplot as plt 
from clusterMod import HT_AGG_custom
import scipy.cluster.hierarchy as sch
from examineMod import examineClusters, checkoutCluster, TopHTSection, CheckoutBy
import os
import labellines
from scipy import stats
import seaborn as sns
from scipy.stats import lognorm
import time
import dateutil.parser
from datetime import datetime 
from datetime import timedelta
from fitRadialCenters import RadialFit, CenterFunc, NearCenters

lines=pd.read_csv('dikedata/crb/AllCRBLinked_24_08_2022.csv')
dikeset=pd.read_csv('dikedata/crb/allCRB_dikes_PreProcessed.csv')
xc,yc=HT_center(dikeset)

fig,ax=plt.subplots()
ax.scatter(lines['Ymid'], lines['AvgTheta'], alpha=0.6)
ax.set_ylabel('Average Theta')
ax.set_xlabel('Y Line Midpoint (Latitude)')


mosaic="AB\nCC\nCC"
fig=plt.figure()
fig.set_size_inches( 8,12)
ax = fig.subplot_mosaic(mosaic)


DotsHT(fig, ax['A'], lines, ColorBy='Ymid', title="All Linked Segments")
mask=lines['TrustFilter']==1

mask= (lines['Size']>2) #& (lines['MaxSegNNDist']<0.7)

mlines=lines[mask]
ax['B'].scatter(lines['AvgTheta'], lines['AvgRho'], c='white')
#DotsHT(fig, ax['B'], lines.loc[mask], ColorBy='Ymid', title="Highest Confidence Linked Segments")



#ax['B'].set_ylim((dikeset['rho'].min(),dikeset['rho'].max() ))


# CentersAll=RadialFit(lines[mask])
# 
# ax['B'].plot(xdata, CenterFunc(xdata, CentersAll['Center'][0][0], CentersAll['Center'][0][1], xc, yc), 'k.-',
#           label='fit: xr=%5.3f, yr=%5.3f' % tuple(CentersAll['Center'][0]), linewidth=5)
# ax['C'].plot(  CentersAll['Center'][0][0], CentersAll['Center'][0][1], '*k', markersize=10, alpha=0.2)

mask2=(mlines['Ymid']>4.98e6) & (mlines['Ymid']<5.02e6) 
xdata=np.linspace(-90,90,100)

Centers24=RadialFit(mlines[mask2])



ax['B'].plot(xdata, CenterFunc(xdata, Centers24['Center'][0][0], Centers24['Center'][0][1], xc, yc), 'r.-',
          label='fit: xr=%5.3f, yr=%5.3f' % tuple(Centers24['Center'][0]), linewidth=5)


mask2=(mlines['Ymid']<4.95e6) & (mlines['Ymid']>4.87e6)

Centers23=RadialFit(mlines[mask2])

ax['B'].plot(xdata, CenterFunc(xdata, Centers23['Center'][0][0], Centers23['Center'][0][1], xc, yc), 'g.-',
          label='fit: xr=%5.3f, yr=%5.3f' % tuple(Centers23['Center'][0]), linewidth=5)


Close23=NearCenters(mlines, Centers23)


Close24=NearCenters(mlines, Centers24)
#ax['B'].scatter(Close23['AvgTheta'], Close23['AvgRho'], c='g', alpha=0.6, edgecolor='k')
#ax['B'].scatter(Close24['AvgTheta'], Close24['AvgRho'], c='r', alpha=0.6, edgecolor='k')
DotsHT(fig, ax['B'], lines, ColorBy='Ymid')

plotlines(lines, 'grey', ax['C'], alpha=0.05)

plotlines(Close24, 'r', ax['C'], alpha=0.6)
plotlines(Close23, 'g', ax['C'], alpha=0.6)

ax['C'].plot(  Centers23['Center'][0][0], Centers23['Center'][0][1], '*g', markersize=10, markeredgecolor='k')

ax['C'].plot(  Centers24['Center'][0][0], Centers24['Center'][0][1], '*r', markersize=10, markeredgecolor='k')
identify_axes(ax, fontsize=36)


mask3=  (mlines['Ymid']>4.706e6) & (mlines['Ymid']<4.76e6)
#(mlines['Ymid']>5.06e6) & (mlines['Ymid']<5.13e6) 
#^ v close to red star

Centers3=RadialFit(mlines[mask3])
Close3=NearCenters(mlines, Centers3)

ax['C'].plot(  Centers3['Center'][0][0], Centers3['Center'][0][1], '*b', markersize=10, markeredgecolor='k')

plotlines(Close3, 'b', ax['C'], alpha=0.6)
ax['B'].plot(xdata, CenterFunc(xdata, Centers3['Center'][0][0], Centers3['Center'][0][1], xc, yc), 'b.-',
          label='fit: xr=%5.3f, yr=%5.3f' % tuple(Centers3['Center'][0]), linewidth=5)



#fig.savefig("/home/akh/myprojects/Linking-and-Clustering-Dikes/dikedata/deccandata/AllDeccanLinked_24_08_2022Images/AllDeccanRadialFitLinked.png")
#fig.savefig("/home/akh/myprojects/Linking-and-Clustering-Dikes/dikedata/deccandata/AllDeccanLinked_24_08_2022Images/AllDeccanRadialFitLinked.pdf", dpi=600)


# mosaic="AB\nCC\nCC"
# fig=plt.figure()
# fig.set_size_inches( 8,12)
# ax = fig.subplot_mosaic(mosaic)


# DotsHT(fig, ax['A'], lines, ColorBy='Ymid', title="All Linked Segments")
# mask=lines['TrustFilter']==1
# mlines=lines[mask]
# DotsHT(fig, ax['B'], lines.loc[mask], ColorBy='Ymid', title="Highest Confidence Linked Segments")


# CentersAll=RadialFit(lines[mask])
# xdata=np.linspace(-90,90,100)
# ax['B'].plot(xdata, CenterFunc(xdata, CentersAll['Center'][0][0], CentersAll['Center'][0][1], xc, yc), 'k.-',
#           label='fit: xr=%5.3f, yr=%5.3f' % tuple(CentersAll['Center'][0]), linewidth=5)
# ax['C'].plot(  CentersAll['Center'][0][0], CentersAll['Center'][0][1], '*k', markersize=10)

# mask2= (mlines['Ymid']>4.97e6) # & (mlines['Ymid']<2.42e6) 


# Centers24=RadialFit(mlines[mask2], ThetaRange=[-50,50])



# ax['B'].plot(xdata, CenterFunc(xdata, Centers24['Center'][0][0], Centers24['Center'][0][1], xc, yc), 'r.-',
#           label='fit: xr=%5.3f, yr=%5.3f' % tuple(Centers24['Center'][0]), linewidth=5)


# mask2=mlines['Ymid']<4.97e6

# Centers23=RadialFit(mlines[mask2])

# ax['B'].plot(xdata, CenterFunc(xdata, Centers23['Center'][0][0], Centers23['Center'][0][1], xc, yc), 'g.-',
#           label='fit: xr=%5.3f, yr=%5.3f' % tuple(Centers23['Center'][0]), linewidth=5)


# Close23=NearCenters(lines, Centers23)


# Close24=NearCenters(lines, Centers24)
# plotlines(lines, 'grey', ax['C'], alpha=0.05)
# plotlines(Close23, 'g', ax['C'], alpha=0.6)
# plotlines(Close24, 'r', ax['C'], alpha=0.6)

# ax['C'].plot(  Centers23['Center'][0][0], Centers23['Center'][0][1], '*g', markersize=10, markeredgecolor='k')

# ax['C'].plot(  Centers24['Center'][0][0], Centers24['Center'][0][1], '*r', markersize=10, markeredgecolor='k')
# identify_axes(ax, fontsize=36)


# #fig.savefig("/home/akh/myprojects/Linking-and-Clustering-Dikes/dikedata/deccandata/AllDeccanLinked_24_08_2022Images/AllDeccanRadialFitLinked_allDikesMatch.png")
# #fig.savefig("/home/akh/myprojects/Linking-and-Clustering-Dikes/dikedata/deccandata/AllDeccanLinked_24_08_2022Images/AllDeccanRadialFitLinked_allDikesMatch.pdf", dpi=600)