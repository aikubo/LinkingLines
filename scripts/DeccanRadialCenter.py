#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 25 12:52:05 2022

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

lines=pd.read_csv('dikedata/deccandata/AllDeccanLinked_euc_18_10_2022.csv')
dikeset=pd.read_csv('/home/akh/myprojects/Linking-and-Clustering-Dikes/dikedata/deccandata/AllDeccan_PreProcessed.csv')
xc,yc=HT_center(dikeset)
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

# mask2=(mlines['Ymid']>2.38e6) & (mlines['Ymid']<2.42e6) 


# Centers24=RadialFit(mlines[mask2], ThetaRange=[-50,50])



# ax['B'].plot(xdata, CenterFunc(xdata, Centers24['Center'][0][0], Centers24['Center'][0][1], xc, yc), 'r.-',
#           label='fit: xr=%5.3f, yr=%5.3f' % tuple(Centers24['Center'][0]), linewidth=5)


# mask2=mlines['Ymid']<2.3e6

# Centers23=RadialFit(mlines[mask2])

# ax['B'].plot(xdata, CenterFunc(xdata, Centers23['Center'][0][0], Centers23['Center'][0][1], xc, yc), 'g.-',
#           label='fit: xr=%5.3f, yr=%5.3f' % tuple(Centers23['Center'][0]), linewidth=5)


# Close23=NearCenters(mlines, Centers23)


# Close24=NearCenters(mlines, Centers24)
# plotlines(lines, 'grey', ax['C'], alpha=0.05)
# plotlines(Close23, 'g', ax['C'], alpha=0.6)
# plotlines(Close24, 'r', ax['C'], alpha=0.6)

# ax['C'].plot(  Centers23['Center'][0][0], Centers23['Center'][0][1], '*g', markersize=10, markeredgecolor='k')

# ax['C'].plot(  Centers24['Center'][0][0], Centers24['Center'][0][1], '*r', markersize=10, markeredgecolor='k')
# identify_axes(ax, fontsize=36)
# fig.savefig("/home/akh/myprojects/Linking-and-Clustering-Dikes/dikedata/deccandata/AllDeccanLinked_24_08_2022Images/AllDeccanRadialFitLinked.png")
# fig.savefig("/home/akh/myprojects/Linking-and-Clustering-Dikes/dikedata/deccandata/AllDeccanLinked_24_08_2022Images/AllDeccanRadialFitLinked.pdf", dpi=600)


# # mosaic="AB\nCC\nCC"
# # fig=plt.figure()
# # fig.set_size_inches( 8,12)
# # ax = fig.subplot_mosaic(mosaic)


# # DotsHT(fig, ax['A'], lines, ColorBy='Ymid', title="All Linked Segments")
# mask=lines['TrustFilter']==1
# mlines=lines[mask]
# # DotsHT(fig, ax['B'], lines.loc[mask], ColorBy='Ymid', title="Highest Confidence Linked Segments")


# CentersAll=RadialFit(lines[mask])
# xdata=np.linspace(-90,90,100)
# # ax['B'].plot(xdata, CenterFunc(xdata, CentersAll['Center'][0][0], CentersAll['Center'][0][1], xc, yc), 'k.-',
# #           label='fit: xr=%5.3f, yr=%5.3f' % tuple(CentersAll['Center'][0]), linewidth=5)
# # ax['C'].plot(  CentersAll['Center'][0][0], CentersAll['Center'][0][1], '*k', markersize=10)

# mask2=(mlines['Ymid']>2.38e6) & (mlines['Ymid']<2.42e6) 


# Centers24=RadialFit(mlines[mask2], ThetaRange=[-50,50])



# # ax['B'].plot(xdata, CenterFunc(xdata, Centers24['Center'][0][0], Centers24['Center'][0][1], xc, yc), 'r.-',
# #           label='fit: xr=%5.3f, yr=%5.3f' % tuple(Centers24['Center'][0]), linewidth=5)


# mask2=mlines['Ymid']<2.3e6

# Centers23=RadialFit(mlines[mask2])

# # ax['B'].plot(xdata, CenterFunc(xdata, Centers23['Center'][0][0], Centers23['Center'][0][1], xc, yc), 'g.-',
# #           label='fit: xr=%5.3f, yr=%5.3f' % tuple(Centers23['Center'][0]), linewidth=5)


# Close23=NearCenters(lines, Centers23)


# Close24=NearCenters(lines, Centers24)
# plotlines(lines, 'grey', ax['C'], alpha=0.05)
# plotlines(Close23, 'g', ax['C'], alpha=0.6)
# plotlines(Close24, 'r', ax['C'], alpha=0.6)

# ax['C'].plot(  Centers23['Center'][0][0], Centers23['Center'][0][1], '*g', markersize=10, markeredgecolor='k')

# ax['C'].plot(  Centers24['Center'][0][0], Centers24['Center'][0][1], '*r', markersize=10, markeredgecolor='k')
# identify_axes(ax, fontsize=36)
# fig.savefig("/home/akh/myprojects/Linking-and-Clustering-Dikes/dikedata/deccandata/AllDeccanLinked_24_08_2022Images/AllDeccanRadialFitLinked_allDikesMatch.png")
# fig.savefig("/home/akh/myprojects/Linking-and-Clustering-Dikes/dikedata/deccandata/AllDeccanLinked_24_08_2022Images/AllDeccanRadialFitLinked_allDikesMatch.pdf", dpi=600)

fig,ax = plt.subplots()
fig.set_size_inches(12,6)
DotsHT(fig, ax, dikeset, ColorBy='Ymid', title="All Segments")


fig = plt.figure(tight_layout=True)
fig.set_size_inches(12,6)
gs = gridspec.GridSpec(6, 6)

cart = fig.add_subplot(gs[:, 1:])
a=fig.add_subplot(gs[0,0])
DotsHT(fig, a, dikeset, ColorBy='Ymid', title="All Segments")

mask=lines['TrustFilter']==1
mlines=lines[mask]
mask2=mlines['Ymid']<2.3e6
b=fig.add_subplot(gs[1,0])
CentersB=RadialFit(mlines[mask2])
CloseB=NearCenters(lines, CentersB)
DotsHT(fig, b, CloseB, color='r', Cbar=False)
plotlines(CloseB, 'r', cart)
cart.plot(  CentersB['Center'][0][0], CentersB['Center'][0][1], '*r', markersize=10, markeredgecolor='k')


mask2=(mlines['Ymid']>2.38e6) & (mlines['Ymid']<2.42e6) 
c=fig.add_subplot(gs[2,0])
CentersC=RadialFit(mlines[mask2])
CloseC=NearCenters(lines, CentersC)
DotsHT(fig, c, CloseC, color='g', Cbar=False)
plotlines(CloseC, 'g', cart)
cart.plot(  CentersC['Center'][0][0], CentersC['Center'][0][1], '*g', markersize=10, markeredgecolor='k')

mask2=(dikeset['Ymid']>2.6) 
d=fig.add_subplot(gs[3,0])
CentersD=RadialFit(dikeset[mask2])
CloseD=NearCenters(lines, CentersD)
DotsHT(fig, d, CloseD, color='b', Cbar=False)
plotlines(CloseD, 'b', cart)
cart.plot(  CentersD['Center'][0][0], CentersD['Center'][0][1], '*b', markersize=10, markeredgecolor='k')

