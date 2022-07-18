#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May  4 10:17:24 2022

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
from examineMod import examineClusters, checkoutCluster, TopHTSection
import os
import labellines
from scipy import stats
import seaborn as sns
from scipy.stats import lognorm


def plotLargeClusters(dikeset,lines,n, path, overide=False):
    if type(n) is int:
        n=float(n)
    largelines=lines.loc[np.greater(lines['Size'].values, n)]
    ImgPath=path[:-4]+"Images/"
    
    isExist = os.path.exists(ImgPath)
    
    if not isExist:
        os.makedirs(ImgPath)
        print("The new directory", ImgPath, " is created!")
    
    for i in largelines['Label'].values:
        ImgName=ImgPath+"Label_"+str(int(i))+".pdf"
        isExist = os.path.exists(ImgName)
        if not isExist or overide:
            fig,ax=checkoutCluster(dikeset, i)
            fig.savefig(ImgName, transparent=True)
            plt.close(fig)


def Histograms(dikeset,lines, path, maxL=150000):
    mosaic="AB\nCB\nDB"
    fig=plt.figure()
    fig.set_size_inches( 12,6)
    ax = fig.subplot_mosaic(mosaic)
    #identify_axes(ax_dict)

    ax['A']=trueDikeLength(lines, dikeset, maxL, axs=ax['A'])
    # shape, loc, scale=lognorm.fit(lines['R_Length'].values)
    # x=np.linspace(0,maxL)
    
    # aa=ax['A'].twinx()
    # aa.plot(x,lognorm.pdf(x,shape, loc, scale), 'r-.')
    # ax['A'].set_ylabel('Counts', color='b')
    # aa.set_ylabel('Normalized Probability', color='r')
    
    # shape, loc, scale=lognorm.fit(dikeset['seg_length'].values)
    # aa.plot(x,lognorm.pdf(x,shape, loc, scale), 'g-.')
    
    # lo = Labeloffset(aa, label='Normalized Probability', axis="y")

    
    ax['C'].hist(lines['R_Width'], bins=np.arange(0,5000,500))
    
    ax['C'].text(.60,.80,'Dike median:'+str(round(np.median(lines['R_Width'].values),0)), transform=ax['C'].transAxes)
    ax['C'].text( .60, .65, 'Dike STD:'+str(round(lines['R_Width'].std(),0)),transform=ax['C'].transAxes)
    mask=np.greater(lines['Size'].values, 3)
    
    AR=lines.loc[mask]['R_Length'].values/(lines['R_Width'].loc[mask].values+1)
    
    ax['D'].hist(AR, bins=np.arange(0,80,10))
    #sns.histplot(AR, bins=np.arange(0,80,10), ax=ax['D'])
    
    ax['D'].text(.70,.80,'Dike median:'+str(round(np.median(AR),2)), transform=ax['D'].transAxes)
    ax['D'].text( .70, 65, 'Dike STD:'+str(round(np.std(AR),2)),transform=ax['D'].transAxes)
    print('Dike STD:'+str(round(np.std(AR),2)))
    ax['B'].scatter(lines[mask]['R_Width'], lines.loc[mask]['R_Length'].values, s=2*lines.loc[mask]['Size'].values**1.5, alpha=0.6, edgecolor='grey')
    
    ax['C'].set_ylabel("Counts")
    ax['C'].set_xlabel("Width (m)")
    
    ax['D'].set_ylabel("Counts")
    ax['D'].set_xlabel("Aspect Ratio (L/W)")
    
    ax['B'].set_ylabel("Length (m)")
    ax['B'].set_xlabel("Width (m))")
    
    x=np.linspace(0,6000,100)
    for AR in [10,25,50]:
        ax['B'].plot(x, x*AR, 'k-.', alpha=0.7, label=str(AR))
    
    labellines.labelLines(ax['B'].get_lines())
    ax['B'].set_ylim([-5,300000])
    ax['B'].set_xlim([-5,6000])
    ImgPath=path[:-4]+"Images/"
    ImgName=ImgPath+"Histograms"+".pdf"
    fig.savefig(ImgName, transparent=True)
    plt.close(fig)

def HT3(dikeset,lines,path):
    
    fig,ax=plt.subplots(1,3)
    fig.set_size_inches( 15,4)
    DotsHT(fig, ax[0], dikeset, ColorBy="seg_length", title="Unlinked Segments")
    DotsHT(fig, ax[1], lines, ColorBy='R_Length', title="Linked Segments")
    mask=np.greater(lines['Size'].values, 4)
    DotsHT(fig, ax[2], lines.loc[mask], ColorBy='R_Length', title="5 or more Linked Segments")
    
    ax[0].text(.10,.10,'n='+str(len(dikeset)), transform=ax[0].transAxes)
    ax[1].text(.10,.10,'n='+str(len(lines)), transform=ax[1].transAxes)
    ax[2].text(.10,.10,'n='+str(len(lines.loc[mask])), transform=ax[2].transAxes)
    plt.tight_layout()
    ImgPath=path[:-4]+"Images/"
    ImgName=ImgPath+"HT3"+".pdf"
    fig.savefig(ImgName, transparent=True)
    plt.close(fig)
    
def Toplines(lines, dikeset, path):
    T, fig, ax=TopHTSection(lines, dikeset, 10000, 4)
    
    ImgPath=path[:-4]+"Images/"
    ImgName=ImgPath+"TopLines"+".pdf"
    fig.savefig(ImgName, transparent=True)
    plt.close(fig)
    
def PairGrid(lines,cols, path,n=2):
    mask=np.greater(lines['Size'].values, n)
    temp=lines.loc[mask][cols]
    ImgPath=path[:-4]+"Images/"
    ImgName=ImgPath+"PairPlot"+".pdf"
    
    g = sns.PairGrid(temp, diag_sharey=False, hue="Size", corner=True)
    g.map_diag(sns.histplot,  hue=None, color=".3")
    g.map_offdiag(sns.scatterplot)
    g.add_legend()
    

    g.savefig(ImgName, transparent=True)
    #plt.close(g)


def allFigures(path1,path2):
    dikeset=pd.read_csv(path1)
    lines=pd.read_csv(path2)
    cols=['Size', 'Overlap', "MaxSegNNDist", "MinSegNNDist"]
    
    ImgPath=path2[:-4]+"_Images/"
    sns.reset_orig()
    isExist = os.path.exists(ImgPath)
    if not isExist:
        os.makedirs(ImgPath)
        print("The new directory", ImgPath, " is created!")
    sns.reset_orig()
    
    g=sns.jointplot(x=lines['Size'].loc[lines['Size']>2], y=lines['R_Length'].loc[lines['Size']>2], alpha=0.6)
    
    
    slope, intercept, r_value, p_value, std_err = stats.linregress(lines['Size'].loc[lines['Size']>2].values, lines['R_Length'].loc[lines['Size']>2].values)
    x=np.linspace(2, lines['Size'].max())
    y=slope*x+intercept 
    
    g.ax_joint.plot(x,y, 'g-.', alpha=0.6)
    g.ax_joint.annotate(f'$r = {r_value:.3f}, p = {p_value:.3f}$',
                    xy=(0.1, 0.9), xycoords='axes fraction',
                    ha='left', va='center',
                    bbox={'boxstyle': 'round', 'fc': 'powderblue', 'ec': 'navy'})
    
    g2=sns.jointplot(x=lines['Size'].loc[lines['Size']>5], y=lines['R_Length'].loc[lines['Size']>5], alpha=0.6)
    g.savefig(ImgPath+"SizebyLength2.pdf", transparent=True)
    g2.savefig(ImgPath+"SizebyLength5.pdf", transparent=True)
    
    
    Histograms(dikeset, lines, path2)
    plotLargeClusters(dikeset, lines, 4, path2, overide=False)
    HT3(dikeset, lines, path2)
    Toplines(lines, path2)
    PairGrid(lines, cols, path2)
    sns.reset_orig()
    
    fig,ax=DotsLinesHist(lines, 5000, 4)
    fig.set_size_inches(10,5)
    
    plt.tight_layout()
    ImgName=ImgPath+"DotsLinesHist"+".pdf"
    fig.savefig(ImgName, transparent=True)
    plt.close(fig)
    
# l=['dikedata/deccandata/Central_3_2000_LINKED.csv', 
#    'dikedata/deccandata/NarmadaTapi_3_2000_LINKED.csv',
#    'dikedata/deccandata/Saurashtra_3_2000_LINKED.csv',
#    'dikedata/crb/CJDS_Lines_3_500_March11.csv',
#    'dikedata/spanish peaks/SpanishPeaks_3857_LINKED_2_2000.csv']

# d=['dikedata/deccandata/Central_preprocesed.csv',
#    'dikedata/deccandata/NarmadaTapi_preprocesed.csv',
#    'dikedata/deccandata/Saurashtra_preprocesed.csv',
#    'dikedata/crb/CJDS_FebStraightened.csv',
#    'dikedata/spanish peaks/SpanishPeaks_3857_preprocessed.csv'
#    ]

# for lines,dikeset in zip(l,d):
#     allFigures(dikeset, lines)