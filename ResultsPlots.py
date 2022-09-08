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

def MakeFigCheck(ImgPath, overide=False):
    
    if os.path.exists(ImgPath):
        mtime=time.ctime(os.path.getmtime(ImgPath))
        mtime=dateutil.parser.parse(mtime)
        now=datetime.now()
        
         
        tdif=mtime + timedelta(days=2) < now
    
    elif not os.path.exists(ImgPath):
        tdif=True
    return tdif
    

def plotLargeClusters(dikeset,lines,n, path, overide=False):
    if type(n) is int:
        n=float(n)
    largelines=lines.loc[np.greater(lines['Size'].values, n)]
    ImgPath=path[:-4]+"Images/LargestClusters/"
    
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
            
    col=['Overlap','MaxSegNNDist', 'nOverlapingSegments', 'MedianSegNNDist','EnEchelonAngleDiff','ThetaRange', 'Size', 'Aspect', 'R_Length', 'R_Width']
    lines=lines.assign(Aspect=lines['R_Length'].values/lines['R_Width'].values)
    
    for i in col:
        
        
        ImgName=ImgPath+"Largest_"+i+".pdf"

        if MakeFigCheck(ImgName):
            fig,ax=CheckoutBy(dikeset, lines, i)
            fig.savefig(ImgName, transparent=True)
            plt.close(fig)
        


def Histograms(dikeset,lines, path, maxL=150000):
    ImgPath=path[:-4]+"Images/"
    ImgName=ImgPath+"Histograms"+".pdf"
    if MakeFigCheck(ImgName):
        mosaic="AB\nCB\nDB"
        fig=plt.figure()
        fig.set_size_inches( 12,6)
        ax = fig.subplot_mosaic(mosaic)
        #identify_axes(ax_dict)
        rlim=lines['Rho_Threshold'][0]
        lines=lines.loc[lines['Linked']==1]
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
    
        
        ax['C'].hist(lines['R_Width'], bins=np.arange(0,rlim*2,500))
        
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
        ax['C'].set_xlabel("Packet Width (m)")
        
        ax['D'].set_ylabel("Counts")
        ax['D'].set_xlabel("Aspect Ratio (L/W)")
        
        ax['B'].set_ylabel("Length (m)")
        ax['B'].set_xlabel("Width (m))")
        ax['B'].set_xlim((0,2*rlim))
        
        x=np.linspace(0,6000,100)
        for AR in [10,25,50]:
            ax['B'].plot(x, x*AR, 'k-.', alpha=0.7, label=str(AR))
        
        labellines.labelLines(ax['B'].get_lines())
        ax['B'].set_ylim([-5,maxL])

        ImgPath=path[:-4]+"Images/"
        
        fig.savefig(ImgName, transparent=True)
        plt.close(fig)

def HT3(dikeset,lines,path):
    ImgPath=path[:-4]+"Images/"
    ImgName=ImgPath+"HT3"+".pdf"
    if MakeFigCheck(ImgName):
        fig,ax=plt.subplots(1,3)
        fig.set_size_inches( 16,6)
        DotsHT(fig, ax[0], dikeset, ColorBy=None, title="Unlinked Segments")
        DotsHT(fig, ax[1], lines, ColorBy=None, title="All Linked Segments")
        mask=lines['TrustFilter']==1
        DotsHT(fig, ax[2], lines.loc[mask], ColorBy=None, title="Highed Confidence Linked Segments")
        
        ax[0].text(.10,.10,'n='+str(len(dikeset)), transform=ax[0].transAxes)
        ax[1].text(.10,.10,'n='+str(len(lines)), transform=ax[1].transAxes)
        ax[2].text(.10,.10,'n='+str(len(lines.loc[mask])), transform=ax[2].transAxes)
        plt.tight_layout()
        
        
        fig.savefig(ImgName, transparent=True)
        plt.close(fig)
    
def Toplines(lines, dikeset, path):
    ImgPath=path[:-4]+"Images/"
    ImgName=ImgPath+"TopLines"+".pdf"
    
    if MakeFigCheck(ImgName):
        rstep=lines['Rho_Threshold'].values[0]*10
        
        T, fig, ax=TopHTSection(lines, dikeset, rstep, 4, n=3)
        
        
        
        fig.savefig(ImgName, transparent=True)
        plt.close(fig)
    
def PairGrid(lines,cols, path,n=2):
    ImgPath=path[:-4]+"Images/"
    ImgName=ImgPath+"PairPlot"+".pdf"
    if MakeFigCheck(ImgName):
        mask=np.greater(lines['Size'].values, n)
        temp=lines.loc[mask][cols]
        
        
        g = sns.PairGrid(temp, diag_sharey=False, hue="Size", corner=True)
        g.map_diag(sns.histplot,  hue=None, color=".3")
        g.map_offdiag(sns.scatterplot)
        g.add_legend()
        
    
        g.savefig(ImgName, transparent=True)
    #plt.close(g)

def SaveFigCheck(fig, ImgPath, overide=False):
    
    if os.path.exists(ImgPath):
        mtime=time.ctime(os.path.getmtime(ImgPath))
        mtime=dateutil.parser.parse(mtime)
        now=datetime.now()
        
         
        tdif=mtime + timedelta(days=7) < now
    
    elif not os.path.exists(ImgPath):
        tdif=True
    
    if tdif:
        fig.savefig(ImgPath, transparent=True, dpi=600)
        print(ImgPath, "saved")
    elif not tdif:
        print(ImgPath, "already was saved < 2 days ago")
        
    if overide:
        fig.savefig(ImgPath, transparent=True, dpi=600)
        print(ImgPath, "saved by overide")

def LocPlots(path, lines, cols):
    
    for i in cols:
        ImgPath=path[:-4]+"Images/"
        ImgName=ImgPath+"ByLoc"+i+".pdf"
        if MakeFigCheck(ImgName):
            if i == "R_Length" or i =="R_Width":
                
                fig,ax=plotByAngle(lines, i, log_scale=(False,True))
            fig,ax=plotByAngle(lines, i)
            
            fig.savefig(ImgName, transparent=True, dpi=600)
    
def AnglePlots(path, lines, cols):
    print("Starting Angle Plots")
    for i in cols:
        ImgPath=path[:-4]+"Images/"
        ImgName=ImgPath+"ByAngle"+i+".pdf"
        print("making image:", ImgName)
        if MakeFigCheck(ImgName):
            if i == "R_Length" or i =="R_Width":
                
                fig=plotByAngle(lines, i, log_scale=True)
            fig=plotByAngle(lines, i)
            fig.savefig(ImgName, transparent=True, dpi=600)
    
        
    
    
def allFigures(path1,path2,w, Large=True, overide=False):
    
    dikeset=pd.read_csv(path1)
    lines=pd.read_csv(path2)
    lines['R_Width']=lines['R_Width']+w
    lines['R_Length']=lines['R_Length']/1000
    
    cols=['Size', 'Overlap', "MaxSegNNDist", "MinSegNNDist"]
    
    ImgPath=path2[:-4]+"Images/"
    sns.reset_orig()
    isExist = os.path.exists(ImgPath)
    if not isExist:
        os.makedirs(ImgPath)
        print("The new directory", ImgPath, " is created!")
    sns.reset_orig()
    m=lines['TrustFilter']==1
    
    g=sns.jointplot(data=lines, x='Size', y='R_Length', alpha=0.6, hue='TrustFilter',  ylim = (0,250))

    slope, intercept, r_value, p_value, std_err = stats.linregress(lines['Size'].values, lines['R_Length'].values)
    x=np.linspace(3, lines['Size'].max())
    y=slope*x+intercept 
    
    g.ax_joint.plot(x,y, 'b-.', alpha=0.8)
    g.ax_joint.annotate(f'$r = {r_value:.3f}, p = {p_value:.3f}$',
                    xy=(0.1, 0.9), xycoords='axes fraction',
                    ha='left', va='center',
                    bbox={'boxstyle': 'round', 'fc': 'powderblue', 'ec': 'navy'})
   
    
    slope, intercept, r_value, p_value, std_err = stats.linregress(lines['Size'].loc[m].values, lines['R_Length'].loc[m].values)
    x=np.linspace(4, lines['Size'].max())
    y=slope*x+intercept 
    
    g.ax_joint.plot(x,y,  '-.',color='orange', alpha=0.9)
    g.ax_joint.annotate(f'$r = {r_value:.3f}, p = {p_value:.3f}$',
                    xy=(0.1, 0.8), xycoords='axes fraction',
                    ha='left', va='center',
                    bbox={'boxstyle': 'round', 'fc': 'peachpuff', 'ec': 'coral'})


    #g.savefig(ImgPath+"SizebyLength2.pdf", transparent=True)
    SaveFigCheck(g, ImgPath+"SizebyLength2.pdf", overide=overide)
    SaveFigCheck(g, ImgPath+"SizebyLength2.png", overide=overide)
    
    
    #dilation 
    ImgName=ImgPath+"DoubleDilation"+".pdf"
    if MakeFigCheck(ImgName):
        fig,axes=DoubleDilationPlot(dikeset, lines, binWidth=1700, averageWidth=w)
        
        #fig.savefig(ImgName, transparent=True)
        SaveFigCheck(fig, ImgName, overide=overide)
        plt.close(fig)
        

    
    ImgName=ImgPath+"TripleDilationTrusted"+".pdf"
    if MakeFigCheck(ImgName):
        fig,axes=TripleDilationPlot(dikeset, lines, binWidth=1700, averageWidth=w)
        
        #fig.savefig(ImgName, transparent=True)
        SaveFigCheck(fig, ImgName, overide=overide)
        plt.close(fig)
        
        
    
    

    Histograms(dikeset, lines, path2)
    

    HT3(dikeset, lines, path2)
        
    Toplines(lines,dikeset, path2)
    PairGrid(lines, cols, path2)
    sns.reset_orig()
    rstep=lines['Rho_Threshold'].values[0]*10
    fig,ax=DotsLinesHist(lines, rstep, 4)
    fig.set_size_inches(15,5)
    
    plt.tight_layout()
    ImgName=ImgPath+"DotsLinesHist"+".pdf"
    ##fig.savefig(ImgName, transparent=True)
    SaveFigCheck(fig, ImgName, overide=overide)
    plt.close(fig)
    
    col1=['Overlap', 'Overlap', 'MaxSegNNDist']
    col2=['EnEchelonAngleDiff', 'MedianSegNNDist', 'EnEchelonAngleDiff' ]
    
    fig,ax=plotScatterHist(lines, "R_Length", "R_Width", log_scale=True)
    fig.set_size_inches(10,5)
    
    plt.tight_layout()
    ImgName=ImgPath+"ScatterHist"+"LengthWidth"+".pdf"
    #fig.savefig(ImgName, transparent=True)
    SaveFigCheck(fig, ImgName, overide=overide)
    plt.close(fig)
    
    for x,y in zip(col1, col2):
        
        if len(lines)>100:
            fig,ax=plotScatterHist(lines, x, y)
            fig.set_size_inches(10,5)
            
            plt.tight_layout()
            ImgName=ImgPath+"ScatterHist"+x+y+".pdf"
            #fig.savefig(ImgName, transparent=True)
            SaveFigCheck(fig, ImgName, overide=overide)
            plt.close(fig)
            
    for x,y in zip(col1, col2):
        if len(lines)>100:
            fig,ax=plotScatterHist(lines, x, y, TrustFilter=False)
            fig.set_size_inches(10,5)
            
            plt.tight_layout()
            ImgName=ImgPath+"ScatterHist_TrustFilter"+x+y+".pdf"
            #fig.savefig(ImgName, transparent=True)
            SaveFigCheck(fig, ImgName, overide=overide)
            plt.close(fig)
                
    colLoc=['AvgTheta', 'Overlap', 'Size', 'R_Length', 'R_Width', 'nOverlapingSegments', 'MedianSegNNDist']
    LocPlots(path2,lines,colLoc)
    
    colLoc=['Overlap', 'Size', 'R_Length', 'R_Width', 'nOverlapingSegments','MedianSegNNDist']
    AnglePlots(path2,lines,colLoc)
    
    if Large:
        plotLargeClusters(dikeset, lines, 10, path2, overide=False)
        
        
        
def reRunTest():
    allFigures("/home/akh/myprojects/Linking-and-Clustering-Dikes/dikedata/crb/CJDS_FebStraightened.csv", '/home/akh/myprojects/Linking-and-Clustering-Dikes/dikedata/crb/CJDS_Complete_2_433.csv', 8)
#def comparisonPlots():
    
    
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