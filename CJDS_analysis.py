#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 28 15:27:35 2022

@author: akh
"""
import pandas as pd 
import matplotlib.pyplot as plt 
import numpy as np 
from htMOD import AKH_HT as HT
from htMOD import MidtoPerpDistance
from clusterMod import HT_AGG_custom as AggHT
from examineMod import examineClusters
from plotmod import *

dikeset=pd.read_csv('/home/akh/myprojects/Linking-and-Clustering-Dikes/dikedata/crb/CJDS_FebStraightened.csv')
# already straightened dikes and did preprecessing

#do the hough transform 
theta,rho,xc,yc=HT(dikeset)
dikeset['theta']=theta
dikeset['rho']=rho
dikeset=MidtoPerpDistance(dikeset, xc, yc)

evals=pd.read_csv('/home/akh/myprojects/Linking-and-Clustering-Dikes/CJDS_sensitivity.csv')


# evals=pd.DataFrame()

# rhos=[200,500,1000,5000]
# #Run sensitivity analysis first
# for dtheta in [1,2,3,5,6]:
#     lengths=[]
#     needikes=[]
#     nclusters=[]
#     for drho in rhos:
#         print("Running ", dtheta, drho)
#         #Run the clustering algorithm
#         temp, clusters, M=AggHT(dikeset, dtheta, drho)
        
#         #Run examinecluster
#         temp,evaluation,EEdikes=examineClusters(temp, ifEE=True)
#         evaluation['Theta_Threshold']=dtheta
#         evaluation['Rho_Threshold']=drho
        
#         evals=evals.append(evaluation, ignore_index=True)

# evals.to_csv('/home/akh/myprojects/Linking-and-Clustering-Dikes/CJDS_sensitivity.csv')


# ax[0].plot(rhos, nclusters, "rs", label=name, markersize=7, alpha=0.7, markeredgecolor="r" )
# ax[1].plot(rhos, lengths, "bp", label=name, markersize=7, alpha=0.7, markeredgecolor="b")
# ax[2].plot(rhos, needikes, "g*", label=name, markersize=7, alpha=0.7, markeredgecolor="g")
dtheta=[1,2,3,5,6]
rhos=[200,500,1000,5000]
x=evals['Rho_Threshold'].to_numpy()
labels=evals['Theta_Threshold']
markers= ["rs-", "yv-", "g*-", "bp-", "kD-" ]
columns=evals.columns[1:20]

xlim=[(150,1050), (4950,5005)]

def plotSensitivityRhoColumns(evals, columns):
    xlim=[(150,1050), (4950,5005)]
    f1,ax11,ax12=breakXaxis(xlim, numAxes=len(columns))
    rhos=np.unique(evals['Rho_Threshold'])
    theta=np.unique(evals['Theta_Threshold'])
    markers= ["rs-", "yv-", "g*-", "bp-", "kD-" ]
    for i,j in zip(columns, range(len(columns))):
        
        for t,m in zip(theta, markers): 
            y1=evals.loc[evals['Theta_Threshold']==t][i]
            name=str(t)+"$^\circ$"
            kwargs=dict( label=name, markersize=7, alpha=0.5)
            plotBreak(xlim, rhos, y1, ax11[j], ax12[j], m, **kwargs)
            ax11[j].set_ylabel(i)
    plt.tight_layout()
    return f1, (ax11, ax12)


def plotSensitivityThetaColumns(evals, columns):
    
    f1,ax=plt.subplots(len(columns), sharex=True)
    rhos=np.unique(evals['Rho_Threshold'])
    theta=np.unique(evals['Theta_Threshold'])
    markers= ["rs-", "yv-", "g*-", "bp-" ]
    
    for i,j in zip(columns, range(len(columns))):
        
        for p,m in zip(rhos, markers): 
            y1=evals.loc[evals['Rho_Threshold']==p][i]
            name=str(p)+"m"
            kwargs=dict( label=name, markersize=7, alpha=0.5)
            ax[j].plot(theta, y1, m, **kwargs)
            ax[j].set_ylabel(i)
    
    ax[j].set_xlabel("Theta ($^\circ$) ")
    plt.tight_layout()
    
    return f1, ax


f1, ax1=plotSensitivityRhoColumns(evals, columns[1:5])
ax=ax1[0]
handles, labels = ax[0].get_legend_handles_labels()
f1.legend(handles, labels, loc='upper left')


f1.savefig("CJDS_sensitivty_rho1.eps",dpi=600)
f2, ax2=plotSensitivityRhoColumns(evals, columns[5:10])
f2.savefig("CJDS_sensitivty_rho2.eps",dpi=600)
f3, ax3=plotSensitivityRhoColumns(evals, columns[10:15])
f3.savefig("CJDS_sensitivty_rho3.eps",dpi=600)
f4, ax4=plotSensitivityRhoColumns(evals, columns[15:19])
f4.savefig("CJDS_sensitivty_rho4.eps",dpi=600)

f5, ax5=plotSensitivityThetaColumns(evals, columns[1:5])
handles, labels = ax5[0].get_legend_handles_labels()
f5.legend(handles, labels, loc='upper left')

f5.savefig("CJDS_sensitivty_rho5.eps",dpi=600)
f6, ax6=plotSensitivityThetaColumns(evals, columns[5:10])
f6.savefig("CJDS_sensitivty_rho6.eps",dpi=600)
f7, ax7=plotSensitivityThetaColumns(evals, columns[10:15])
f7.savefig("CJDS_sensitivty_rho7.eps",dpi=600)
f8, ax8=plotSensitivityThetaColumns(evals, columns[15:19])
f8.savefig("CJDS_sensitivty_rho8.eps",dpi=600)

# fig1, ax1=plt.subplots(10) 
# fig2, ax2=plt.subplots(9) 

# for i,j in zip(columns[0:10], range(len(columns[0:10]))):
    
#     for t,m in zip(dtheta, markers): 
#         y1=evals.loc[evals['Theta_Threshold']==t][i]
#         name=str(t)+"$^\circ$"
#         ax1[j].plot(rhos, y1, m, label=name, markersize=7, alpha=0.7)
#         ax1[j].set_title(i)
        
# for i,j in zip(columns[10:19], range(len(columns[10:19]))):
    
#     for t,m in zip(dtheta, markers): 
#         y1=evals.loc[evals['Theta_Threshold']==t][i]
#         name=str(t)+"$^\circ$"
#         ax2[j].plot(rhos, y1, m, label=name, markersize=7, alpha=0.7)
#         ax2[j].set_title(i)
