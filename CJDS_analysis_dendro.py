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
from clusterMod import HT_AGG_custom
from examineMod import examineClusters
from plotmod import *
import scipy.cluster.hierarchy as sch

dikeset=pd.read_csv('/home/akh/myprojects/Linking-and-Clustering-Dikes/dikedata/crb/CJDS_FebStraightened.csv')
# already straightened dikes and did preprecessing

#do the hough transform 
xc,yc=HT_center(dikeset)
# theta,rho,xc,yc=HT(dikeset)
# dikeset['theta']=theta
# dikeset['rho']=rho
# dikeset=MidtoPerpDistance(dikeset, xc, yc)

dtheta=3 
drho=500
        
#Run examinecluster

"""
are the scipy and sklearn  the same?
lines1,evaluation,EEdikes=examineClusters(dikeset, ifEE=True)
evaluation['Theta_Threshold']=dtheta
evaluation['Rho_Threshold']=drho

dikeset2, Z=HT_AGG_custom(dikeset, dtheta, drho)

lines2,evaluation2,EEdikes2=examineClusters(dikeset2, ifEE=True)
dikeset['TrueLabel']=dikeset["Labels"]
from examineMod import checkLineClusterChange
checkLineClusterChange(lines1, lines2)
yes
"""



# dikeset, Z=HT_AGG_custom2(dikeset, dtheta, drho)
# from examineMod import overlapSegments
# lines,evaluation,EEdikes=examineClusters(dikeset)


#lines=pd.read_csv('/home/akh/myprojects/Linking-and-Clustering-Dikes/dikedata/crb/CJDS_500m_3deg_level1.csv')



import scipy.cluster.hierarchy as sch

''' for loop to check 
# for i in np.arange(1,25,7):
    
#     labels=sch.fcluster(Z, t=i, criterion='distance')
    # dikeset5=dikeset
    # dikeset5['Labels']=labels
    # lines5,evaluation=examineClusters(dikeset5)
    # fig,ax=plt.subplots(); c=ax.scatter(lines5['Overlap'], lines5['EnEchelonAngleDiff'], c=lines5['AvgTheta'], alpha=0.6, edgecolor='k', cmap="turbo")
    # ax.set_xlabel("Segment Overlap (%)")
    # ax.set_title("Level"+str(i))
    # fig.colorbar(c, label='Average Angle ($^\circ$)')
    # ax.set_ylabel("En Echelon Angle Difference ($^\circ$)")
    
    # fig,ax=plt.subplots(1,4)
    # ax[0].scatter(b,p, edgecolor='k', alpha=0.7)
    # ax[0].set_xlabel('Birth')
    # ax[0].set_ylabel('Persistance')
    # ax[0].set_xlim((0,40))
    
    # ax[0].axvline(i,0, np.max(p), color="r")
    # ax[2].scatter(lines5['Size'], lines5['R_Length']/1000)
    # ax[2].set_xlabel('Length (km)')
    # ax[2].set_ylabel('Size')
    # k=str(evaluation['nClusters'].values)+" clustered dikes"
    # ax[1].text(0.1, 0.90, k, transform=ax[1].transAxes)
    # DotsHT(fig, ax[3], lines5, ColorBy='StdTheta')
    # plotlines(lines5, 'k', ax[1], alpha=0.7)
    # name="CJDS_AG_level"+str(i)+".png"
    # plt.tight_layout()
    # fig.savefig(name, dpi=600)
'''

Z=np.load('/home/akh/myprojects/Linking-and-Clustering-Dikes/dikedata/crb/CJDS_MarchLines_3_500_distances.npy')

#labels=sch.fcluster(Z, t=1, criterion='distance')
#dikeset5=dikeset
#dikeset['Labels']=labels
#lines,evaluation=examineClusters(dikeset)
lines= pd.read_csv("/home/akh/myprojects/Linking-and-Clustering-Dikes/dikedata/crb/CJDS_500m_3deg_level1.csv")

EEDikes_lines=lines.loc[lines['EnEchelonAngleDiff']>5]



fig,ax=plt.subplots(); c=ax.scatter(lines['Overlap'], lines['EnEchelonAngleDiff'], c=lines['AvgTheta'], alpha=0.6, edgecolor='k', cmap="turbo")
ax.set_xlabel("Segment Overlap (%)")
fig.colorbar(c, label='Average Angle ($^\circ$)')
ax.set_ylabel("En Echelon Angle Difference ($^\circ$)")


fig, ax, p, b=persistancePlot(Z)

fig,ax=plt.subplots(2,2)

DotsHT(fig, ax[0,0], lines, ColorBy='PerpOffsetDist')
ax[0,1].hist(lines['PerpOffsetDist'])
ax[0,1].set_xlabel("PerpOffsetDist (m)")
ax[0,1].set_ylabel("Counts")
psplit= np.mean(lines['PerpOffsetDist'].values)

lines1=lines[ lines['PerpOffsetDist'] > psplit]
lines2=lines[ lines['PerpOffsetDist'] < psplit]

DotsHT(fig, ax[1,0], lines1, ColorBy='PerpOffsetDist')
DotsHT(fig, ax[1,1], lines2, ColorBy='PerpOffsetDist')


x=np.linspace(-90,90,10)
slope1=(-2.8e4-2.9e4)/(-87.9-87.9)
b=2.9e4-87.9*slope1
ax[1,0].plot([-87.9, 87.9],[-2.8e4, 2.9e4], 'r*-', linewidth=2)

def mysigmoid(x, c=1, a=1, b=0):
    return a/(1+np.e**( -c*(x-c)))+b

ax[1,1].plot(x, mysigmoid(x, 0.06, -100000, 50000), 'y', linewidth=2)
x=np.random.normal(0, 30, 200) #np.linspace(-90,90,100)
y1=x*slope1+b
y2=mysigmoid(x, 0.06, -100000, 50000)

from synthetic import fromHT, makeRadialSwarmdf

df1=fromHT(x,y1, length=50000, xc=xc, yc=yc, CartRange=300000)
df2=fromHT(x,y2, length=50000, xc=xc, yc=yc, CartRange=300000)
label=[1]*len(df1)+[2]*len(df2)

linefit=df1.append(df2)
linefit['Label']=label

fig2,ax2=DotsLines(linefit, ColorBy='Label', cmap='viridis')
ax2[1].scatter(lines['AvgTheta'], lines['AvgRho'], alpha=0.1)
radEx=makeRadialSwarmdf(50000, center=[xc,yc])

#ig2,ax2=DotsLines(radEx, ColorBy='theta', cmap='viridis')



# labels15=sch.fcluster(Z, t=15, criterion='distance')
# dikeset15=dikeset
# dikeset15['Labels']=labels15
# lines15,evaluation=examineClusters(dikeset15)

# labels25=sch.fcluster(Z, t=25, criterion='distance')
# dikeset25=dikeset
# dikeset25['Labels']=labels25
# lines25,evaluation=examineClusters(dikeset25)

# labels35=sch.fcluster(Z, t=35, criterion='distance')
# dikeset35=dikeset
# dikeset35['Labels']=labels35
# lines35,evaluation=examineClusters(dikeset35)


# fig,ax=plt.subplots()
# ax.scatter(lines['Overlap'], lines['EnEchelonAngleDiff'], c="red", alpha=0.6, edgecolor='k', label="1")
# ax.scatter(lines5['Overlap'], lines5['EnEchelonAngleDiff'],  c="yellow", alpha=0.6, edgecolor='k', label="5")
# ax.scatter(lines15['Overlap'], lines15['EnEchelonAngleDiff'], c="green", alpha=0.6, edgecolor='k', label="15")
# ax.scatter(lines25['Overlap'], lines25['EnEchelonAngleDiff'], c="blue", alpha=0.6, edgecolor='k', label="25")
# ax.scatter(lines35['Overlap'], lines35['EnEchelonAngleDiff'], c="purple", alpha=0.6, edgecolor='k', label="35")
# ax.set_xlabel("Segment Overlap (%)")
# #fig.colorbar(c, label='Average Angle ($^\circ$)')
# ax.set_ylabel("En Echelon Angle Difference ($^\circ$)")
# ax.legend()

# fig,ax=plt.subplots()
# plotlines(lines, 'r', ax, alpha=0.6)
# plotlines(lines15, 'green', ax)
# plotlines(lines35, 'purple', ax)

# evals=pd.read_csv('/home/akh/myprojects/Linking-and-Clustering-Dikes/CJDS_sensitivity.csv')


# # evals=pd.DataFrame()

# # rhos=[200,500,1000,5000]
# # #Run sensitivity analysis first
# # for dtheta in [1,2,3,5,6]:
# #     lengths=[]
# #     needikes=[]
# #     nclusters=[]
# #     for drho in rhos:
# #         print("Running ", dtheta, drho)
# #         #Run the clustering algorithm
# #         temp, clusters, M=AggHT(dikeset, dtheta, drho)
        
# #         #Run examinecluster
# #         temp,evaluation,EEdikes=examineClusters(temp, ifEE=True)
# #         evaluation['Theta_Threshold']=dtheta
# #         evaluation['Rho_Threshold']=drho
        
# #         evals=evals.append(evaluation, ignore_index=True)

# # evals.to_csv('/home/akh/myprojects/Linking-and-Clustering-Dikes/CJDS_sensitivity.csv')


# # ax[0].plot(rhos, nclusters, "rs", label=name, markersize=7, alpha=0.7, markeredgecolor="r" )
# # ax[1].plot(rhos, lengths, "bp", label=name, markersize=7, alpha=0.7, markeredgecolor="b")
# # ax[2].plot(rhos, needikes, "g*", label=name, markersize=7, alpha=0.7, markeredgecolor="g")
# dtheta=[1,2,3,5,6]
# rhos=[200,500,1000,5000]
# x=evals['Rho_Threshold'].to_numpy()
# labels=evals['Theta_Threshold']
# markers= ["rs-", "yv-", "g*-", "bp-", "kD-" ]
# columns=evals.columns[1:20]

# xlim=[(150,1050), (4950,5005)]

# def plotSensitivityRhoColumns(evals, columns):
#     xlim=[(150,1050), (4950,5005)]
#     f1,ax11,ax12=breakXaxis(xlim, numAxes=len(columns))
#     rhos=np.unique(evals['Rho_Threshold'])
#     theta=np.unique(evals['Theta_Threshold'])
#     markers= ["rs-", "yv-", "g*-", "bp-", "kD-" ]
#     for i,j in zip(columns, range(len(columns))):
        
#         for t,m in zip(theta, markers): 
#             y1=evals.loc[evals['Theta_Threshold']==t][i]
#             name=str(t)+"$^\circ$"
#             kwargs=dict( label=name, markersize=7, alpha=0.5)
#             plotBreak(xlim, rhos, y1, ax11[j], ax12[j], m, **kwargs)
#             ax11[j].set_ylabel(i)
#     plt.tight_layout()
#     return f1, (ax11, ax12)


# def plotSensitivityThetaColumns(evals, columns):
    
#     f1,ax=plt.subplots(len(columns), sharex=True)
#     rhos=np.unique(evals['Rho_Threshold'])
#     theta=np.unique(evals['Theta_Threshold'])
#     markers= ["rs-", "yv-", "g*-", "bp-" ]
    
#     for i,j in zip(columns, range(len(columns))):
        
#         for p,m in zip(rhos, markers): 
#             y1=evals.loc[evals['Rho_Threshold']==p][i]
#             name=str(p)+"m"
#             kwargs=dict( label=name, markersize=7, alpha=0.5)
#             ax[j].plot(theta, y1, m, **kwargs)
#             ax[j].set_ylabel(i)
    
#     ax[j].set_xlabel("Theta ($^\circ$) ")
#     plt.tight_layout()
    
#     return f1, ax


# f1, ax1=plotSensitivityRhoColumns(evals, columns[1:5])
# ax=ax1[0]
# handles, labels = ax[0].get_legend_handles_labels()
# f1.legend(handles, labels, loc='upper left')


# f1.savefig("CJDS_sensitivty_rho1.eps",dpi=600)
# f2, ax2=plotSensitivityRhoColumns(evals, columns[5:10])
# f2.savefig("CJDS_sensitivty_rho2.eps",dpi=600)
# f3, ax3=plotSensitivityRhoColumns(evals, columns[10:15])
# f3.savefig("CJDS_sensitivty_rho3.eps",dpi=600)
# f4, ax4=plotSensitivityRhoColumns(evals, columns[15:19])
# f4.savefig("CJDS_sensitivty_rho4.eps",dpi=600)

# f5, ax5=plotSensitivityThetaColumns(evals, columns[1:5])
# handles, labels = ax5[0].get_legend_handles_labels()
# f5.legend(handles, labels, loc='upper left')

# f5.savefig("CJDS_sensitivty_rho5.eps",dpi=600)
# f6, ax6=plotSensitivityThetaColumns(evals, columns[5:10])
# f6.savefig("CJDS_sensitivty_rho6.eps",dpi=600)
# f7, ax7=plotSensitivityThetaColumns(evals, columns[10:15])
# f7.savefig("CJDS_sensitivty_rho7.eps",dpi=600)
# f8, ax8=plotSensitivityThetaColumns(evals, columns[15:19])
# f8.savefig("CJDS_sensitivty_rho8.eps",dpi=600)

# # fig1, ax1=plt.subplots(10) 
# # fig2, ax2=plt.subplots(9) 

# # for i,j in zip(columns[0:10], range(len(columns[0:10]))):
    
# #     for t,m in zip(dtheta, markers): 
# #         y1=evals.loc[evals['Theta_Threshold']==t][i]
# #         name=str(t)+"$^\circ$"
# #         ax1[j].plot(rhos, y1, m, label=name, markersize=7, alpha=0.7)
# #         ax1[j].set_title(i)
        
# # for i,j in zip(columns[10:19], range(len(columns[10:19]))):
    
# #     for t,m in zip(dtheta, markers): 
# #         y1=evals.loc[evals['Theta_Threshold']==t][i]
# #         name=str(t)+"$^\circ$"
# #         ax2[j].plot(rhos, y1, m, label=name, markersize=7, alpha=0.7)
# #         ax2[j].set_title(i)
