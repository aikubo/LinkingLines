#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 21 14:07:06 2021

@author: akh
"""
import pandas as pd
from htMOD import AKH_HT as HT
from sklearn.preprocessing import scale
import numpy as np 
from clusterMod import *
import matplotlib.pyplot as plt
from plotmod import plotlines, labelcolors, plotbyAngle, BA_HT, HThist
from examineMod import examineClusters
import seaborn as sns
from jitteringHTcenter import moveHTcenter, rotateHT
from matplotlib import cm

# def fragmentDikes(df, maxL=20000, ndikesMax=None, distortion=0):
#     np.random.seed(5)
#     dfFrag=pd.DataFrame(columns=df.columns)
#     ndikes=0
#     if ndikesMax is None: 
#         ndikesMax=2000
        
       
#     for i in range(len(df)):
#         nSegments=np.random.randint(3)
#         if ndikes < ndikesMax: 
#             for j in range(nSegments):
#                 high=max(0, df['Xend'].iloc[i])
#                 low=min(0, df['Xend'].iloc[i])
#                 xrange=np.random.randint(low,high, size=2)
#                 m=df['Slope'].iloc[i]*(1+np.random.rand()*distortion)
#                 yrange=m*xrange
                
#                 L=np.sqrt((xrange[0]-xrange[1])**2+(yrange[0]-yrange[1])**2)
#                 if L > maxL: 
#                     continue
                
#                 if max(abs(yrange)) > df['Yend'].max() :
#                     continue 
#                 dfFrag=dfFrag.append(pd.DataFrame({'Xstart':xrange[0], 'Xend': xrange[1], 
#                                                    'Ystart': yrange[0], 'Yend':yrange[1], 
#                                                    'Slope':m, 'Length':L
#                                                        }, index=[0]), ignore_index=True)
#                 ndikes=ndikes+1
                
#     # a=np.full(len(dfFrag), False)
#     # a[:int(len(dfFrag)/mask)]=True 
#     # np.random.shuffle(a)
#     # dfFrag=dfFrag.iloc[a]
    
    
#     return dfFrag 


def makeRadialSwarm(radius, doubled=True, anglestart=-90, anglestop=90, ndikes=50, center=[0,0]):
    
    #center=np.array([0,0])
    angles=np.linspace(anglestart, anglestop, ndikes)
    m=np.tan(np.deg2rad(angles))
    
    Xstart=np.zeros(ndikes)+center[0]
    Ystart=np.zeros(ndikes)+center[1]
    Xend=radius/np.sqrt(1+m**2)+center[0]
    Yend=Xend*m+center[1]
    
    if doubled:
        
        Xstart=np.append(Xstart, np.zeros(ndikes)+center[0])
        Ystart=np.append(Ystart, np.zeros(ndikes)+center[1])
        Yend=np.append(Yend, -1*Xend*m+center[1])
        Xend=np.append(Xend, -1*radius/np.sqrt(1+m**2)+center[0])
        
        
    Xstart[abs(Xstart)<10**-6]=1
    print(center)
    return Xstart, Ystart, Xend, Yend

def makeRadialSwarmdf(radius, doubled=True, anglestart=-90, anglestop=90, ndikes=50, center=[0,0]):
    
    #center=np.array([0,0])
    angles=np.linspace(anglestart, anglestop, ndikes)
    m=np.tan(np.deg2rad(angles))
    
    Xstart=np.zeros(ndikes)+center[0]
    Ystart=np.zeros(ndikes)+center[1]
    Xend=radius/np.sqrt(1+m**2)+center[0]
    Yend=Xend*m+center[1]
    
    if doubled:
        
        Xstart=np.append(Xstart, np.zeros(ndikes)+center[0])
        Ystart=np.append(Ystart, np.zeros(ndikes)+center[1])
        Yend=np.append(Yend, -1*Xend*m+center[1])
        Xend=np.append(Xend, -1*radius/np.sqrt(1+m**2)+center[0])
        
        
    Xstart[abs(Xstart)<10**-6]=1
    
    df=pd.DataFrame({'Xstart':Xstart, 'Xend': Xend, 'Ystart': Ystart, 'Yend':Yend})
                     
    return df

def addSwarms(dflist):
    dfSwarm=pd.DataFrame()
    
    for i, b in zip(dflist, range(len(dflist))):
        i['Label']=[b]*len(i)
        dfSwarm=dfSwarm.append(i)

    return dfSwarm

def makeLinearSwarm(length, slope,  ndikes=20):
    
    Xstart=np.linspace(-1*length,length, ndikes)
    b=np.random.rand(ndikes)*2*length
    Ystart=Xstart*slope+b
    Xend=Xstart+length/np.sqrt(1+slope**2)
    Yend=Xend*slope+b 
    
    
    return Xstart, Ystart, Xend, Yend

# make linear dike swarms of angle and rho distributions 
def makeLinear2(length, angle, angleSTD, rho, rhoSTD, ndikes=100, CartRange=300000):
    angles=np.random.normal(angle, angleSTD, ndikes)
    rhos=np.random.normal(rho, rhoSTD, ndikes)

    b=rhos/np.sin(np.deg2rad(angles))
    slopes=-1/(np.tan(np.deg2rad(angles))+0.000000001)
    Xstart=np.random.normal(rho, rhoSTD, ndikes)
    
    Ystart=slopes*Xstart+b
    Xend=Xstart-length/np.sqrt(1+slopes**2)
    Yend=slopes*Xend+b

    df=pd.DataFrame({'Xstart':Xstart, 'Xend': Xend, 'Ystart': Ystart, 'Yend':Yend, 'theta':angles, 'rho':rhos})
    
    df=df.drop(df[ abs(df['Ystart']) > CartRange].index)

    return df


# def makeLinearSwarmdf(length, angle, rho=0, ndikes=20, angleSTD=1, rhoSTD=10000, CartRange=300000):
    
    
#     angles=np.random.normal(angle, angleSTD, ndikes)
#     slopes=-1/np.tan(np.deg2rad(angles))
#     Xstart= np.random.random_sample(0, rhoSTD*5,ndikes)#(3*length)*np.random.random_sample(ndikes) + 2*length#np.linspace(-1*length,length, ndikes)  
#     b=np.random.normal(rho, rhoSTD, ndikes)*np.sin(np.deg2rad(angles))
    
#     Ystart=Xstart*slopes+b
#     Xend=Xstart-length**2/(1+slopes**2)
#     Yend=Xend*slopes+b 
    
#     df=pd.DataFrame({'Xstart':ystart, 'Xend':Yend 'Ystart': Xstart, 'Yend':Yend})
    
#     largenss= (abs(Xstart-Xend) > CartRange) | (abs(Ystart-Yend) > CartRange)
#     #df.drop(largenss)
                     
#     return df


def OverLappingSwarms(nlinear, nradial, A1,A2, slope, center):
    length=45000
    
    X1, Y1, X2, Y2=makeRadialSwarm(length, doubled=False, anglestart=A1, anglestop=A2, ndikes=nradial, center=center)
    X3, Y3, X4, Y4=makeLinearSwarm(length,slope,ndikes=nlinear)
    Xstart=np.append(X1,X3)
    Ystart=np.append(Y1,Y3)
    
    Xend=np.append(X2,X4)
    Yend=np.append(Y2,Y4)
    
    #print(len(Xstart), len(Ystart))
    
    iddike=np.arange(0, len(Xstart))
    dlength=np.sqrt( (Xstart-Xend)**2 + (Ystart-Yend)**2)
    reg=np.ones(len(X1))
    
    reg=np.append(reg, np.zeros(nlinear))
    """ Regime label, 1 for radial 0 for linear"""
    #print(len(iddike), len(reg), len(dlength))
    df=pd.DataFrame({'Xstart':Xstart, 'Xend': Xend, 'Ystart': Ystart, 'Yend':Yend, 'ID': iddike, 'seg_length':dlength, 'Regime': reg})


    return df

def ManyRadial(nradial):
    length=5000
    ndikes=50
    
    Xstart, Ystart, Xend, Yend=makeRadialSwarm(length, doubled=True, anglestart=30, anglestop=90, ndikes=ndikes, center=[0,0])
    
    for i in range(nradial-1):
        center= [5000*(i+1),5000*(i+1)]
        X1=Xstart+center[0]
        Y1=Ystart+center[1]
        X2=Xend+center[0]
        Y2=Yend+center[1]
        print(center)
        Xstart=np.append(Xstart,X1)
        Ystart=np.append(Ystart,Y1)
    
        Xend=np.append(Xend,X2)
        Yend=np.append(Yend,Y2)

    
    #print(len(Xstart), len(Ystart))
    
    iddike=np.arange(0, len(Xstart))
    dlength=np.sqrt( (Xstart-Xend)**2 + (Ystart-Yend)**2)

    """ Regime label, 1 for radial 0 for linear"""
    #print(len(iddike), len(reg), len(dlength))
    df=pd.DataFrame({'Xstart':Xstart, 'Xend': Xend, 'Ystart': Ystart, 'Yend':Yend, 'ID': iddike, 'seg_length':dlength})


    return df

    
def RandomDikes(n, CartRange=100000):

    Xstart=np.random.uniform(-1*CartRange,CartRange,n)
    Ystart=np.random.uniform(-1*CartRange,CartRange,n)
    Xend=np.random.uniform(-1*CartRange,CartRange,n)
    Yend=np.random.uniform(-1*CartRange,CartRange,n)

    angles=np.arctan2(Yend-Ystart, Xend-Xstart)
   


    df=pd.DataFrame({'Xstart':Xstart, 'Xend': Xend, 'Ystart': Ystart, 'Yend':Yend})
    
    df=df.drop(df[ abs(df['Ystart']) > CartRange].index)
    return df

def standardSpacing(ndikes, angle, RhoStart, RhoSpacing, CartRange=100000 ):
    angles=np.ones(ndikes)*angle #np.random.normal(angle, angleSTD, ndikes)
    RhoEnd=RhoStart+ndikes*RhoSpacing
    rhos= np.arange(RhoStart,RhoEnd, RhoSpacing)
    length=3*RhoSpacing
    b=rhos/np.sin(np.deg2rad(angles))
    slopes=-1/(np.tan(np.deg2rad(angles))+0.000000001)
    Xstart=RhoStart #np.random.normal(rho, rhoSTD, ndikes)
    
    Ystart=slopes*Xstart+b
    Xend=Xstart-length/np.sqrt(1+slopes**2)
    Yend=slopes*Xend+b

    df=pd.DataFrame({'Xstart':Xstart, 'Xend': Xend, 'Ystart': Ystart, 'Yend':Yend, 'theta':angles, 'rho':rhos, 'Label': np.ones(ndikes)})
    
    df=df.drop(df[ abs(df['Ystart']) > CartRange].index)
    
    return df 

def EnEchelonSynthetic(ndikes, angle, RhoStart, RhoSpacing, Overlap=0, CartRange=100000 ):
    angles=np.ones(ndikes)*angle #np.random.normal(angle, angleSTD, ndikes)
    RhoEnd=RhoStart+ndikes*RhoSpacing
    rhos= np.arange(RhoStart,RhoEnd, RhoSpacing)
    length=3*RhoSpacing
    b=rhos/np.sin(np.deg2rad(angles))
    slopes=-1/(np.tan(np.deg2rad(angles))+0.000000001)
    #
    Xstart=np.linspace(RhoStart,RhoEnd,ndikes) #np.random.normal(rho, rhoSTD, ndikes)
    
    Ystart=slopes*Xstart+b
    Xend=Xstart-length/np.sqrt(1+slopes**2)
    Yend=slopes*Xend+b
    l=np.sqrt( (Xstart-Xend)**2 + (Ystart-Yend)**2)
    df=pd.DataFrame({'Xstart':Xstart, 'Xend': Xend, 'Ystart': Ystart, 'Yend':Yend, 'theta':angles, 'seg_length':l, 'rho':rhos, 'Label': np.ones(ndikes)})
    
    df=df.drop(df[ abs(df['Ystart']) > CartRange].index)
    
    return df 

def fromHT(angles, rhos, length=10000, xc=0, yc=0, CartRange=100000, label=None, xrange=None):
    ndikes=len(angles)
    b=rhos/np.sin(np.deg2rad(angles))
    slopes=-1/(np.tan(np.deg2rad(angles))+0.000000001)
    
    if xrange is None:
        Xstart=np.random.normal(-np.max(abs(rhos)), np.max(abs(rhos)), ndikes)
        Xend=Xstart-length/np.sqrt(1+slopes**2)
    else:
        Xstart=xc-xrange/2
        Xend=xc+xrange/2
    
    
    Ystart=slopes*Xstart+b
    
    Yend=slopes*Xend+b
    
    l=np.sqrt( (Xstart-Xend)**2 + (Ystart-Yend)**2)
    df=pd.DataFrame({'Xstart':Xstart, 'Xend': Xend, 'Ystart': Ystart, 'Yend':Yend, 'theta':angles, 'seg_length':l, 'rho':rhos, 'Label': np.ones(ndikes)})
    
    df=df.drop(df[ abs(df['Ystart']) > CartRange].index)
    
    
    return df 

def fragmentDikes(df):
    m=(df['Ystart'].values-df['Yend'].values)/(df['Xstart'].values-df['Xend'].values)
    df['Slope']=m
    bint=df['Ystart'].values-m*df['Xstart'].values
    dfFrag=pd.DataFrame()
    print(len(m), len(bint))
    for i in range(len(df)):
        nSegments=5 #np.random.randint(2,7)
        high=max(df['Xend'].iloc[i], df['Xstart'].iloc[i])
        low=min(df['Xend'].iloc[i], df['Xstart'].iloc[i])
        xrange1=np.random.randint(low,high, size=nSegments)
        xrange2=np.random.randint(low,high, size=nSegments)
        m2=(m[i]+np.random.rand()*np.deg2rad(10))

        yrange1=m2*xrange1+bint[i]
        yrange2=m2*xrange2+bint[i]
        L=np.sqrt((xrange1-xrange2)**2+(yrange1-yrange2)**2)

        dfFrag=dfFrag.append(pd.DataFrame({'Xstart':xrange1, 'Xend': xrange2, 
                                           'Ystart': yrange1, 'Yend':yrange2,
                                           'seg_length':L, 'Xmid': (xrange1+xrange2)/2, 'Ymid': (yrange1+yrange2)/2 ,
                                           'Label': np.ones(nSegments)*i, 
                                           'Original Label':np.ones(nSegments)*df['Label'].iloc[i]
                                               }), ignore_index=True)
    return dfFrag


# CRB 500-300km
# min(lines['Ystart'])-max(lines['Ystart'])
# Out[7]: -530604.0

# min(lines['Xstart'])-max(lines['Xstart'])
# Out[8]: -325692.0

# radius=10000
# df1=makeRadialSwarmdf(radius, doubled=False, anglestart=-90, anglestop=90, ndikes=200, center=[0,0])
# df2=makeRadialSwarmdf(radius, doubled=False, anglestart=-40, anglestop=40, ndikes=200, center=[40000,40000])
# df3=makeRadialSwarmdf(radius, doubled=False, anglestart=-20, anglestop=50, ndikes=200, center=[-100000,100000])
# df4=makeRadialSwarmdf(radius, doubled=False, anglestart=-10, anglestop=20, ndikes=200, center=[200000,-300000])
# df5=makeRadialSwarmdf(radius, doubled=False, anglestart=-10, anglestop=15, ndikes=200, center=[10000,500000])
# #df6=makeLinearSwarmdf(radius, 2, rho=500000, ndikes=1000, angleSTD=10, rhoSTD=20000)
# df6, theta1, rho1=makeLinear2(radius, 6,30, 3200, 5000, ndikes=1000)
# df=addSwarms([df1,df2,df3,df4,df5,df6])

# lines=df6
# fig,ax=plt.subplots(1,2)
# plotlines(df6, 'k', ax[0], center=True)
# theta, rho, xc, yc= HT(lines, yc=0, xc=0)
# lines['AvgRho']=rho
# lines['AvgTheta']=theta
# print("mean",np.mean(theta), "std", np.std(theta) )
# print("mean",np.mean(rho), "std", np.std(rho) )
# print( np.sum(np.int64(theta)==np.int64(theta1)), np.sum(np.int64(rho)==np.int64(rho1)))
# ax[1].scatter(theta,rho, c=lines['label'], cmap=cm.turbo)
# ax[0].axis('equal')
# ax[1].set_xlim([-90,90])
# df.to_csv('/home/akh/myprojects/Linking-and-Clustering-Dikes/syntheticRadial_testfindRadial_doubledFalse.csv')
#df=OverLappingSwarms(500,500, 20, 40, 0.2, center=[0,0])

#df2=OverLappingSwarms(500,500, -60, 40, 4, center=[40000,20000])


# df=df.append(df2,  ignore_index=True)



# #dfFrag1=fragmentDikes(df, distortion=1)


# #dfFrag1.to_csv('/home/akh/myprojects/Linking-and-Clustering-Dikes/syntheticFragmented.csv')
# #df=rotateData2(df, 50)

# dfFrag=fragmentDikes(df)

# plotlines(dfFrag, 'r', ax[0], center=True)
# theta1,rho1, xc, yc=AKH_HT(df.astype(float), xc=0, yc=0)
# df['theta']=theta1
# df['rho']=rho1
# df['Labels']=df['ID']
# #lines=examineClusters(df)
# dfFrag.to_csv('/home/akh/myprojects/Linking-and-Clustering-Dikes/syntheticFrag_large.csv', index=False)
# df.to_csv('/home/akh/myprojects/Linking-and-Clustering-Dikes/syntheticTrue_large.csv', index=False)

# rdf=OverLappingSwarms(500,500, 20, 40, 0.2, center=[0,0]) #ManyRadial(2)
# rdfFrag=fragmentDikes(rdf)
# #rdfFrag.to_csv('/home/akh/myprojects/Linking-and-Clustering-Dikes/syntheticManyRadial.csv', index=False)
# #rdf.to_csv('/home/akh/myprojects/Linking-and-Clustering-Dikes/syntheticTrue_ManyRadial.csv', index=False)
# fig, ax =plt.subplots(1,3)
# theta1,rho1, xc, yc=AKH_HT(rdf.astype(float), xc=0, yc=0)
# rdf['theta']=theta1
# rdf['rho']=rho1

# plotlines(rdf, 'k', ax[0], center=True)
# theta, rho, xc, yc=AKH_HT(rdf, xc=0, yc=0)
# xr1, yr1=[0, 0]
# xr2, yr2=[5000, 5000]
# A1=np.sqrt( (xr1-xc)**2 + (yr1-yc)**2)
# A2=np.sqrt( (xr2-xc)**2 + (yr2-yc)**2)
# phi1=0#np.arctan( (yr1-yc)/(xr1-xc) )
# phi2=np.arctan( (yr2-yc)/(xr2-xc) )
# rhoRadial1=(xr1-xc)*np.cos(np.deg2rad(theta))+(yr1-yc)*np.sin(np.deg2rad(theta))#  A1*np.sin(theta-np.rad2deg(phi1))#(xr1-xc)*np.cos(theta)+(yr1-yc)*np.sin(theta)
# rhoRadial2=(xr2-xc)*np.cos(np.deg2rad(theta))+(yr2-yc)*np.sin(np.deg2rad(theta))
# #A2*np.sin(theta-np.rad2deg(phi2))#

# c1=ax[1].scatter( theta, rho, c=(rho-rhoRadial1)**2, cmap=cm.Reds, edgecolor="black")
# ax[1].plot( theta, rhoRadial1, "k")
# cbar=fig.colorbar(c1, ax=ax[1])

# c2=ax[2].scatter( theta, rho, c=(rho-rhoRadial2)**2, cmap=cm.Reds,edgecolor="black")
# ax[2].plot( theta, rhoRadial2, "k")
# cbar=fig.colorbar(c2, ax=ax[2])




# ax[1][0].scatter(theta1, rho1, edgecolor='black')
# ax[0][0].plot(0,0, "*", mec='black', markersize=20)
# ax[1][1].set_xlabel('Theta (deg)')
# ax[1][0].set_xlabel('Theta (deg)')

# ax[0][1].set_label('Perfect Radial Swarm')
# ax[0][0].set_label('With Noise')

# ax[1][0].set_ylabel('Rho (m)')

#dfFrag1=fragmentDikes(df, distortion=1)
d=20000

#rotateHT(df, 20)

#l, evaluation=examineClusters(df)

# plotlines(dfFrag1, 'r', ax[0][1], linewidth=2)
# plotlines(df, 'k', ax[0][1], alpha=0.1)
# theta2, rho2, xc2, yc2=AKH_HT(dfFrag1.astype(float), xc=0, yc=0)
#ax[1].scatter(l['AvgTheta'], l['AvgRho'], edgecolor='black', c=l['KNN2'])
# ax[0][1].plot(0,0, "*", mec='black', markersize=20)

""" make gif of center change """
# fig,ax=plt.subplots(1,2)
# fig.set_size_inches(12,6)
# rold=rho1
# told=theta1
# n=0

# plotlines(df, 'grey', ax[0], center=True)
# plotlines(dfFrag1, 'r', ax[0], linewidth=2)

# ax[0].plot(0,0, "*", mec='black', markersize=20)
# ax[1].scatter(theta1,rho1, edgecolor='black')
# ax[1].set_ylim([-30000, 30000,])
# ax[1].set_xlabel('Theta (deg)')
# ax[1].set_ylabel('Rho (m)')

#title="Center at ["+str(0) +"m ,"+str(0)+" m]"
#ax[0].set_title(title)
# n=0
# plt.tight_layout()
# name="radial"+str(n)+".png"
# fig.savefig(name, dpi=600)
# plt.tight_layout()
# xcs=np.array([0,d,.5*d, -0.2*d])
# ycs=np.array([d,d,.7*d, -d])

# fig,ax=plt.subplots(1,2)
# fig.set_size_inches(12,6)

# plotlines(df, 'grey', ax[0], center=True)



#plotlines(dfFrag1, 'r', ax[0], linewidth=2)
# for ic in range(len(xcs)):

#     i=xcs[ic]
#     j=ycs[ic]

#     #ax[1].scatter(told,rold, c='grey', s=20)
    

#     theta2, rho2, xc2, yc2=AKH_HT(df.astype(float), xc=i, yc=j)
#     ax[0].plot(i,j, "*", mec='black', markersize=20)
#     ax[1].scatter(theta2,rho2, edgecolor='black')
#     ax[1].set_ylim([-30000, 30000,])
#     ax[1].set_xlabel('Theta (deg)')
#     ax[1].set_ylabel('Rho (m)')
#     dist=np.sqrt(i**2+j**2)
#     print(dist, max(abs(rho2)))
#     #title="Center at ["+str(i) +"m ,"+str(j)+" m]|"
#     #ax[0].set_title(title)
#     rold=np.append(rold, rho2)
#     told=np.append(told, theta2)
#     n=n+1
#     plt.tight_layout()
#     #name="radial"+str(n)+".png"

""" """

#fig.savefig(name, dpi=600)
    
# dfFrag2=fragmentDikes(df, ndikesMax=40, distortion=0)
# dfFrag3=fragmentDikes(df, ndikesMax=20, distortion=5)


# 
# theta2, rho2, xc2, yc2=AKH_HT(dfFrag1.astype(float), xc=20000, yc=20000)
# theta3, rho3, xc3, yc3=AKH_HT(dfFrag1.astype(float), xc=-20000, yc=20000)
# theta4, rho4, xc4, yc4=AKH_HT(dfFrag1.astype(float), xc=20000, yc=-20000)
# theta5, rho5, xc5, yc5=AKH_HT(dfFrag1.astype(float), xc=-20000, yc=-20000)

# plotlines(dfFrag1, 'r', ax[0][1], linewidth=3, center=True, xc=xc2, yc=yc2)

# plotlines(df, 'grey', ax[0][1], alpha=0.1)

# plotlines(dfFrag2, 'r', ax[0][2], linewidth=3)
# plotlines(df, 'grey', ax[0][2], alpha=0.1, center=True, xc=xc3, yc=yc3)
# plotlines(dfFrag3, 'r', ax[0][3], linewidth=3)
# plotlines(df, 'grey', ax[0][3], alpha=0.1, center=True, xc=xc4, yc=yc4)
                    
# 
# ax[1][1].scatter(theta2, rho2)
# ax[1][2].scatter(theta3, rho3)
# ax[1][3].scatter(theta4, rho4)
# ax[1][0].set_ylabel('Rho (m) ')

# for i in range(4):
#     ax[1][i].set_xlabel('Theta (deg)')

# ax[0][0].set_title('Centered')
# ax[0][1].set_title('[20000,20000]')
# ax[0][2].set_title('[-20000,20000]')
# ax[0][3].set_title('[20000,-20000]')
# ax[0][4].set_title('[-20000,-20000]')

# fig, ax =plt.subplots(2,4)
# # 1 complete info
# # 2 Short dikes only 
# # 3 Little Dikes only 
# # 4 Short few dikes 
# plotlines(df, 'k', ax[0][0])

# dfFrag1=fragmentDikes(df, maxL=5000)
# dfFrag2=fragmentDikes(df, maxL=10000, ndikesMax=20)
# dfFrag3=fragmentDikes(df, maxL=1000, ndikesMax=20)


# theta1,rho1, xc, yc=AKH_HT(df.astype(float))
# theta2, rho2, xc, yc=AKH_HT(dfFrag1.astype(float))
# theta3, rho3, xc, yc=AKH_HT(dfFrag2.astype(float))
# theta4, rho4, xc, yc=AKH_HT(dfFrag3.astype(float))


# plotlines(dfFrag1, 'r', ax[0][1], linewidth=3)
# plotlines(df, 'grey', ax[0][1], alpha=0.1)

# plotlines(dfFrag2, 'r', ax[0][2], linewidth=3)
# plotlines(df, 'grey', ax[0][2], alpha=0.1)
# plotlines(dfFrag3, 'r', ax[0][3], linewidth=3)
# plotlines(df, 'grey', ax[0][3], alpha=0.1)
                    
# ax[1][0].scatter(theta1, rho1)
# ax[1][1].scatter(theta2, rho2)
# ax[1][2].scatter(theta3, rho3)
# ax[1][3].scatter(theta4, rho4)
# ax[1][0].set_xlabel('Rho (m) ')
# for i in range(4):
#     ax[1][i].set_xlabel('Theta (deg)')

# ax[0][0].set_title('Complete Info')
# ax[0][1].set_title('Short Dikes only')
# ax[0][2].set_title('Few Dikes')
# ax[0][3].set_title('Few Short Dikes')


# dikelength=500000

# center=np.array([0,0])
# ndikes=100
# angles=np.random.normal(45, 2, ndikes)
# m=np.tan(angles)
# Xstart=np.random.randint(0,high=dikelength, size=ndikes)
# Ystart=np.random.randint(0,high=dikelength, size=ndikes)
# Xend=dikelength/np.sqrt(1+m**2)+Xstart
# Yend=m*Xend

# df=pd.DataFrame({'Xstart':Xstart, 'Xend': Xend, 'Ystart': Ystart, 'Yend':Yend, 'Slope':m})

# fig, ax =plt.subplots(2,2)
# plotlines(df, 'k', ax[0][0])

# dfFrag=fragmentDikes(df)

# theta1,rho1, xc, yc=AKH_HT(df.astype(float))
# theta2, rho2, xc, yc=AKH_HT(dfFrag.astype(float))

# plotlines(dfFrag, 'r', ax[0][1])
# ax[1][0].scatter(rho1, theta1)
# ax[1][1].scatter(rho2, theta2)
