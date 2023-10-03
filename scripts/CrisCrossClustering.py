#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 26 15:23:51 2022

@author: akh
"""
from synthetic import fromHT
from plotmod import * 
from clusterMod import * 
from examineMod import * 


# t=[-85,85,-85,85]
# r=[100,100, -100, -100]
# df=fromHT(t,r, scale=100, length=100)

# DotsLines(df, ColorBy='Label')

#How do i make the yellow and cyan cluster together and the purple and red?
#Cluster 1
#[-0.100, 0.100] [-85, 85]  labels=[1,2]
# cluster 2
#[-0.100, 0.100] [85, -85] labels=[0,3]

def ToroidalDist(u,v, xbound, ybound):
    """
    Based on https://stackoverflow.com/questions/4940636/how-can-i-calculate-the-distance-between-two-points-in-cartesian-space-while-res

    Parameters
    ----------
    u : TYPE
        DESCRIPTION.
    v : TYPE
        DESCRIPTION.
    xbound : TYPE
        DESCRIPTION.
    ybound : TYPE
        DESCRIPTION.

    Returns
    -------
    TYPE
        DESCRIPTION.

    """
    raw_dx=abs( (u[0]+90)-(v[0]+90)) #angle treatment, wraps between 0-180 or -90,90
    raw_dy=abs(u[1]-v[1])
    print(raw_dx, raw_dy)
    if raw_dx < xbound/2:
        dx=raw_dx
    else:
        dx=xbound-raw_dx
    if raw_dy < ybound/2:
        dy=raw_dy
    else:
        dy=ybound-raw_dy
    print(dx,dy)
    #dx= [i if i<xbound/2 else xbound-i for i in raw_dx]
    #dy= [i if i<ybound/2 else ybound-i for i in raw_dy]
    
    return np.sqrt(dx**2+dy**2)

def TwistToroidalDistTest(u,v):
    dtheta = min(( (u[0]+90)-(v[0]+90)) % 180, ( (v[0]+90)-(u[0]+90)) % 180)
    drho= (np.sign(u[0])+int(u[0]==0))*u[1]-(np.sign(v[0])+int(v[0]==0))*v[1] #u[1]/(np.sin(np.deg2rad(u[0]))+0.000000001)- v[1]/(np.sin(np.deg2rad(v[0]))+0.000000001)
    return np.sqrt( (drho)**2+ (dtheta)**2)

def ScaledTwistToroidalDistTest(u,v):
    trho=1
    ttheta=2
    
    dtheta = min(( (u[0]+90)-(v[0]+90)) % 180, ( (v[0]+90)-(u[0]+90)) % 180)
    #drho=u[1]-v[1]
    drho= (np.sign(u[0])+int(u[0]==0))*u[1]-(np.sign(v[0])+int(v[0]==0))*v[1]
    #print(dtheta, drho)
    #if np.sum(np.sign([u[0]+0.000000001,v[0]+0.000000001]))==0:
    #drho= np.sign(u[0])*u[1]-np.sign(v[0])*v[1] #u[1]/(np.sin(np.deg2rad(u[0]))+0.000000001)- v[1]/(np.sin(np.deg2rad(v[0]))+0.000000001)
        #print('changed drho', drho)
    return np.sqrt( (drho)**2+ (dtheta/2)**2)


def tryTwistDist():
    u=[0, 100]
    v=[10, 100]
    xbound=180
    ybound=2000
    if TwistToroidalDistTest(u, v) == 10:
        print("Test 1 Passed, angle dist")
    else:
        print("Test 1 failed, angle dist")
        
    u=[0, 100]
    v=[0, 100]
    if TwistToroidalDistTest(u, v) == 0:
        print("Test 2 Passed, 0 dist")
    else:
        print("Test 2 failed, 0 dist")
        
    u=[0, 100]
    v=[0, 0]
    if TwistToroidalDistTest(u, v) == 100:
        print("Test 3 Passed, rho dist")
    else:
        print("Test 3 failed, rho dist")
    
    u=[-85, 0]
    v=[85, 0]
    if TwistToroidalDistTest(u, v) == 10:
        print("Test 4 Passed, angle crossing dist")
    else:
        print("Test 4 failed, angle crossing dist")
    
    u=[5, 1000]
    v=[5, -1000]
    if TwistToroidalDistTest(u, v) == 2000:
        print("Test 5 Passed, rho crossing dist 1")
    else:
        print("Test 5 failed, rho crossing dist 1")
            
        
    u=[-90, 1000]
    v=[90, -1000]
    if TwistToroidalDistTest(u, v) == 0:
        print("Test 6 Passed, rho + angle crossing dist 1" )
    else:
        print("Test 6 failed, rho + angle crossing dist 1")
        
    u=[-85, 1000]
    v=[90, -1000]
    if TwistToroidalDistTest(u, v) == 5:
        print("Test 7 Passed, rho + angle crossing dist 2" )
    else:
        print("Test 7 failed, rho + angle crossing dist 2")
            
    u=[-5, 1000]
    v=[5, 1000]
    if TwistToroidalDistTest(u, v) == np.sqrt(10**2+2000**2):
        print("Test 8 Passed, rho crossing dist 3")
    else:
        print("Test 8 failed, rho crossing dist 3")
        
        
    u=[5, 1000]
    v=[5, -1000]
    if TwistToroidalDistTest(u, v) == 2000:
        print("Test 9 Passed, rho crossing dist 3")
    else:
        print("Test 9 failed, rho crossing dist 3")
            
                    
tryTwistDist()
""" Synthetic Tests
theta1=np.random.normal(0,5, 10)
rho1=np.random.normal(0,100, 10)

df=fromHT(theta1, rho1, scale=1, length=1000)
DotsLines(df, ColorBy='theta')
metric= TwistToroidalDist
X=(np.vstack((theta1, rho1)).T)
M= pdist(X, metric)
threshold=1
Z=sch.complete(M)
labels=sch.fcluster(Z, t=threshold, criterion='distance')
#rootnode, nodelist=sch.to_tree(Z)
df['Labels']=labels
DotsLines(df, ColorBy='Labels')
"""


import pandas as pd

theta=[-90,90,-90,90, -1,.1,-1.1,0.5]
rho=[100,100,-100,-100, 100,-100,-100,-100]
df=fromHT(theta,rho, scale=1, length=10000)
df=DikesetReProcess(df, HTredo=False)
df, Z=HT_AGG_custom(df, 2, 100, linkage='complete', metric='Twist')
lines,evaluation=examineClusters(df, xc=0, yc=0)
DotsLines(df, ColorBy='Labels')
dtheta=2
drho=100

X=(np.vstack((df['theta'], df['rho'])).T)


X=X /[1, drho]
metric= ScaledTwistToroidalDistTest
M= squareform(pdist(X, metric))


df = pd.read_csv('/home/akh/myprojects/Linking-and-Clustering-Dikes/dikedata/deccandata/NarmadaTapi_preprocesed.csv')
l= pd.read_csv('/home/akh/myprojects/Linking-and-Clustering-Dikes/dikedata/deccandata/NarmadaTapi_Complete_2_1730.csv')
xc,yc=HT_center(df)
labels=l[l['ClusterCrossesZero']==1]['Label'].values
dfSub=df #.loc[df['Labels'].isin(labels)]
dfSub=DikesetReProcess(dfSub)
print("test subset is ", len(dfSub), "segments")
DotsLines(dfSub)
drho=df['seg_length'].mean()
dtheta=2

dfSub, Z=HT_AGG_custom(dfSub, dtheta, drho, linkage='complete', metric='Twist')
lines,evaluation=examineClusters(dfSub, xc=xc, yc=yc)

DotsLines(lines, ColorBy='R_Width')

# d2, Z=HT_AGG_custom(dfSub, dtheta, drho, linkage='complete')
# lines2,evaluation=examineClusters(d2)

# DotsLines(lines2, ColorBy='R_Width')

# checkAllClusterChange(lines,lines2)

# df = pd.read_csv('/home/akh/myprojects/Linking-and-Clustering-Dikes/dikedata/deccandata/NarmadaTapi_preprocesed.csv')
# print("test subset is ", len(df), "segments")
# DotsLines(df)
# drho=df['seg_length'].mean()

# dikeset, Z=HT_AGG_custom(df, dtheta, drho, linkage='complete', metric='Twist')
# lines,evaluation=examineClusters(dikeset)

# DotsLines(lines, ColorBy='R_Width')

# dikeset, Z=HT_AGG_custom(df, dtheta, drho, linkage='complete')
# lines2,evaluation=examineClusters(dikeset)

# DotsLines(lines2, ColorBy='R_Width')

# checkAllClusterChange(lines,lines2)