#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 30 14:08:27 2021

@author: akh
"""
import pandas as pd
from htMOD import AKH_HT as HT
from htMOD import rotateData2 as Rotate
import numpy as np 
import matplotlib.pyplot as plt
import seaborn as sns 
np.random.seed(30) #20
sns.set_context("talk")
df=pd.DataFrame({ 'Xstart': np.random.randint(-100, 100, 4), 'Ystart': np.random.randint(-100, 100, 4), 
                 'Xend': np.random.randint(-100, 100, 4),
                 'Yend': np.random.randint(-100, 100, 4)})

theta,rho,xc,yc=HT(df, 0,0)
df["Rho"]=rho
df['Theta']=theta
m=-1*np.cos(np.deg2rad(df['Theta']))/np.sin(np.deg2rad(df['Theta']))
df['m']=m

def setupPlot():
    f,ax=plt.subplots(1,2)
    ax[0].set_ylabel("Y")
    ax[0].set_xlabel("X")
    ax[0].set_title("Cartesian Space")
    ax[1].set_ylabel("Rho (l)")
    ax[1].set_xlabel("Theta ($^\circ$)")
    ax[1].set_title("Hough Space")
    ax[0].axhline(linewidth=1, color='k')
    ax[0].axline((0,0), (0,1), linewidth=1, color='k')
    ax[0].plot( 0,0,"w*", markersize=20)
    f.set_size_inches(6, 4)
    
    ax[1].set_ylim([-100,100])
    ax[1].set_xlim([-90,90])
    ax[0].set_ylim([-130,130])
    ax[0].set_xlim([-130,130])
    plt.tight_layout()

    return f,ax

def rhoPlot(df,i,ax, colors):
    xp=np.cos(np.deg2rad(df.iloc[i]['Theta']))*(df.iloc[i]['Rho'])
    yp=np.sin(np.deg2rad(df.iloc[i]['Theta']))*(df.iloc[i]['Rho'])
    ax.axline( (df.iloc[i]['Xstart'], df.iloc[i]['Ystart']), slope=df.iloc[i]['m'], linestyle=':', color=np.append(colors[i], 0.7))
    ax.arrow( xc, yc, xp,yp, head_width=5, head_length=-5, ec=colors[i], fc=colors[i])
'''
f,ax=setupPlot()


frames=0


# for i in range(4):
#     ax[0].plot( [df.iloc[i]['Xstart'], df.iloc[i]['Xend']], [df.iloc[i]['Ystart'], df.iloc[i]['Yend']], color=colors[i])
#     ax[1].scatter(df.iloc[i]['Theta'], df.iloc[i]['Rho'], color=colors[i], s=400)
#     rhoPlot(df,i,ax[0],colors)
#     name="HTvisframe"+str(frames)+".png"
#     f.savefig(name)
#     frames=frames+1
    
    '''
    
f,ax=setupPlot()

x1=np.random.randint(-100,100,1)
y1=np.random.randint(-100,100,1)
x2=np.random.randint(-100,100,1)
y2=np.random.randint(-100,100,1)
rotateBy=20
colors=sns.color_palette("tab10")
frames=0
df2=pd.DataFrame({ 'Xstart': x1, 'Ystart': y1, 
                 'Xend': x2,
                 'Yend': y2})
theta,rho,xc,yc=HT(df2, 0,0)
df2["Rho"]=rho
df2['Theta']=theta
m=-1*np.cos(np.deg2rad(df['Theta']))/np.sin(np.deg2rad(df['Theta']))
df2['m']=m
#ax[1].set_xlim([-90,90])
for i in range(5): 
    print(df2['Theta'], df2['Rho'])
    ax[0].plot( [df2['Xstart'], df2['Xend']], [df2['Ystart'], df2['Yend']], color=colors[0])

    
    ax[1].scatter(df2['Theta'], df2['Rho'], color=colors[0], s=400)
    
    name="HTvis2frame"+str(frames)+".png"
    f.savefig(name, dpi=600)
    frames=frames+1
    ax[0].plot( [df2['Xstart'], df2['Xend']], [df2['Ystart'], df2['Yend']], color="grey")
    ax[1].scatter(df2['Theta'], df2['Rho'], color="grey", s=400)
    
    if i < 4:
        df2=Rotate(df2, 20, xc=0, yc=0)
        theta,rho,xc,yc=HT(df2, 0,0)
        df2["Rho"]=rho
        df2['Theta']=theta
    
ax[0].plot( [df2['Xstart'], df2['Xend']], [df2['Ystart'], df2['Yend']], color="grey")
ax[1].scatter(df2['Theta'], df2['Rho'], color="grey", s=400)


x1=df2['Xstart']
y1=df2['Ystart']
x2=df2['Xend']
y2=df2['Yend']

p=df2['Rho'].values

deltaRho=10

for i in range(5): 
    p=p-deltaRho
    
    xp=p[0]*np.cos(np.deg2rad(df2['Theta'].values))
    yp=p[0]*np.sin(np.deg2rad(df2['Theta'].values))
    
    print(p[0],xp[0],yp[0])
    
    ax[0].axline( (xp[0],yp[0]), slope=df2['m'].values, color=colors[2])
    ax[1].scatter(df2['Theta'], p, color=colors[2], s=400, marker='*')
    name="HTvis2frame"+str(frames)+".png"
    f.savefig(name, dpi=600)
    frames=frames+1
    ax[0].axline((xp[0],yp[0]), slope=df2['m'].values, color="grey")
    ax[1].scatter(df2['Theta'], p,  color="grey", s=400, marker='*')
    
    
    #rhoPlot(df,0,ax[0],colors)
    
    
    
    
    
