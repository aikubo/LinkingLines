#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 10 13:30:34 2021

@author: akh
"""

import numpy as np 

import matplotlib.pyplot as plt


def rotateXYShift(ang,x,y,h,k):
    xp= (x-h)*np.cos(ang)-(y-k)*np.sin(ang)
    yp= (x-h)*np.sin(ang)+(y-k)*np.cos(ang)
    return xp, yp


def unrotate(ang,x,y,h,k):
    xr= (x)*np.cos(ang)+(y)*np.sin(ang)+h
    yr= -1*(x)*np.sin(ang)+(y)*np.cos(ang)+k
    return xr, yr


def testRotate(ang,x,y,h,k):
    plt.plot(x,y, 'b')
    xp,yp=rotateXYShift(ang,x,y,h,k)
    plt.plot(xp,yp, 'g')
    xr,yr=unrotate(ang,xp,yp,h,k)
    plt.plot(xr,yr, 'r.-')
    plt.plot(h,k,'r*')

def inRectangle(a,b,xp,yp):
    xtop=b; ytop=a
    insideX=True
    insideY=True

    if any(abs(xp) > xtop):
        insideX=False
    if any(abs(yp) > ytop):
        insideY=False
        
    return insideX, insideY

def endpoints(lines):
    xlist=np.array([])
    ylist=np.array([])
    
    for i in range(0,len(lines)):
        x1=lines['Xstart'].iloc[i]
        y1=lines['Ystart'].iloc[i]
        y2=lines['Yend'].iloc[i]
        x2=lines['Xend'].iloc[i]
        
        xlist=np.append(xlist,x1)
        xlist=np.append(xlist, x2)
        ylist=np.append(ylist, y1)
        ylist=np.append(ylist, y2)
    return xlist, ylist

def endpoints2(lines):
    xstart=lines['Xstart'].to_numpy()
    ystart=lines['Ystart'].to_numpy()
    xend=lines['Xend'].to_numpy()
    yend=lines['Yend'].to_numpy()
    
    xlist=np.append(xstart,xend)
    ylist=np.append(ystart,yend)
    
    return xlist,ylist

def midpoint(lines):
    xstart=lines['Xstart'].to_numpy()
    ystart=lines['Ystart'].to_numpy()
    xend=lines['Xend'].to_numpy()
    yend=lines['Yend'].to_numpy()
    xlist=np.array([])
    ylist=np.array([])
    
    for i in range(len(xstart)):
        xs=np.mean([xstart[i], xend[i]])
        ys=np.mean([ystart[i], yend[i]])
    
        xlist=np.append(xlist, xs)
        ylist=np.append(ylist, ys)
    
    return xlist,ylist

def allpoints(lines):
    xstart=lines['Xstart'].to_numpy()
    ystart=lines['Ystart'].to_numpy()
    xend=lines['Xend'].to_numpy()
    yend=lines['Yend'].to_numpy()
    xlist=np.array([])
    ylist=np.array([])
    
    for i in range(len(xstart)):
        xs=np.linspace(xstart[i], xend[i])
        ys=np.linspace(ystart[i], yend[i])
    
        xlist=np.append(xlist, xs)
        ylist=np.append(ylist, ys)
    
    return xlist,ylist

def testRec(lines, xc, yc): 
    fig, a=plt.subplots()
    
    xi,yi=endpoints2(lines)
    x0=xc
    y0=xc
    ang=-np.deg2rad(np.mean(lines['theta']))
    xp, yp= rotateXYShift(ang, xi,yi, x0,y0)
    #plotlines(lines, 'k', a)
    
    for i in range(0,len(lines)):
        a.plot( [xp[i], xp[i+len(lines)]],  [yp[i], yp[i+len(lines)]], 'r')
    width=np.ptp(xp)
    length=np.ptp(yp)
    
    # if width>length :
    #     length=width
    #     width=length
    xc=(max(xp)-min(xp))/2 + min(xp)
    yc=(max(yp)-min(yp))/2 + min(yp)
    
    xr=xc+width/2
    xl=xc-width/2
    yu=yc+length/2
    yd=yc-length/2
    xs=np.append(xr,xl)
    ys=np.append(yu,yd)
    
    
    xpi, ypi=unrotate(ang, xp, yp, x0, y0)
    # for i in range(0,len(lines)):
    #     a.plot( [xpi[i], xpi[i+9]],  [ypi[i], ypi[i+9]], 'b')
        
    Xedges=np.array([xs[0], xs[0], xs[1], xs[1], xs[0]])
    Yedges=np.array([ys[1], ys[0], ys[0], ys[1], ys[1]])
    a.plot(Xedges, Yedges, 'r.-')
    
    xs,ys=unrotate(ang,Xedges,Yedges,x0,y0)

    a.plot(xs, ys, 'r.-')
    

def fit_Rec(lines,xc,yc):
    xi,yi=endpoints2(lines)
    x0=xc
    y0=xc
    ang=-np.deg2rad(np.mean(lines['theta']))
    xp, yp= rotateXYShift(ang, xi,yi, x0,y0)
    
    width=np.ptp(xp)
    length=np.ptp(yp)
    
    # if width>length :
    #     length=width
    #     width=length
    xr=xc+width/2
    xl=xc+-width/2
    yu=yc+length/2
    yd=yc-length/2
    xs=np.append(xr,xl)
    ys=np.append(yu,yd)
    xsi,ysi=unrotate(ang,xs,ys,x0,y0)
    
    return width, length

def RecEdges(lines, xc,yc):
    xi,yi=endpoints2(lines)
    x0=xc
    y0=xc
    ang=-np.deg2rad(np.mean(lines['theta']))

    xp, yp= rotateXYShift(ang, xi,yi, x0,y0)
    
    xc=(max(xp)-min(xp))/2
    yc=(max(yp)-min(yp))/2
    
    width=np.ptp(xp)
    length=np.ptp(yp)
    xr=xc+width/2
    xl=xc+-width/2
    yu=yc+length/2
    yd=yc-length/2
    xs=np.append(xr,xl)
    ys=np.append(yu,yd)
    xsi,ysi=unrotate(ang,xs,ys,x0,y0)
    
    return xsi, ysi


def pltLine(lines, xc, yc, ax):
    avgtheta=np.deg2rad(np.average(lines['theta']))
    avgrho=np.average(lines['rho'])
    xs,ys=allpoints(lines)
    
    m=-np.cos(avgtheta)/np.sin(avgtheta)
    b=avgrho/np.sin(avgtheta)
    
    ax.plot( [min(xs), max(xs)], [m*min(xs-xc)+b+yc, m*max(xs-xc)+b+yc], 'orange', '..')
    

def W_L(Clusters):
    width=np.array([])
    length=np.array([])
    
    labels=np.unique(Clusters['Labels'])
    for i in labels:
        mask=(Clusters['Labels']==i)
        lines=Clusters[mask]
        w,l=fit_Rec(lines)
        width=np.append(width,w)
        length=np.append(length,l)

    return width, length


def squaresError(lines, xc, yc):
    
    avgtheta=np.deg2rad(np.average(lines['theta']))
    #print(avgtheta)
    avgrho=np.average(lines['rho'])
    xs,ys=midpoint(lines)
    
    m=-np.cos(avgtheta)/np.sin(avgtheta)
    b=avgrho/np.sin(avgtheta)
    
    r=np.sum((ys-(m*(xs-xc)+b+yc))**2)/lines['seg_length'].sum() #len(lines)
    
    # just do b? 
    
    return r