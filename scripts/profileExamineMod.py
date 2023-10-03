#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 11 14:32:30 2022

@author: akh
"""
import cProfile
from examineMod import overlapSegments
from fitRectangle import endpoints2 
from synthetic import EnEchelonSynthetic
import numpy as np
X =EnEchelonSynthetic(10, 45, 100, 5)
x,y=endpoints2(X)
size=5
Xstart=x[0:size]
Ystart=y[0:size]
Xend=x[size:size*2]
Yend=y[size:size*2]
slope=np.mean( (Ystart-Yend)/(Xstart-Xend))

from viztracer import VizTracer

#cProfile.run('')

def overlapSegments2(Xstart, Ystart, Xend, Yend, slope):
    """
    Calculate overlap between lines 
    

    """
    
    xs=[np.arange(min(x,y), max(x,y), 1) for x,y in zip(Xstart,Xend)]
    #totalL=int(np.sum( np.sqrt( (Xstart-Xend)**2 + (Ystart-Yend)**2)))
    
    arr=np.zeros( (len(xs),)+xs[0].shape)
    for i,v in enumerate(xs):
        arr[i]=np.floor(v).astype(np.uint8)
    
    u,xcounts=np.unique(arr, return_counts=True)
    
    overlapx=np.sum(xcounts[xcounts>1])
    overlapy=slope*overlapx
    
    overlap=np.sqrt(overlapx**2 +overlapy**2)
    #for i in len(Xstart):
        
    
    return overlap



#cProfile.run('overlapSegments2(X)')
import timeit 
 
def runTimes(): 
    for n in [2,5,10,15]:
        X =EnEchelonSynthetic(n, 45, 100, 5)
        x,y=endpoints2(X)
        size=n
        Xstart=x[0:size]
        Ystart=y[0:size]
        Xend=x[size:size*2]
        Yend=y[size:size*2]
        slope=np.mean( (Ystart-Yend)/(Xstart-Xend))
        print("For", n, "segments")
        #%timeit overlapSegments(X)
        #%timeit overlapSegments2(Xstart, Ystart, Xend, Yend, slope)
        print("              ")
        
x=np.array([2,5,10,15])
y1=[152, 183, 197, 211]
y2=[50.7, 64.1, 77.8, 93.9]
y1er=[6.33, 28.8, 24.7, 1.36]
y2er=[2.84, 9.45,4.46, 4.5]
y3=[38.5, 48.8, 66.6, 81.8]
y3er=[1.2, .529, 1.93, 1]
fig,ax=plt.subplots()
ax.errorbar(x,y1,y1er, label="Original")
ax.errorbar(x,y2,y2er, label="First Optimization")
ax.errorbar(x,y3,y3er, label="Second Optimization")
ax.set_ylabel('Time $\micro s$')
ax.set_ylabel('Time $\mu s$')
ax.set_ylabel('Time ($\mu s$)')
ax.set_xlabel('Number of Segments')
ax.legend(loc="upper left")