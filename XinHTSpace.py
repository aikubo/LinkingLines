#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr  1 18:58:10 2022

@author: akh
"""

import numpy as np
from synthetic import * 
from plotmod import *

angles=np.random.normal(0, scale=30, size=100)
rho=np.arctan(angles)*10000
anglesLin=np.linspace(-90,90, 100)

df1=fromHT(angles, rho)
fig1,ax1=DotsLinesHT(df1)



rho=angles*100
df1=fromHT(angles, rho)
fig2,ax2=DotsLinesHT(df1, ColorBy='theta')

rho2=np.append(rho, angles*-100)
angles2=np.append(angles, angles)
df1=fromHT(angles2, rho2)
fig2,ax2=DotsLinesHT(df1, ColorBy='theta')

def mysigmoid(x, c=1):
    return 1/(1+np.e**( -c*(x-c)))

n=1000
rho=mysigmoid(anglesLin, c=0.1)*n-n/2

df1=fromHT(anglesLin, rho, xrange=1000)
fig2,ax2=DotsLinesHT(df1, ColorBy='rho')


rho=1000*np.cos( np.deg2rad(angles))+1000*np.sin(np.deg2rad(angles))

df1=fromHT(angles, rho)
fig2,ax2=DotsLinesHT(df1, ColorBy='theta')



rho=angles*100
df11=fromHT(angles, rho)

rho=angles*500
df12=fromHT(angles, rho)

label=[1]*len(df11)+[2]*len(df12)
straightlines=df11.append(df12)
straightlines['label']=label
fig2,ax2=DotsLinesHT(straightlines, ColorBy='label', cmap='viridis')


