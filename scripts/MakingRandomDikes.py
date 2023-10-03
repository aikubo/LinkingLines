#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 10 15:49:33 2023

@author: akh
"""
from synthetic import * 
from PrePostProcess import DikesetReProcess
from plotmod import * 
import numpy as np 
import matplotlib.pyplot as plt
from fitRadialCenters import *

n=1000
tries=10

theta=np.random.rand(tries*n)
rho=np.random.rand(tries*n)

for i in range(tries):
    popt, pcov=curve_fit( lambda angle, xr,yr: CenterFunc(angle, xr, yr, 0, 0), theta[(i-1)*n: n*i-1], rho[(i-1)*n: (n)*i-1] )
    perr = np.sqrt(np.diag(pcov))
    residuals=rho[(i-1)*n: n*i-1]-CenterFunc(theta[(i-1)*n: n*i-1], *popt, 0,0)
    ss_res=np.sum(residuals**2)
    ss_tot=np.sum( (rho-np.mean(rho))**2)
    r_sq=1-(ss_res/ss_tot)
    print(r_sq)