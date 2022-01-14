#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 19 11:50:49 2021

@author: akh
"""
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
from openmod import *
from entrainmod import *
from pltfunc import *
from curv import curvat
from normalize import labelparam
import os 

def paramlists(labels, Xval, Yval):
    Y = []
    X = []
    for sim in labels:
        param = labelparam(sim)
        Y.append(param.at[0, Yval])
        X.append(param.at[0, Xval])
    return X, Y

path = "/home/akh/myprojects/channelized-pdcs/graphs/processed/"
os.chdir('/home/akh/myprojects/channelized-pdcs/graphs/')

alllabels = [
    # 'AVX4', 'AVZ4', 'AVY4', 'AWX4', 'AWZ4', 'AWY4', 
    # 'BVX4', 'BVZ4', 'BVY4', 'BWY4', 'BWX4', 'BWZ4', 
    # 'CVX4', 'CVZ4', 'CWY4', 'CVY4', 'CWX4', 'CWZ4',
    # 'DVX4', 'DVY4', 'DVZ4', 'DWX4', 'DWY4', 'DWZ4',
    'AVX7', 'AVZ7', 'AVY7', 'AWY7', 'AWX7', 'AWZ7',
    'BVX7', 'BVY7', 'BVZ7', 'BWX7', 'BWY7', 'BWZ7',
    'CVX7', 'CVY7', 'CWX7', 'CVZ7', 'CWZ7', 'CWY7',
    'DVX7', 'DVY7', 'DVZ7', 'DWY7', 'DWZ7', 'DWX7']

IMAX=400
ZMAX=300
YMAX=150
dx=3
center=ZMAX/2*3
slope=0.18
clear=50

amprat, wave = paramlists(alllabels, 'Amprat', 'Wave')
amp, vol = paramlists(alllabels, 'Amp', 'Vflux')
width, depth = paramlists(alllabels, 'Width', 'Depth')
inlet, inletrat = paramlists(alllabels, 'Inlet', 'Inletrat')



tot, avulsed, buoyant, massout, area, areaout= openmassdist(alllabels, path)


avul = []
for sim in alllabels:
    avul.append(avulsed[sim].max())

mass_norm=avul/np.array(width)*np.array(wave)



curv=np.array(amp)/np.array(width)
cdict=pd.DataFrame({"Labels":alllabels, "Curv":curv, "amprat": amprat, "wave":wave, "width":width, 
                    "mass":mass_norm})
cdict=cdict.sort_values(by=["mass"], ignore_index=True)

def setupFig():
    fig,[ax, ax2]=plt.subplots(1,2, figsize=[16,9], dpi=300)
    ax.set_xlim([0,900])
    ax.set_ylim([0,1200])
    ax.tick_params(
        axis='x',          # changes apply to the x-axis
        which='both',      # both major and minor ticks are affected
        bottom=False,      # ticks along the bottom edge are off
        top=False,         # ticks along the top edge are off
        labelbottom=False)
    
    ax.tick_params(
        axis='y',          # changes apply to the x-axis
        which='both',      # both major and minor ticks are affected
        left=False,      # ticks along the bottom edge are off
        right=False,         # ticks along the top edge are off
        labelleft=False)
    
    # ax2.tick_params(
    #     axis='x',          # changes apply to the x-axis
    #     which='both',      # both major and minor ticks are affected
    #     bottom=False,      # ticks along the bottom edge are off
    #     top=False,         # ticks along the top edge are off
    #     labelbottom=False)
    
    ax2.tick_params(
        axis='y',          # changes apply to the x-axis
        which='both',      # both major and minor ticks are affected
        left=False,      # ticks along the bottom edge are off
        right=False,         # ticks along the top edge are off
        labelleft=False)
    plt.tight_layout()
    ax2.set_ylim([0,1.5])
    ax.set_xlabel("Channel Curvature")
    
    return fig,ax,ax2

font = {'family' : 'Helvetica',
        'size'   :18}
matplotlib.rc('font', **font)

for l in range(len(cdict["Labels"])):
    fig, ax, ax2=setupFig()
    amp= cdict['wave'].iloc[l]*cdict['amprat'].iloc[l]/dx
    aspect = 8;

    width = cdict['width'].iloc[l]
    
    #if W> l and l > 0: 
        #continue
    i=np.arange(1,IMAX*3)
    sinuous1 = amp*np.sin(np.deg2rad(360*((i)/cdict['wave'].iloc[l])))+ center-width/2
    sinuous2 = amp*np.sin(np.deg2rad(360*((i)/cdict['wave'].iloc[l])))+ center+width/2
    ax.plot(sinuous1,i, 'g')
    ax.plot(sinuous2,i, 'g')
    ax2.bar(['Mass Overspilled'], cdict["mass"].iloc[l], color='orange')
    name="topoanimation"+str(l)+".png"
    fig.savefig(name)
    plt.close()
#            for i in np.arange(1,IMAX):
#                if l == 0:
#                    sinuous= center;
#                else:
#                    sinuous = amp*np.sin(np.deg2rad(360*((i*3)/l)))+ center;


#                for z in np.arange(1,ZMAX):
#                    hill=slope*dx*(IMAX-i)+clear
                    
#                    if z >= sinuous-width/2 and z<=sinuous+width/2:
#                        Y[z,i]=hill-depth
#                    else:
#                        Y[z,i]=hill
#            curv=4 * ((np.pi)**2) * amp / (l+0.0000000001)
# fig = plt.figure()
# ax = plt.axes(projection='3d')
# x=np.arange(0,IMAX); z=np.arange(0,ZMAX)
# Z,X=np.meshgrid(x,z)
# ax.plot_surface(X,Z,Y, rstride=1, cstride=1,
#                 cmap='viridis', edgecolor='none')

