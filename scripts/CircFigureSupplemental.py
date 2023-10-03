#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 27 14:31:21 2023

@author: akh
"""
from synthetic import * 
from PrePostProcess import DikesetReProcess
from plotmod import * 
fig,ax=plt.subplots(4,2)
fig.set_size_inches(6,8)
for i in ax[:,0]:
    i.set_facecolor('lightgrey')
df3=makeCircumfrentialSwarmdf(10000, lenf=100, ndikes=20, center=[10000,10000], label=1)

df3=DikesetReProcess(df3, xc=0, yc=0)
DotsLines(df3, ColorBy='Label', Cbar=False, CbarLabels=False, fig=fig, ax=ax[0,:], cmap='viridis')



df3=makeCircumfrentialSwarmdf(10000, ndikes=5, center=[10000,10000], label=1)

df3=DikesetReProcess(df3, xc=0, yc=0)
DotsLines(df3, ColorBy='Label', Cbar=False, CbarLabels=False, fig=fig, ax=ax[1,:], cmap='viridis')


df3=makeCircumfrentialSwarmdf(10000, ndikes=5, center=[10000,10000], label=1)
df3=df3.append(makeCircumfrentialSwarmdf(10000, ndikes=5,center=[15000,15000], label=2))
df3=df3.append(makeCircumfrentialSwarmdf(10000, ndikes=5, center=[20000,20000], label=3))

df3=DikesetReProcess(df3, xc=0, yc=0)
DotsLines(df3, ColorBy='Label', Cbar=False, CbarLabels=False, fig=fig, ax=ax[2,:], cmap='viridis')


df3=makeCircumfrentialSwarmdf(10000, ndikes=5, center=[10000,10000], label=1)
df3=df3.append(makeCircumfrentialSwarmdf(20000, ndikes=5,center=[10000,10000], label=2))
df3=df3.append(makeCircumfrentialSwarmdf(30000, ndikes=5, center=[10000,10000], label=3))

df3=DikesetReProcess(df3, xc=0, yc=0)
DotsLines(df3, ColorBy='Label', Cbar=False, CbarLabels=False, fig=fig, ax=ax[3,:], cmap='viridis')
for i in ax.flatten():
    if i == ax[0,1]:
        continue
    FixAxisAspect(ax[0,1], i)



for i in ax.flatten():
    i.tick_params(
        axis='x',          # changes apply to the x-axis
        which='both',      # both major and minor ticks are affected
        bottom=False,      # ticks along the bottom edge are off
        top=False,         # ticks along the top edge are off
        labelbottom=False)
    i.tick_params(
        axis='y',          
        which='both',      
        right=False,      
        left=False,         
        labelleft=False) 

    
    if i == ax[0,1]:
        continue
    FixAxisAspect(ax[0,1], i)


for i in ax[:,0]:
    i.set_facecolor('lightgrey')
plt.tight_layout()
labelSubplots(ax.flatten(), ['A', 'a',  'D', 'd', 'B', 'b', 'E', 'e'], fontsize=18)

fig.savefig('/home/akh/myprojects/Linking-and-Clustering-Dikes/Publication Figures/CircFig.png', dpi=600)
fig.savefig('/home/akh/myprojects/Linking-and-Clustering-Dikes/Publication Figures/CircFig.pdf', dpi=600)

