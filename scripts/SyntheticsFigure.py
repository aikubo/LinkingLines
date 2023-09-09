#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun  1 10:11:15 2022

@author: akh
"""
from synthetic import * 
from PrePostProcess import DikesetReProcess
from plotmod import * 

"""First Synthetic Figure"""

SMALL_SIZE = 8
MEDIUM_SIZE = 8
BIGGER_SIZE = 12

plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
plt.rc('figure', titlesize=MEDIUM_SIZE)  # fontsize of the figure title
plt.rcParams['legend.title_fontsize'] = SMALL_SIZE


fig,ax=plt.subplots(5,2)
fig.set_size_inches((100/25.4, 170/25.4))

kw = dict(ha="center", va="center", fontsize=12, color="black", fontstyle='italic')

l=50000
df1=makeLinearDf(5000, 30, 2, 1000, 500, ndikes= 20, label=1)
#df1=df1.append(makeLinear2(1000,90, 2, -1000, 500,ndikes=20, label=2))
#df1['Label']=df1['Label'].astype(str)
DotsLines(df1, ColorBy='Label', Cbar=False, linewidth=2,CbarLabels=False, fig=fig, ax=ax[0,:2], cmap='viridis')
#ax[0,0].text(0.1, 0.1, "Linear", transform=ax[0,0].transAxes, **kw)

l=10000
df2=makeLinear2(l, 75, 2, 0, 1000, ndikes=10, label=1)
df2=df2.append(makeLinear2(l, 30, 2, 0, 700, ndikes=10, label=2))
df2=df2.append(makeLinear2(l, -30, 2, 0, 500, ndikes=10, label=3))
for i in ['Xstart', 'Ystart', 'Xend', 'Yend']:
    df2[i]=df2[i].values+10000
#df2['Label']=df2['Label'].astype(str)
df2=DikesetReProcess(df2)
DotsLines(df2, ColorBy='Label', Cbar=False,linewidth=2, CbarLabels=False, fig=fig, ax=ax[1,:2], cmap='viridis')
#ax[1,0].text(0.1, 0.1, "Ray Radial", transform=ax[0,1].transAxes, **kw)



df3=makeCircumfrentialSwarmdf(10000, lenf=100, ndikes=20, center=[10000,10000], label=1)

df3=DikesetReProcess(df3, xc=0, yc=0)
DotsLines(df3, ColorBy='Label', Cbar=False, CbarLabels=False, fig=fig, ax=ax[3,:2], cmap='viridis')


l=5000
df3=makeRadialSwarmdf(l, center=[10000,10000], ndikes=10)
df3=df3.append(makeRadialSwarmdf(l, center=[15000,15000], label=2, ndikes=10))
df3=df3.append(makeRadialSwarmdf(l, center=[20000,20000], label=3, ndikes=10))

#df3['Label']=df3['Label'].astype(str)
df3=DikesetReProcess(df3)
DotsLines(df3, ColorBy='Label', Cbar=False, linewidth=4,CbarLabels=False, fig=fig, ax=ax[2,:2], cmap='viridis')
#ax[2,0].text(0.1, 0.1, "Overlaping Radial", transform=ax[0,2].transAxes, **kw)

df3=makeCircumfrentialSwarmdf(10, ndikes=5, center=[10,10], label=1)
df3=df3.append(makeCircumfrentialSwarmdf(10, ndikes=5,center=[15,15], label=2))
df3=df3.append(makeCircumfrentialSwarmdf(10, ndikes=5, center=[20,20], label=3))

df3=DikesetReProcess(df3, xc=15, yc=15)
DotsLines(df3, ColorBy='Label', Cbar=False, linewidth=4, CbarLabels=False, fig=fig, ax=ax[4,:2], cmap='viridis')



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

    if i==ax[4,1]:
        continue
    i.axes.get_xaxis().set_visible(False)
    if i == ax[0,1]:
        continue
    FixAxisAspect(ax[0,1], i)

ax[4,1].set_xlabel('Theta ($^\circ$)')

for i in ax[:,0]:
    i.set_facecolor('lightgrey')
plt.tight_layout()

labelSubplots(ax.flatten(), ['A', 'a',  'C', 'c', 'B', 'b', 'D', 'd', 'E', 'e'], fontsize=16)



fig.savefig('/home/akh/myprojects/Linking-and-Clustering-Dikes/Publication Figures/Synthetics1.pdf', dpi=600)




fig,ax=plt.subplots(2,2)
fig.set_size_inches((100/25.4, 100/25.4))

angles=np.linspace(-90,90, 50)
rhos=mysigmoid(angles, 0.06, -1000, 5000)
# df4=fromHT(angles, rhos, scale=1000, length=10000, label=1)
           
             
# df4=DikesetReProcess(df4, HTredo=False)
# DotsLines(df4, ColorBy='Label',Cbar=False, CbarLabels=False, fig=fig, ax=ax[0,2:4], cmap='viridis')
# ax[0,2].text(0.1, 0.1, "Sigmoid", transform=ax[2,0].transAxes, **kw)


df4=pd.DataFrame()

for c,l in zip([0.06, 0.15], [1,2]):
    rhos=mysigmoid(angles, c, -1000, 5000)
    df4=df4.append(fromHT(angles, rhos, label=l))
    

df4=DikesetReProcess(df4, HTredo=False)
DotsLines(df4, ColorBy='Label', Cbar=False, CbarLabels=False, fig=fig, ax=ax[0,:2], cmap='viridis')
#ax[1,2].text(0.1, 0.1, "Overlapping Sigmoid", transform=ax[2,1].transAxes, **kw)


rhos=angles*-500
rhos2=angles*1000

df4=fromHT(angles, rhos, scale=1000, length=10000, CartRange=50000, label=1)
df4=df4.append(fromHT(angles, rhos2, scale=1000, length=10000, label=2, CartRange=90000))

df4=DikesetReProcess(df4, HTredo=False)
DotsLines(df4, ColorBy='Label', Cbar=False, CbarLabels=False, fig=fig, ax=ax[1,:2], cmap='viridis')
#ax[2,2].text(0.1, 0.1, "Crossing Lines", transform=ax[2,2].transAxes, **kw)


for i in ax.flatten():
    if i == ax[0,1]:
        continue
    FixAxisAspect(ax[0,1], i)

plt.tight_layout()


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
    if i==ax[1,1]:
        continue
    i.axes.get_xaxis().set_visible(False)
    
    if i == ax[0,1]:
        continue
    FixAxisAspect(ax[0,1], i)


for i in ax[:,0]:
    i.set_facecolor('lightgrey')
plt.tight_layout()
labelSubplots(ax.flatten(), ['A', 'a', 'B', 'b'], fontsize=16)


fig.savefig('/home/akh/myprojects/Linking-and-Clustering-Dikes/Publication Figures/Synthetics2.pdf', dpi=600)