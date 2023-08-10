#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun  1 10:11:15 2022

@author: akh
"""
from synthetic import * 
from PrePostProcess import DikesetReProcess
from plotmod import * 
fig,ax=plt.subplots(3,4)
fig.set_size_inches(12,6)
kw = dict(ha="center", va="center", fontsize=12, color="black", fontstyle='italic')

l=50000
df1=makeLinear2(1000, 30, 5, 1000, 500, ndikes= 25, label=1)
#df1=df1.append(makeLinear2(1000,90, 2, -1000, 500,ndikes=20, label=2))
#df1['Label']=df1['Label'].astype(str)
DotsLines(df1, ColorBy='Label', Cbar=False, CbarLabels=False, fig=fig, ax=ax[0,:2], cmap='viridis')
#ax[0,0].text(0.1, 0.1, "Linear", transform=ax[0,0].transAxes, **kw)

l=10000
df2=makeLinear2(l, 75, 2, 0, 1000, ndikes=15, label=1)
df2=df2.append(makeLinear2(l, 30, 2, 0, 700, ndikes=15, label=2))
df2=df2.append(makeLinear2(l, -30, 2, 0, 500, ndikes=15, label=3))
for i in ['Xstart', 'Ystart', 'Xend', 'Yend']:
    df2[i]=df2[i].values+10000
#df2['Label']=df2['Label'].astype(str)
df2=DikesetReProcess(df2)
DotsLines(df2, ColorBy='Label', Cbar=False, CbarLabels=False, fig=fig, ax=ax[1,:2], cmap='viridis')
#ax[1,0].text(0.1, 0.1, "Ray Radial", transform=ax[0,1].transAxes, **kw)


l=5000
df3=makeRadialSwarmdf(l, center=[10000,10000], ndikes=10)
df3=df3.append(makeRadialSwarmdf(l, center=[15000,15000], label=2, ndikes=10))
df3=df3.append(makeRadialSwarmdf(l, center=[20000,20000], label=3, ndikes=10))

#df3['Label']=df3['Label'].astype(str)
df3=DikesetReProcess(df3)
DotsLines(df3, ColorBy='Label', Cbar=False, CbarLabels=False, fig=fig, ax=ax[2,:2], cmap='viridis')
#ax[2,0].text(0.1, 0.1, "Overlaping Radial", transform=ax[0,2].transAxes, **kw)




angles=np.linspace(-90,90, 50)
rhos=mysigmoid(angles, 0.06, -1000, 5000)
# df4=fromHT(angles, rhos, scale=1000, length=10000, label=1)
           
             
# df4=DikesetReProcess(df4, HTredo=False)
# DotsLines(df4, ColorBy='Label',Cbar=False, CbarLabels=False, fig=fig, ax=ax[0,2:4], cmap='viridis')
# ax[0,2].text(0.1, 0.1, "Sigmoid", transform=ax[2,0].transAxes, **kw)

df3=makeCircumfrentialSwarmdf(10, ndikes=5, center=[10,10], label=1)
df3=df3.append(makeCircumfrentialSwarmdf(10, ndikes=5,center=[15,15], label=2))
df3=df3.append(makeCircumfrentialSwarmdf(10, ndikes=5, center=[20,20], label=3))

df3=DikesetReProcess(df3)
DotsLines(df3, ColorBy='Label', Cbar=False, CbarLabels=False, fig=fig, ax=ax[0,2:4], cmap='viridis')

df4=pd.DataFrame()

for c,l in zip([0.06, 0.15], [1,2]):
    rhos=mysigmoid(angles, c, -1000, 5000)
    df4=df4.append(fromHT(angles, rhos, label=l))
    

df4=DikesetReProcess(df4, HTredo=False)
DotsLines(df4, ColorBy='Label', Cbar=False, CbarLabels=False, fig=fig, ax=ax[1,2:4], cmap='viridis')
#ax[1,2].text(0.1, 0.1, "Overlapping Sigmoid", transform=ax[2,1].transAxes, **kw)


rhos=angles*500
rhos2=angles*1000

df4=fromHT(angles, rhos, scale=1000, length=10000, CartRange=50000, label=1)
df4=df4.append(fromHT(angles, rhos2, scale=1000, length=10000, label=2, CartRange=90000))

df4=DikesetReProcess(df4, HTredo=False)
DotsLines(df4, ColorBy='Label', Cbar=False, CbarLabels=False, fig=fig, ax=ax[2,2:4], cmap='viridis')
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

    
    if i == ax[0,1]:
        continue
    FixAxisAspect(ax[0,1], i)


for i in ax[:,0]:
    i.set_facecolor('lightgrey')
for i in ax[:,2]:
    i.set_facecolor('lightgrey')
labelSubplots(ax.flatten(), ['A', 'a',  'D', 'd', 'B', 'b', 'E', 'e', 'C', 'c',  'F', 'f'], fontsize=18)


fig.savefig('/home/akh/myprojects/Linking-and-Clustering-Dikes/Publication Figures/SyntheticsSept.png', dpi=600)
fig.savefig('/home/akh/myprojects/Linking-and-Clustering-Dikes/Publication Figures/SyntheticsSept.pdf', dpi=600)