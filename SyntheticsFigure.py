#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun  1 10:11:15 2022

@author: akh
"""
from synthetic import * 
from PrePostProcess import DikesetReProcess
from plotmod import * 
fig,ax=plt.subplots(3,2)
fig.set_size_inches(6,12)

l=50000
df1=makeLinear2(l, 30, 5, 10000, 100000, label=1)
#df1=df1.append(makeLinear2(l, 30, 2, -100000, 10000,ndikes=10, label=2))
df1['Label']=df1['Label'].astype(str)
DotsLines(df1, ColorBy='Label', CbarLabels=False, fig=fig, ax=ax[0,:])
FixAxisAspect(ax[0,1], ax[0,0])
l=10000
df2=makeLinear2(l, 75, 2, 0, 1000, ndikes=15)
df2=df2.append(makeLinear2(l, 30, 2, 0, 700, ndikes=15, label=2))
df2=df2.append(makeLinear2(l, -30, 2, 0, 500, ndikes=15, label=3))
for i in ['Xstart', 'Ystart', 'Xend', 'Yend']:
    df2[i]=df2[i].values+10000
df2['Label']=df2['Label'].astype(str)
df2=DikesetReProcess(df2)
DotsLines(df2, ColorBy='Label', CbarLabels=False, fig=fig, ax=ax[1,:])

FixAxisAspect(ax[1,1], ax[1,0])
l=5000
df3=makeRadialSwarmdf(l, center=[10000,10000], ndikes=10)
df3=df3.append(makeRadialSwarmdf(l, center=[15000,15000], label=2, ndikes=10))
df3=df3.append(makeRadialSwarmdf(l, center=[20000,20000], label=3, ndikes=10))
df3=df3.append(makeRadialSwarmdf(l, center=[21000,21000], label=4, ndikes=10))
df3['Label']=df3['Label'].astype(str)
df3=DikesetReProcess(df3)
DotsLines(df3, ColorBy='Label', CbarLabels=False, fig=fig, ax=ax[2,:])
FixAxisAspect(ax[2,1], ax[2,0])

fig2,ax2=plt.subplots(3,2)
fig.set_size_inches(6,12)


angles=np.linspace(-90,90, 50)
rhos=mysigmoid(angles, 0.06, -1000, 5000)
df4=fromHT(angles, rhos)
df4=DikesetReProcess(df4, HTredo=False)
DotsLines(df4, ColorBy='Label', CbarLabels=False, fig=fig, ax=ax2[0,:])
FixAxisAspect(ax2[0,1], ax2[0,0])


df4=pd.DataFrame()

for c,l in zip([0.06, 0.15], [1,2]):
    rhos=mysigmoid(angles, c, -1000, 5000)
    df4=df4.append(fromHT(angles, rhos, label=l))
    

df4=DikesetReProcess(df4, HTredo=False)
DotsLines(df4, ColorBy='Label', CbarLabels=False, fig=fig, ax=ax2[1,:])
FixAxisAspect(ax2[1,1], ax2[1,0])


rhos=angles*500
rhos2=angles*1000

df4=fromHT(angles, rhos, CartRange=50000)
df4=df4.append(fromHT(angles, rhos2, label=2, CartRange=90000))

df4=DikesetReProcess(df4, HTredo=False)
DotsLines(df4, ColorBy='Label', CbarLabels=False, fig=fig, ax=ax2[2,:])
FixAxisAspect(ax2[2,1], ax2[2,0])


