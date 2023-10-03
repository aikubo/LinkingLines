# LinkingLines Package 
 # Written by aikubo 
 # Version: 2.0.0
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 21 14:07:06 2021

This module contains functions for generating synthetic dike data and performing various operations on dike datasets.
These functions are primarily used for creating, manipulating, and analyzing synthetic dike datasets for geological studies
and modeling purposes.

@author: akh
"""
import pandas as pd
from htMOD import HoughTransform
from sklearn.preprocessing import scale
import numpy as np 
from clusterMod import *
import matplotlib.pyplot as plt
from plotmod import plotlines, labelcolors, BA_HT, HThist, DotsLines, labelSubplots
from examineMod import examineClusters
import seaborn as sns
from matplotlib import cm
from PrePostProcess import transformXstart, DikesetReProcess




def makeRadialSwarmdf(radius, doubled=True, anglestart=-90, anglestop=90, ndikes=50, center=[0,0], label=1, CartRange=100000):
    """
    Generate a DataFrame containing radial swarm dike data.

    Parameters:
        radius (float): The radius of the radial swarm.
        doubled (bool): Whether to create a doubled radial swarm.
        anglestart (float): The starting angle for dike generation (in degrees).
        anglestop (float): The stopping angle for dike generation (in degrees).
        ndikes (int): The number of dikes to generate.
        center (list): The center coordinates of the radial swarm.
        label (int): The label to assign to the generated dikes.
        CartRange (float): The Cartesian range to filter dikes based on coordinates.

    Returns:
        DataFrame: A DataFrame containing radial swarm dike data.
    """

    #center=np.array([0,0])
    angles=np.linspace(anglestart, anglestop, ndikes)
    m=-1/np.tan(np.deg2rad(angles))
    
    Xstart=np.zeros(ndikes)+center[0]
    Ystart=np.zeros(ndikes)+center[1]
    Xend=radius*np.cos(np.deg2rad((angles)))+center[0]#/np.sqrt(1+m**2)+center[0]
    Yend=radius*np.sin(np.deg2rad((angles)))+center[1] #Xend*m+center[1]
    rho=center[0]*np.cos(np.deg2rad(angles))+center[1]*np.sin(np.deg2rad(angles))
    if doubled:
        
        Xstart=np.append(Xstart, np.zeros(ndikes)+center[0])
        Ystart=np.append(Ystart, np.zeros(ndikes)+center[1])
        Yend=np.append(Yend, -radius*np.cos(np.deg2rad(90-(angles)))+center[0])
        Xend=np.append(Xend, -radius*np.sin(np.deg2rad(90-(angles)))+center[1])
        #theta2=
        #angles=np.append(angles, )
     
    
    df=pd.DataFrame({'Xstart':Xstart, 'Xend': Xend, 
                     'Ystart': Ystart, 'Yend':Yend})
                    # 'theta': angles, 'rho':rho})
    df=df.drop(df[ abs(df['Ystart']) > CartRange].index)
    df=df.drop(df[ abs(df['Yend']) > CartRange].index)
    df=df.drop(df[ abs(df['Xstart']) > CartRange].index)
    df=df.drop(df[ abs(df['Xend']) > CartRange].index)
    labels=[label]*len(df)
    df['Label']=labels
    return df


def makeCircumfrentialSwarmdf(radius, lenf=1, anglestart=-90, anglestop=90, ndikes=50, center=[0,0], label=1, CartRange=100000):
    """
    Generate a DataFrame containing circumferential swarm dike data.

    Parameters:
        radius (float): The radius of the circumferential swarm.
        lenf (float): The length factor for dike segments.
        anglestart (float): The starting angle for dike generation (in degrees).
        anglestop (float): The stopping angle for dike generation (in degrees).
        ndikes (int): The number of dikes to generate.
        center (list): The center coordinates of the circumferential swarm.
        label (int): The label to assign to the generated dikes.
        CartRange (float): The Cartesian range to filter dikes based on coordinates.

    Returns:
        DataFrame: A DataFrame containing circumferential swarm dike data.
    """

    #center=np.array([0,0])
    angles=np.linspace(anglestart, anglestop, ndikes)
    length=(radius/ndikes)*lenf
    #m=-1/np.tan(np.deg2rad(angles))
    t=np.tile(np.deg2rad(angles),2)
    radius=np.concatenate([np.ones(ndikes)*-radius,np.ones(ndikes)*radius])
    x0=radius*np.cos(t)+center[0]
    y0=radius*np.sin(t)+center[1]
    
    a = np.cos(t)
    b = np.sin(t)
    
    x1 = (x0 + length/2 * (-b))
    y1 = (y0 + length/2 * (a))
    x2 = (x0 - length/2 * (-b))
    y2 = (y0 - length/2 * (a))

    df=pd.DataFrame({'Xstart':x1, 'Xend': x2, 'Ystart':y1, 'Yend':y2})


    df=df.drop(df[ abs(df['Ystart']) > CartRange].index)
    df=df.drop(df[ abs(df['Yend']) > CartRange].index)
    df=df.drop(df[ abs(df['Xstart']) > CartRange].index)
    df=df.drop(df[ abs(df['Xend']) > CartRange].index)
    labels=[label]*len(df)
    df['Label']=labels
    
    return df

def addSwarms(dflist):
    """
    Combine multiple swarm DataFrames into a single DataFrame.

    Parameters:
        dflist (list of DataFrames): A list of DataFrames containing swarm dike data.

    Returns:
        DataFrame: A combined DataFrame containing swarm dike data.
    """

    for i, b in zip(dflist, range(len(dflist))):
        i['Label']=[b]*len(i)
    dfSwarm=pd.concat(dflist)

    return dfSwarm

def makeLinear2(length, angle, angleSTD, rho, rhoSTD, ndikes=100, CartRange=300000, label=None):
    """
    Generate a DataFrame containing linear dike data with angle and rho distributions.

    Parameters:
        length (float): The length of the dike segments.
        angle (float): The mean angle of dike orientation (in degrees).
        angleSTD (float): The standard deviation of the angle distribution.
        rho (float): The mean rho value of dike segments.
        rhoSTD (float): The standard deviation of the rho distribution.
        ndikes (int): The number of dikes to generate.
        CartRange (float): The Cartesian range to filter dikes based on coordinates.
        label (int or None): The label to assign to the generated dikes. If None, labels will be assigned automatically.

    Returns:
        DataFrame: A DataFrame containing linear dike data.
    """

    angles=np.random.normal(angle, angleSTD, ndikes)
    rhos=np.random.normal(rho, rhoSTD, ndikes)

    b=rhos/np.sin(np.deg2rad(angles))
    slopes=-1/(np.tan(np.deg2rad(angles))+0.000000001)
    Xstart=np.random.normal(0, rhoSTD, ndikes)
    
    Ystart=slopes*Xstart+b

    Xend=Xstart-length/np.sqrt(1+slopes**2)
    Yend=slopes*Xend+b
    
    if type(label) is int:
        labels=np.ones(ndikes)*label
    else:
        labels=np.arange(0,ndikes)
    
    df=pd.DataFrame({'Xstart':Xstart, 'Xend': Xend, 'Ystart': Ystart, 'Yend':Yend, 'theta':angles, 'rho':rhos, 'Labels':labels})
    
    df=df.drop(df[ abs(df['Ystart']) > CartRange].index)
    labels=[label]*len(df)
    df['Label']=labels
    
    return df

def makeLinearDf(length, angle, angleSTD, rho, rhoSTD, ndikes=100, CartRange=300000, label=None):
    """
    Generate a DataFrame containing linear dike data with angle and rho distributions.

    Parameters:
        length (float): The length of the dike segments.
        angle (float): The mean angle of dike orientation (in degrees).
        angleSTD (float): The standard deviation of the angle distribution.
        rho (float): The mean rho value of dike segments.
        rhoSTD (float): The standard deviation of the rho distribution.
        ndikes (int): The number of dikes to generate.
        CartRange (float): The Cartesian range to filter dikes based on coordinates.
        label (int or None): The label to assign to the generated dikes. If None, labels will be assigned automatically.

    Returns:
        DataFrame: A DataFrame containing linear dike data.
    """

    angles=np.random.normal(angle, angleSTD, ndikes)
    rhos=np.random.normal(rho, rhoSTD, ndikes)

    b=rhos/np.sin(np.deg2rad(angles))
    slopes=-1/(np.tan(np.deg2rad(angles))+0.000000001)
    Xstart=np.random.normal(0, rhoSTD, ndikes)
    
    Ystart=slopes*Xstart+b
    a=np.array([1,-1])
    c=np.random.choice(a, ndikes)
    Xend=Xstart-length/np.sqrt(1+slopes**2)*c
    Yend=slopes*Xend+b
    
    if type(label) is int:
        labels=np.ones(ndikes)*label
    else:
        labels=np.arange(0,ndikes)
    
    df=pd.DataFrame({'Xstart':Xstart, 'Xend': Xend, 'Ystart': Ystart, 'Yend':Yend, 'theta':angles, 'rho':rhos, 'Labels':labels})
    
    df=df.drop(df[ abs(df['Ystart']) > CartRange].index)
    labels=[label]*len(df)
    df['Label']=labels
    
    return df

def EnEchelonSynthetic(ndikes, angle, RhoStart, RhoSpacing, Overlap=0, CartRange=100000):
    """
    Generate a DataFrame containing en echelon synthetic dike data.

    Parameters:
        ndikes (int): The number of en echelon dikes to generate.
        angle (float): The angle of dike orientation (in degrees).
        RhoStart (float): The starting rho value.
        RhoSpacing (float): The spacing between rho values.
        Overlap (float): The overlap between en echelon dikes.
        CartRange (float): The Cartesian range to filter dikes based on coordinates.

    Returns:
        DataFrame: A DataFrame containing en echelon synthetic dike data.
    """

    angles=np.ones(ndikes)*angle #np.random.normal(angle, angleSTD, ndikes)
    RhoEnd=RhoStart+ndikes*RhoSpacing
    rhos= np.arange(RhoStart,RhoEnd, RhoSpacing)
    length=3*RhoSpacing
    b=rhos/np.sin(np.deg2rad(angles))
    slopes=-1/(np.tan(np.deg2rad(angles))+0.000000001)
    #
    Xstart=np.linspace(RhoStart,RhoEnd,ndikes) #np.random.normal(rho, rhoSTD, ndikes)
    
    Ystart=slopes*Xstart+b
    Xend=Xstart-length/np.sqrt(1+slopes**2)
    Yend=slopes*Xend+b
    l=np.sqrt( (Xstart-Xend)**2 + (Ystart-Yend)**2)
    df=pd.DataFrame({'Xstart':Xstart, 'Xend': Xend, 'Ystart': Ystart, 'Yend':Yend, 'theta':angles, 'seg_length':l, 'rho':rhos, 'Label': np.ones(ndikes)})
    
    df=df.drop(df[ abs(df['Ystart']) > CartRange].index)
    
    return df 

def fromHT(angles, rhos, scale=10000, length=10000, xc=0, yc=0, CartRange=100000, label=None, xrange=None, test=False):
    """
    Generate a DataFrame containing dike data from angles and rhos using the HT method.

    Parameters:
        angles (array-like): Array of dike angles (in degrees).
        rhos (array-like): Array of rho values.
        scale (float): Scaling factor for the generated dike coordinates.
        length (float): Length of the dike segments.
        xc (float): X-coordinate of the dike center.
        yc (float): Y-coordinate of the dike center.
        CartRange (float): The Cartesian range to filter dikes based on coordinates.
        label (int or None): The label to assign to the generated dikes. If None, labels will be assigned automatically.
        xrange (float or None): Optional X-coordinate range for dike generation.
        test (bool): Whether to perform a test to check the generated angles and rhos.

    Returns:
        DataFrame: A DataFrame containing dike data generated from HT parameters.
    """
    if len(angles) is not len(rhos):
        raise ValueError('Angles and Rhos arrays are not the same length')
        

    ndikes=len(angles)
    
    slopes=-1/(np.tan(np.deg2rad(angles))+0.000000001)

    b=rhos/(np.sin(np.deg2rad(angles))+0.000000001)+yc-xc*slopes


    if xrange is None:
        a=np.array([1,-1])
        c=np.random.choice(a, ndikes)
        Xstart= c*scale*np.random.randn(ndikes)+xc
        Xend=Xstart+ c*length*np.cos(np.arctan(slopes))
        Ystart=slopes*Xstart+b+yc
    
        Yend=Ystart+c*length*np.sin(np.arctan(slopes))#Ystart+length*np.sin(np.deg2rad((angles)))
        
    else:
        Xstart=xc-xrange/2
        Xend=xc+xrange/2
        Ystart=Xstart*slopes+b
        Yend=Xend*slopes+b
    
    
    if type(label) is int:
        labels=np.ones(ndikes)*label
    else:
        labels=np.arange(0,ndikes)
    
    l=np.sqrt( (Xstart-Xend)**2 + (Ystart-Yend)**2)
    
    df=pd.DataFrame({'Xstart':Xstart, 'Xend': Xend, 'Ystart': Ystart, 'Yend':Yend, 'Label': labels})
    
    df=DikesetReProcess(df, xc=xc, yc=yc)

    
    return df 

def fragmentDikes(df, nSegments=5):
    """
    Fragment dike segments into smaller segments.

    Parameters:
        df (DataFrame): The input DataFrame containing dike data.

    Returns:
        DataFrame: A DataFrame containing fragmented dike segments.
    """

    m=(df['Ystart'].values-df['Yend'].values)/(df['Xstart'].values-df['Xend'].values)
    df['Slope']=m
    bint=df['Ystart'].values-m*df['Xstart'].values
    dfFrag=pd.DataFrame()
    print(len(m), len(bint))
    for i in range(len(df)):
        high=max(df['Xend'].iloc[i], df['Xstart'].iloc[i])
        low=min(df['Xend'].iloc[i], df['Xstart'].iloc[i])
        xrange1=np.random.randint(low,high, size=nSegments)
        xrange2=np.random.randint(low,high, size=nSegments)
        m2=(m[i]+np.random.rand()*np.deg2rad(10))

        yrange1=m2*xrange1+bint[i]
        yrange2=m2*xrange2+bint[i]
        L=np.sqrt((xrange1-xrange2)**2+(yrange1-yrange2)**2)

        dfFrag=pd.concat( [dfFrag, (pd.DataFrame({'Xstart':xrange1, 'Xend': xrange2, 
                                           'Ystart': yrange1, 'Yend':yrange2,
                                           'seg_length':L, 'Xmid': (xrange1+xrange2)/2, 'Ymid': (yrange1+yrange2)/2 ,
                                           'Label': np.ones(nSegments)*i, 
                                           'Original Label':np.ones(nSegments)*df['Label'].iloc[i]
                                               }))], ignore_index=True)
                         
    return dfFrag

