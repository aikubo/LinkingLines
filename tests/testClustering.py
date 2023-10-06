# #!/usr/bin/env python3
# # -*- coding: utf-8 -*-
# """
# Created on Tue Feb  1 16:27:58 2022

# @author: akh
# """

def testValidity(lines, dtheta, drho):
    
    if any(lines['ThetaRange'] > dtheta):
        print("Failed Theta Validity Check 1")
    else: 
        print("Passed Theta Validity Check 1")
        
    if any(lines['RhoRange'] > drho):
        print("Failed Theta Validity Check 1")
    else: 
        print("Passed Theta Validity Check 1")

def Run3Times(dikeset, dtheta, drho, plot=False, **kwargs):
    xc1,yc1=HT_center(dikeset)
    # run the first time 
    dikeset, Z1=HT_AGG_custom(dikeset, dtheta, drho, **kwargs)
    lines1, IC1=examineClusters(dikeset)
    
    # rotate 45 deg 
    dikeset2=rotateData2(dikeset, 45)
    theta2, rho2, xc2, yc2=AKH_HT(dikeset2)
    dikeset2['theta']=theta2
    dikeset2['rho']=rho2
    dikeset2, Z2=HT_AGG_custom(dikeset2, dtheta, drho, **kwargs)
    dikeset2=rotateData2(dikeset2,-45)
    lines2, IC2=examineClusters(dikeset2)
    
    #move to lower left
    dx=np.min([dikeset['Xstart'].min(), dikeset['Xend'].min()])
    dy=np.min([dikeset['Ystart'].min(), dikeset['Yend'].min()])
    theta3, rho3, xc3, yc3=AKH_HT(dikeset, xc=dx,yc=dy)
    dikeset3=dikeset.copy()
    dikeset3['theta']=theta3
    dikeset3['rho']=rho3
    
    dikeset3, Z3=HT_AGG_custom(dikeset3, dtheta, drho, **kwargs)
    
    lines3, IC3=examineClusters(dikeset3)
    
    #check the changes 
    eq12, diff12=checkClusterChange(lines1, lines2)
    eq13, diff13=checkClusterChange(lines1, lines3)
    eq32, diff32=checkClusterChange(lines3, lines2)
    
    Flines1=FilterLines(lines1)
    Flines2=FilterLines(lines2)
    Flines3=FilterLines(lines3)
    
    print("Comparing the filtered lines")
    eq12, diff12=checkClusterChange(Flines1, Flines2)
    eq13, diff13=checkClusterChange(Flines1, Flines3)
    eq32, diff32=checkClusterChange(Flines3, Flines2)
    
    if plot: 
        
        fig,ax=plt.subplots(1,2)
        plotlines(lines1, 'k', ax[0])
        plotlines(lines2, 'b', ax[0])
        plotlines(lines3, 'r', ax[0])
        
        ax[1].scatter(dikeset['theta'], dikeset['rho'], c=dikeset['Labels'].values, alpha=0.6, cmap=cm.Greys, marker='p')
        ax[1].scatter(dikeset2['theta'], dikeset2['rho'], c=dikeset2['Labels'].values, alpha=0.6, cmap=cm.Blues, marker='*')
        ax[1].scatter(dikeset3['theta'], dikeset3['rho'], c=dikeset3['Labels'].values, alpha=0.6, cmap=cm.Reds, marker='^')
        
        fig,ax=plt.subplots(3,2)
        
        plotlines(lines1.iloc[diff12], 'k', ax[0,0])
        plotlines(lines2.iloc[diff12], 'b', ax[0,1])
        
        plotlines(lines1.iloc[diff13], 'k', ax[1,0])
        plotlines(lines3.iloc[diff13], 'r', ax[1,1])
        
        plotlines(lines2.iloc[diff32[diff32<len(lines2)]], 'b', ax[2,0])
        plotlines(lines3.iloc[diff32], 'r', ax[2,1])
        
        
        
    return lines1, lines2, lines3

def RotateAndCluster(dikeset,dtheta, drho,**kwargs):
    dikesetOrig, ZOrig=HT_AGG_custom(dikeset, dtheta, drho)
    
    dikeset45=rotateData(dikeset, 45)
    dikeset45=DikesetReProcess(dikeset45, HTredo=True)
    dikeset45, Z45=HT_AGG_custom(dikeset45, dtheta, drho)
    linesOrig=examineClusterShort(dikesetOrig)
    lines45=examineClusterShort(dikeset45)
    
    changedect=checkAllClusterChange(linesOrig, lines45)
    if changedect:
        return dikesetOrig
    if ~changedect:
        eqLabels, diffLabels=checkIndividualClusterChange(linesOrig, lines45)
        dikesetFinal=dikesetOrig.iloc[eqLabels[:len(dikesetOrig)]]

        return dikesetFinal
    
# from scipy.cluster import hierarchy
# from scipy.spatial.distance import squareform
# import pandas as pd
# from htMOD import AKH_HT as HT
# from sklearn.preprocessing import scale
# import numpy as np 
# from clusterMod import *
# import matplotlib.pyplot as plt
# from plotmod import plotlines, labelcolors, plotbyAngleBin, BA_HT, plotResults
# from examineMod import examineClusters,extendLines
# from PrePostProcess import *
# from matplotlib import cm

# from matplotlib.lines import Line2D  

# from matplotlib import cm

# import numpy as np
# import matplotlib.pyplot as plt

# from plotmod import plotlines, HT3D

# from skimage import feature
# from htMOD import HT_center
# from findRadialCenters import sweepCenters, detectLocalPeaks
# from scipy.spatial.distance import pdist, squareform
# import scipy.cluster.hierarchy as sch

# from examineMod import persitance

# def HashLabels(df):
    
#     l=[]
#     for i in df['Labels']:
#         h=df[df['Labels']==i]['HashID'].sum()
#         l.append(h)
        
#     df['ClusterHashLabel']=l
        
#     return df

# def checkClusterHash(df1, df2):
#     return df1['ClusterHashLabel']==df2['ClusterHashLabel']

# def cluster(df, metric):
#     X=(np.vstack((df['theta'], df['ScaledRho'])).T)
#     M=squareform(pdist(X, metric=metric))
#     threshold =1 
#     linkage='single'
    
#     Y=sch.single(squareform(M))
#     clusters=sch.fcluster(Y, t=.5, criterion='distance')    
#     # ag = AGG(n_clusters=None, affinity='precomputed', distance_threshold=threshold, linkage=linkage)
#     # clusters=ag.fit(squareform(M))
#     # df['Labels']=clusters.labels_
    
#     return M, clusters, df 

# def cluster(df, metric):
#     X=(np.vstack((df['theta'], df['ScaledRho'])).T)
#     M=squareform(pdist(X, metric=metric))
#     threshold =1 
#     linkage='avera'
    
#     Y=sch.average(squareform(M))
#     clusters=sch.fcluster(Y, t=1, criterion='distance')    
#     # ag = AGG(n_clusters=None, affinity='precomputed', distance_threshold=threshold, linkage=linkage)
#     # clusters=ag.fit(squareform(M))
#     df['Labels']=clusters
    
#     return M, clusters, df


# testData=pd.read_csv('/home/akh/myprojects/Linking-and-Clustering-Dikes/dikedata/testDatafmCRB.csv')

# metric1= lambda u, v: CyclicEuclideanScaled(u,v, 2,572)
# metric2= lambda u, v: CyclicEuclideanScaled(u,v, 5,1000)
# metric3= lambda u, v: CyclicEuclideanScaled(u,v,.2,.1)
# testData=HashLabels(testData)
# trueData1=HashLabels(testData)
# h=testData['HashID'].values
# trueLabels2=pd.DataFrame({'HashID':h, 'Labels':[1,2,3,4,5,5,5,5,6,6,7,7]})
# trueLabels2=HashLabels(trueLabels2)

# trueLabels3=pd.DataFrame({'HashID':h, 'Labels':[1,2,3,4,5,6,7,8,9,9,10,11, ]})
# trueLabels3=HashLabels(trueLabels3)


# print( "Test 0: does the expected hashes line up with testData and testData")

# if all(checkClusterHash(testData, testData)):
#     print ("passed test 0")

# print('Test 1, dtheta=2, drho=572')
# M1,cluster1,df1=cluster(testData.copy(), metric1)
# df1=HashLabels(df1)

# if all(checkClusterHash(trueData1, df1)):
#     print ("passed test 1")
    
# else:
#     print('Failed test 1')
#     print(checkClusterHash(trueData1, df1))


# print('Test 2, dtheta=10, drho=1000')
# M2,cluster2,df2=cluster(testData.copy(), metric2)
# df2=HashLabels(df2)
# if all(checkClusterHash(trueLabels2, df2)):
#     print ("passed test 2")
# else:
#     print('Failed test 2')
#     print(checkClusterHash(trueLabels2, df2))
#     Z=sch.average(squareform(M2))
#     fig,ax=plt.subplots()
#     Y=sch.dendrogram(Z)
#     plt.show()

# print('Test 3, dtheta=.1, drho=10')
# M3,cluster3,df3=cluster(testData.copy(), metric3)
# df3=HashLabels(df3)
# if all(checkClusterHash(trueLabels3, df3)):
#     print ("passed test 3")
# else:   
#     print('Failed test 3')
    
#     print(checkClusterHash(trueLabels3, df3))
#     Z=sch.average(squareform(M3))
#     fig,ax=plt.subplots()
#     Y=sch.dendrogram(Z)
#     plt.show()

# from synthetic import *
# from plotmod import plotlines, DotsHT

# TrueLength=5000
# df1 = makeLinear2(TrueLength, 20, 5, 5000, 1000,10, CartRange=100000)
# df2 = makeLinear2(TrueLength, 45, 5, -5000, 1000, 10, CartRange=100000)

# fig,ax=plt.subplots(3,2)


# df = addSwarms([df1, df2])
# df['Label']=np.arange(0,len(df))
# df=completePreProcess(df)
# dikes = fragmentDikes(df)
# theta, rho, xc, yc= HT(dikes, xc=0, yc=0)
# dikes['theta']=theta
# dikes['rho']=rho
# dikes= MidtoPerpDistance(dikes, xc, yc)
# midist=dikes['PerpOffsetDist'].values

# plotlines(df, 'r', ax[0,0], ColorBy='Label')
# plotlines(dikes, 'r', ax[1,0], ColorBy='Label')
# fig,ax[0,1]=DotsHT(fig, ax[0,1], df, ColorBy='Label')
# fig,ax[1,1]=DotsHT(fig, ax[1,1], dikes, ColorBy='Label')


# dikes['TrueLabel']=dikes['Label']
# from clusterMod import HT_AGG_custom
# from examineMod import checkClusterChange2

# dtheta=1
# drho=500#dikes['seg_length'].mean()
# dikes,clusters, M=HT_AGG_custom(dikes, dtheta, drho, linkage='average')
# lines,IC=examineClusters(dikes)

# plotlines(lines, 'r', ax[2,0], ColorBy='Label')
# fig,ax[2,1]=DotsHT(fig, ax[2,1], lines, ColorBy='Label')


# #checkClusterChange2(dikes,lines, df)
# #fig,ax=plt.subplots(3)
# if len(lines)==len(df): 
#     print("Test 4.1 same number of dikes identified")    
# else: 
#     print("Test 4.1 Failed")    
        
# for ind, row in lines.iterrows():
#     print(ind)
#     terror=np.min( abs(df['theta']-row['AvgTheta']))
#     perror=np.min( abs(df['rho']-row['AvgRho']))
#     lerror=abs(5000-row['R_Length'])
#     print("theta error", terror )
#     print('rho error', perror)
#     print('length error', )
    
#     if (row['StdTheta'] > dtheta):
#         print(ind, "Test 4.2 Failed")
    
    
    
#     # little correlation between types of error
#     #ax[0].plot(perror, terror, 'r*')
#     #ax[1].plot(lerror, terror, 'r*')
#     #ax[2].plot(lerror, perror, 'r*')

# print("Mean length:", lines['R_Length'].mean(), "Max:",lines['R_Length'].max())
# if lines['R_Length'].max() > TrueLength: 
#     print('Test 4.3 Max length greater than true length')
#     print('Failed')
# else: 
#     print('Test 4.3 Max length greater than true length')
#     print("Passed")

# test5=False
# if test5:
#     print('entering test 5 ')
#     spacing=[200, 500, 1000, 5000]
#     drhos=[200,350, 500,1000, 2000, 5000, 6000, 7000]
#     cs=[]
#     cs2=[]
#     cs3=[]
#     T=[]
#     ndikes=10
#     for i in spacing:
#         df=standardSpacing(ndikes, 45, 10000, i)
#         df=completePreProcess(df)
#         dikes = fragmentDikes(df)
#         theta, rho, xc, yc= HT(dikes, xc=0, yc=0)
#         dikes['theta']=theta
#         dikes['rho']=rho
#         dikes= MidtoPerpDistance(dikes, xc, yc)
#         midist=dikes['PerpOffsetDist'].values
#         dtheta=1
#         for drho in drhos:
            
            
#             dikes,clusters, M=HT_AGG_custom(dikes, dtheta, drho, linkage='single')
#             lines,IC=examineClusters(dikes)
#             #print ("for true spacing", i, "and clustering threshold", drho)
#             #print( len(lines), "clusters idenitified of ", len(df), "true dikes")
#             cs.append(len(lines))
            
            
#             dikes,clusters, M=HT_AGG_custom(dikes, dtheta, drho, linkage='complete')
#             lines,IC=examineClusters(dikes)
#             #print ("for true spacing", i, "and clustering threshold", drho)
#             #print( len(lines), "clusters idenitified of ", len(df), "true dikes")
#             cs2.append(len(lines))
            
            
#             dikes,clusters, M=HT_AGG_custom(dikes, dtheta, drho, linkage='average')
#             lines,IC=examineClusters(dikes)
#             #print ("for true spacing", i, "and clustering threshold", drho)
#             #print( len(lines), "clusters idenitified of ", len(df), "true dikes")
#             cs3.append(len(lines))
            
#             T.append(drho/i)
    
    
#     x=np.array(T)
    
    
#     cs=np.array(cs)
#     cs2=np.array(cs2)
#     cs3=np.array(cs3)
    
#     f,ax=plt.subplots()
#     ax.plot(x, cs, "rs", label='Single Linkage', markersize=7, alpha=0.7, markeredgecolor="r" )
#     ax.plot(x, cs2, "bp", label='Complete Linkage', markersize=7, alpha=0.7, markeredgecolor="b")
#     ax.plot(x,cs3, "g*", label='Average Linkage', markersize=7, alpha=0.7, markeredgecolor="g")
    
#     #fit exponential curve to cs data using scipy curve_fit
#     import scipy.optimize as opt
    
#     def expFit(x,y):
    
    
#         def exp(x,a,b,c):
#             return a*np.exp(b*-1*x)+c
    
#         popt, pcov = opt.curve_fit(exp, x,y)
        
#         #calculate the r^2 value
#         residuals = np.sum(np.square(y - exp(x, *popt)))
#         squaredDiff = np.sum(np.square(y - np.mean(y)))
#         r2 = 1 - (residuals / squaredDiff)
#         x2=np.linspace(min(x), max(x), 100)
#         y2=exp(x2, *popt)
        
    
#         return popt, pcov, r2, x2, y2
    
#     popt3, pcov3, r2_3, x2, y2_3=expFit(x, cs3)
#     popt2, pcov2, r2_2, x2, y2_2=expFit(x, cs2)
#     popt1, pcov1, r2, x2, y2=expFit(x, cs)
    
#     ax.plot(x2, y2, "r")
#     ax.plot(x2, y2_2, "b")
#     ax.plot(x2, y2_3, "g")
    
#     ax.legend(loc='upper right')
#     ax.set_ylabel('Number of Dikes')
#     ax.set_xlabel('RhoThreshold/TrueSpacing')

# from examineMod import enEchelon
# test7=False
# if test7:
#     print('entering test 7')
    
#     AvgTheta=80
#     df=EnEchelonSynthetic(5, AvgTheta, 10000, 200)
#     df=completePreProcess(df)
#     dikes = fragmentDikes(df)
#     theta, rho, xc, yc= HT(dikes, xc=0, yc=0)
#     dikes['theta']=theta
#     dikes['rho']=rho
#     dikes= MidtoPerpDistance(dikes, xc, yc)
#     midist=dikes['PerpOffsetDist'].values
#     df['Label']=np.ones(len(df))
#     Xmid=df['Xmid'].to_numpy()
#     Ymid=df['Ymid'].to_numpy()
    
#     slope, intercept, r_value, p_value, std_err = stats.linregress(Xmid, Ymid)
#     thetaM=np.rad2deg(np.arctan(-1/(slope+ 0.0000000000000001)))
#     if np.isnan(slope):
#         thetaM=0
#         p_value=0.05
    
    
#     tdiff=CyclicAngleDist([AvgTheta], [thetaM])
#     print(tdiff, thetaM, p_value)
#     #if p_value > 0.05:
#     #    tdff=np.nan
    
#     #tdiff=enEchelon(df,45)
#     print(tdiff)
    
#     tdiff2, EXstart, EXend, EYstart, EYend=enEchelon(df, AvgTheta)
    
#     fig,ax=plt.subplots()
#     plotlines(df, 'k', ax)
#     ax.plot(Xmid, Ymid, 'r*-')
#     ax.axis('equal')
#     ax.plot([EXstart, EXend], [EYstart, EYend])


# print ("Test 8 ")

# print("Are sklearn and scipy the same?")

# from clusterMod import HT_AGG_custom as AggHT
# from clusterMod import HT_AGG_custom2
# drho=572
# dtheta=2

# test1, clusters, M=AggHT(testData, dtheta, drho)

# test2, Z=HT_AGG_custom2(testData, dtheta, drho)

# print("Sklearn", test1['Labels'])
# print("Scipy", test2['Labels'])

# from examineMod import checkClusterChange

# print(checkClusterChange(test1, test2))

