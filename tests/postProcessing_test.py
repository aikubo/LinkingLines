# #!/usr/bin/env python3
# # -*- coding: utf-8 -*-
# """
import pytest
import pandas as pd
import numpy as np
from fitRectangle import rotateXYShift

class TestRotate:
    
    @staticmethod
    def test_rotate90():
        # Test case 1: Rotate 90 degrees counterclockwise about the origin
        x, y = 1, 0
        ang = np.pi / 2
        xp, yp = rotateXYShift(ang, x, y, 0, 0)
        assert np.isclose(xp, 0)
        assert np.isclose(yp, 1)
        
    @staticmethod
    def test_rotate45():
        # Test case 2: Rotate 45 degrees clockwise about (1, 1)
        x, y = 2, 1
        ang = -np.pi / 4
        h, k = 1, 1
        xp, yp = rotateXYShift(ang, x, y, h, k)
        assert np.isclose(xp, 1 + np.sqrt(2))
        assert np.isclose(yp, 1 - np.sqrt(2))
        
# Created on Sat Sep  9 09:44:53 2023

# @author: akh
# """
# #!/usr/bin/env python3
# # -*- coding: utf-8 -*-
# """
# Created on Fri Jul 15 12:28:04 2022

# @author: akh
# """
# from examineMod import * 
# from PrePostProcess import * 
# from plotmod import * 

# import pandas as pd
# import numpy as np 
# from htMOD import rotateData2

# #!/usr/bin/env python3
# # -*- coding: utf-8 -*-
# """
# Created on Fri Jul 15 12:28:04 2022

# @author: akh
# """
# from examineMod import * 
# from PrePostProcess import * 
# from plotmod import * 
# from dilationCalculation import dilation
# import pandas as pd
# import numpy as np 
# from htMOD import rotateData2
# from synthetic import makeLinear2
# from fitRadialCenters import * 

# from synthetic import fromHT 
    


# def SimpleTestFunc(Value, Expected, TestName):
#     if type(Value) is list and type(Expected) is list and type(TestName) is list: 
#         passed=0
#         for Value1, Expected1, TestName1 in zip(Value, Expected, TestName):
            
#             if np.isclose(Expected1, Value1):
#                 print("Pass", TestName1)
#                 passed=passed+1
#             else:
#                 print("Fail", TestName1)
#                 print("Expected value", Expected1, "Returned", Value1)
#         return passed
#     else:
#         if np.isclose(Expected, Value):
#             print("Pass", TestName)
#         else:
#             print("Fail", TestName)
#             print("Expected value", Expected, "Returned", Value)
#         return np.isclose(Expected, Value)*1
    
    
# def TestDilation():
#     #df=makeLinear2(1000, 1, 0, 5000, 500, ndikes=5)
#     df2=makeLinear2(1000, 90, 0, 5000, 500, ndikes=5)
#     df3=makeLinear2(1000, 45, 0, 5000, 500, ndikes=5)
    
    
#     #DotsLines(df, ColorBy='Labels')
#     fig, ax=DotsLines(df2, ColorBy='Labels')
#     fig, ax=DotsLines(df3, ColorBy='Labels')

#     EWDilation, NSDilation, _, _=dilationPlot(df2, binWidth=1)
#     EWDilation2, NSDilation2, _, _=dilationPlot(df3, binWidth=1)

#     print("Testing 90 Degrees")
#     passed=0
#     if np.isclose(np.max(EWDilation), 0):
#         print("Pass EW 1")
#         passed=passed+1
#     else:
#         print("Fail EW 1")
#         print("Expected value", 0, "Returned", np.max(EWDilation))

#     if np.isclose(np.max(NSDilation), len(df2)):
#         print("Pass NS 1")
#         passed=passed+1
#     else:
#         print("Fail NS 1")
#         print("Expected value", len(df2), "Returned", np.max(NSDilation))

#     print("Testing 45 Degrees")
#     EW45=np.cos(np.deg2rad(45))*2 #acounts for where tips interact
#     if np.isclose(np.max(EWDilation2), EW45):
#         print("Pass EW 2")
#         passed=passed+1
#     else:
#         print("Fail EW 2")
#         print("Expected value", EW45, "Returned", np.max(EWDilation2))

#     NS45=np.sin(np.deg2rad(45))*5
#     if np.isclose(np.max(NSDilation2), NS45):
#         print("Pass NS 2")
#         passed=passed+1
#     else:
#         print("Fail NS 2")
#         print("Expected value", NS45, "Returned", np.max(NSDilation2))
        
        
        
#     print('Testing non-1 averageWidth')
#     EWDilation, NSDilation, _, _=dilationPlot(df2, binWidth=1, averageWidth=10)
#     Value=[np.max(NSDilation), np.max(EWDilation)]
#     Expected=[10*len(df2), 0]
#     TestName=['NS 3', 'EW 3']
#     p=SimpleTestFunc(Value, Expected, TestName)
        
#     passed=passed+p
    
#     EWDilation, NSDilation, _, _=dilationPlot(df3, binWidth=1, averageWidth=10)
#     Value=[np.max(NSDilation), np.max(EWDilation)]
#     Expected=[NS45*10, EW45*10]
#     TestName=['NS 4', 'EW 4']
#     p=SimpleTestFunc(Value, Expected, TestName)
        
#     passed=passed+p
    
    
#     print('Testing Expanded Method')
#     EWDilation, NSDilation, _, _=dilationPlot(df2, binWidth=1, averageWidth=10, method='Expanded')
#     Value=[np.max(NSDilation), np.max(EWDilation)]
#     Expected=[10*len(df2), 0]
#     TestName=['NS 5', 'EW 5']
#     p=SimpleTestFunc(Value, Expected, TestName)
        
#     passed=passed+p
    
#     EWDilation, NSDilation, _, _=dilationPlot(df3, binWidth=1, averageWidth=10, method='Expanded')
#     Value=[np.max(NSDilation), np.max(EWDilation)]
#     Expected=[NS45*10, EW45*10]
#     TestName=['NS 6', 'EW 6']
#     p=SimpleTestFunc(Value, Expected, TestName)
#     passed=passed+p
#     total=12
#     print( passed, "of total tests", total)

# ''' Test Overlap Calculation'''
# def RotateOverlapTest(lines):
#     theta=np.mean(lines['theta'].values)
    
#     dfRotated=rotateData2(lines, (90-theta))
#     dfRotated=transformXstart(dfRotated)
#     fig,ax=plt.subplots()

#     plotlines(lines, 'k', ax)
#     plotlines(dfRotated, 'r', ax)
    
#     Xstart=dfRotated['Xstart'].to_numpy()
#     Ystart=dfRotated['Ystart'].to_numpy()
#     Xend=dfRotated['Xend'].to_numpy()
#     Yend=dfRotated['Yend'].to_numpy()
#     step=1
#     totalL=np.sum( np.sqrt( (Xstart-Xend)**2 + (Ystart-Yend)**2))
#     xs=[np.arange(min(x,y), max(x,y), step) for x,y in zip(np.floor(Xstart),np.ceil(Xend))]
    
#     l=np.max([len(xi) for xi in xs])
#     xs_sameLength=[np.append(xi,[np.nan]*(l-len(xi))) for xi in xs]

#     arr=np.vstack(xs_sameLength) #better
#     u,xcounts=np.unique(arr[~np.isnan(arr)], return_counts=True)
    
#     overlapx=np.float64(np.sum(xcounts[xcounts>1]-1)*step)
    
#     print(overlapx)
#     print(totalL)
#     overlap=overlapx/(totalL-overlapx+0.00001)
#     if overlapx > totalL:
#         overlap=overlapx/totalL
    
#     return overlap
# print('                    ')

# def AllOverlapTests():
#     print("Entering OverlapTests")
    
#     df=pd.read_csv('/home/akh/myprojects/Linking-and-Clustering-Dikes/dikedata/test_overlap.csv')
#     df=DikesetReProcess(df)
#     fig,ax=plt.subplots()
    
#     plotlines(df, 'k', ax, ColorBy='label')
    
#     # for i in df['label'].unique():
#     #     lines=df[ df['label']==i]
#     #     o=overlapSegments(lines)
        
#     #     if ~(np.isclose(o, lines['trueoverlap'].iloc[0])):
            
#     #         print('Test', str(i),' Failed')
#     #         print("Overlap calcuation error!")
#     #         print(o, lines['trueoverlap'].iloc[0])
#     #         print('                    ')
#     #     else: 
#     #         print('Test', str(i), 'passed')
#     #         print('                    ')
            
#     for i in df['label'].unique():
#         lines=df[ df['label']==i]
#         o=RotateOverlapTest(lines)
    
#         if ~(np.isclose(o, lines['trueoverlap'].iloc[0])):
#             print('Test', str(i),' Failed')
#             print("Overlap calcuation error!")
#             print(o, lines['trueoverlap'].iloc[0])
#             print('                    ')
            
    
    
#         else: 
#             print('Test', str(i), 'passed')
#             print('                    ')
            
            

# #test radial fit 
# def TestRadialFit():

#     xdata=np.linspace(-90,90,100)
#     ydata0=np.cos(np.deg2rad(xdata))+np.sin(np.deg2rad(xdata))
#     ydata=100*np.cos(np.deg2rad(xdata))-900*np.sin(np.deg2rad(xdata))
#     xc=0
#     yc=0
    
#     test0=fromHT(xdata,ydata0, test=True)
#     test1=fromHT(xdata,ydata, test=True)
    
    
#     for i,C in zip([test0,test1], ['0,0', '100,-900']):
#         Centers=RadialFit(i, plot=True)
        
#         print( "True Value", C)
#         print("Calulated", Centers['Center'])
    


# """ fit rec """
# def testRec(lines, xc, yc): 
#     fig, a=plt.subplots()
    
#     xi,yi=endpoints2(lines)
#     x0=xc
#     y0=xc
#     ang=-np.deg2rad(np.mean(lines['theta']))
#     xp, yp= rotateXYShift(ang, xi,yi, x0,y0)
#     #plotlines(lines, 'k', a)
    
#     for i in range(0,len(lines)):
#         a.plot( [xp[i], xp[i+len(lines)]],  [yp[i], yp[i+len(lines)]], 'r')
#     width=np.ptp(xp)
#     length=np.ptp(yp)
    
#     # if width>length :
#     #     length=width
#     #     width=length
#     xc=(max(xp)-min(xp))/2 + min(xp)
#     yc=(max(yp)-min(yp))/2 + min(yp)
    
#     xr=xc+width/2
#     xl=xc-width/2
#     yu=yc+length/2
#     yd=yc-length/2
#     xs=np.append(xr,xl)
#     ys=np.append(yu,yd)
    
    
#     xpi, ypi=unrotate(ang, xp, yp, x0, y0)
#     # for i in range(0,len(lines)):
#     #     a.plot( [xpi[i], xpi[i+9]],  [ypi[i], ypi[i+9]], 'b')
        
#     Xedges=np.array([xs[0], xs[0], xs[1], xs[1], xs[0]])
#     Yedges=np.array([ys[1], ys[0], ys[0], ys[1], ys[1]])
#     a.plot(Xedges, Yedges, 'r.-')
    
#     xs,ys=unrotate(ang,Xedges,Yedges,x0,y0)

#     a.plot(xs, ys, 'r.-')


# """Test Overlap
# """

# df=pd.read_csv('/home/akh/myprojects/Linking-and-Clustering-Dikes/dikedata/test_overlap.csv')
# df=DikesetReProcess(df)
# fig,ax=plt.subplots()

# plotlines(df, 'k', ax, ColorBy='label')

# # for i in df['label'].unique():
# #     lines=df[ df['label']==i]
# #     o=overlapSegments(lines)
    
# #     if ~(np.isclose(o, lines['trueoverlap'].iloc[0])):
        
# #         print('Test', str(i),' Failed')
# #         print("Overlap calcuation error!")
# #         print(o, lines['trueoverlap'].iloc[0])
# #         print('                    ')
# #     else: 
# #         print('Test', str(i), 'passed')
# #         print('                    ')
        
# def RotateOverlapTest(lines):
#     theta=np.mean(lines['theta'].values)
    
#     dfRotated=rotateData2(lines, (90-theta))
#     dfRotated=transformXstart(dfRotated)
#     fig,ax=plt.subplots()

#     plotlines(lines, 'k', ax)
#     plotlines(dfRotated, 'r', ax)
    
#     Xstart=dfRotated['Xstart'].to_numpy()
#     Ystart=dfRotated['Ystart'].to_numpy()
#     Xend=dfRotated['Xend'].to_numpy()
#     Yend=dfRotated['Yend'].to_numpy()
#     step=1
#     totalL=np.sum( np.sqrt( (Xstart-Xend)**2 + (Ystart-Yend)**2))
#     xs=[np.arange(min(x,y), max(x,y), step) for x,y in zip(np.floor(Xstart),np.ceil(Xend))]
    
#     l=np.max([len(xi) for xi in xs])
#     xs_sameLength=[np.append(xi,[np.nan]*(l-len(xi))) for xi in xs]

#     arr=np.vstack(xs_sameLength) #better
#     u,xcounts=np.unique(arr[~np.isnan(arr)], return_counts=True)
    
#     overlapx=np.float64(np.sum(xcounts[xcounts>1]-1)*step)
    
#     print(overlapx)
#     print(totalL)
#     overlap=overlapx/(totalL-overlapx+0.00001)
#     if overlapx > totalL:
#         overlap=overlapx/totalL
    
#     return overlap

# print('                    ')

# for i in df['label'].unique():
#     lines=df[ df['label']==i]
#     o=RotateOverlapTest(lines)

#     if ~(np.isclose(o, lines['trueoverlap'].iloc[0])):
#         print('Test', str(i),' Failed')
#         print("Overlap calcuation error!")
#         print(o, lines['trueoverlap'].iloc[0])
#         print('                    ')
        


#     else: 
#         print('Test', str(i), 'passed')
#         print('                    ')
        
        

