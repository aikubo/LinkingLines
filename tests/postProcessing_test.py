# #!/usr/bin/env python3
# # -*- coding: utf-8 -*-
# """
import pytest
import pandas as pd
import numpy as np
from linkinglines import dilation, rotateXYShift, DikesetReProcess, RadialFit
from linkinglines import fromHT, makeLinearDataFrame


class TestRotate:
    
    @staticmethod
    def test_rotate90():
        # Test case 1: Rotate 90 degrees counterclockwise about the origin
        x, y = 1, 0
        ang = np.pi / 2
        xp, yp = rotateXYShift(ang, x, y, 0, 0)
        assert np.isclose(xp, 0, atol=1e7)
        assert np.isclose(yp, 1)
        
    @staticmethod
    def test_rotate45():
        # Test case 2: Rotate 45 degrees clockwise about (1, 1)
        x, y = 2, 1
        ang = -np.pi / 4
        h, k = 1, 1
        xp, yp = rotateXYShift(ang, x, y, h, k)
        assert np.isclose(xp,0.7071067811865475)
        assert np.isclose(yp,-0.7071067811865475)
    
    @staticmethod
    def test_arrayRotate():
        #test that an array put in is an array output
        x=np.linspace(0,100)
        y=np.linspace(0,100)
       
        ang = np.pi / 2
        xp, yp = rotateXYShift(ang, x, y, 0, 0)
        
        assert len(x)==len(xp)
        assert len(y)==len(yp)
        
        

class TestDilation():
    #Create some fixtures of linear line dataframes
    @pytest.fixture(autouse=True)
    def df2(self):
        return makeLinearDataFrame(1000, 90, 0, 5000, 500, ndikes=5)

    @pytest.fixture(autouse=True)
    def df3(self):
        return makeLinearDataFrame(1000, 45, 0, 5000, 500, ndikes=5)
        
    
    def test_dilation_90_degrees(self, df2):
        EWDilation, NSDilation, EWbin, NSbin = dilation(df2, binWidth=1)
        
        assert np.isclose(np.sum(np.abs(EWDilation)), 0.0)
        assert np.isclose(np.max(NSDilation),4.0)
        assert len(EWDilation)==len(NSbin)
        assert len(NSDilation)==len(EWbin)
        

    def test_dilation_45_degrees(self, df3):
        EW45 = np.cos(np.deg2rad(45)) * 2
        NS45 = np.sin(np.deg2rad(45)) * 4
        
        EWDilation, NSDilation, EWbin, NSbin = dilation(df3, binWidth=1)
        
        print(np.max(EWDilation), EW45)
        print((np.max(NSDilation), NS45))
        
        assert np.isclose(np.max(EWDilation), EW45)
        assert np.isclose(np.max(NSDilation), NS45)
        assert len(EWDilation)==len(NSbin)
        assert len(NSDilation)==len(EWbin)
        
    
    def test_dilation_non_1_average_width(self, df2):
        EWDilation2, NSDilation2, _, _ = dilation(df2, binWidth=1, averageWidth=10)
        Value = [np.max(NSDilation2), np.max(EWDilation2)]
        Expected = [30, 0]
        
        assert np.allclose(Value, Expected)
        


class TestOverlap():
    
    def testfileOverlap(self):
        data = {
                'Xstart': [0, 7, 3, 3, 3, 2, 9, -3, -4, -4, -1, -1, -1],
                'Ystart': [5, 8, 2, 0, 2, 3, 9, 3, 3, 3, 2, 3, 4],
                'Xend': [5, 10, 6, 6, 10, 5, 12, -3, -4, -9, 5, 5, 5],
                'Yend': [1, 11, 5, 0, 2, 9, 15, 10, 6, 9, 2, 3, 4],
                'label': [1, 1, 1, 2, 2, 3, 3, 4, 4, 5, 6, 6, 6],
                'trueoverlap': [0.25234649409924387,0.25234649409924387,0.25234649409924387,
                                0.66666, 0.66666, 
                                0, 0, 
                                0.666666, 0.6666666, 
                                0, 
                                2, 2, 2]
            }

        df = pd.DataFrame(data)
        df=DikesetReProcess(df, xc=0, yc=0)
 
        for i in df['label'].unique():
            lines=df[ df['label']==i]
            o=RotateOverlap(lines)
            print(o[0] , df[ df['label']==i]['trueoverlap'].values[0])
            assert np.isclose(o[0],df[ df['label']==i]['trueoverlap'].values[0])

class TestFitRadial():
    
    @pytest.fixture(autouse=True)
    def df1(self):
        xdata=np.linspace(-90,90,100)
        ydata0=np.cos(np.deg2rad(xdata))+np.sin(np.deg2rad(xdata))

        test0=fromHT(xdata,ydata0)
        return test0

    @pytest.fixture(autouse=True)
    def df2(self):
        xdata=np.linspace(-90,90,100)
        ydata0=100*np.cos(np.deg2rad(xdata))-900*np.sin(np.deg2rad(xdata))

        test0=fromHT(xdata,ydata0)
        return test0


    def test_outputs(self, df1):
        Centers=RadialFit(df1)
        assert isinstance(Centers, pd.DataFrame)
        assert len(Centers["Center"][0])==2
        
    
    def test_tworadial(self, df1, df2):
        
        true=np.array([[0,0], [100,-900]])
        
        for i,C in zip([df1, df2], true):
            Centers=RadialFit(i, plot=False)
            
            print( "True Value", C)
            print("Calulated", Centers['Center'])
            assert np.allclose(Centers["Center"][0], C, atol=1)



