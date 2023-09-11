#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Sep 9 10:01:18 2023

@author: akh
"""

# Import necessary modules and libraries
from htMOD import HoughTransform, CyclicAngleDist, rotateData, HT_center
from syntheticMod import fromHT
from PrePostProcess import * 
from syntheticMod import * 
import pytest
import pandas as pd
import numpy as np


class TestHTCenter:
    # Test case 1: Test with a DataFrame containing horizontal line segments
    def test_horizontal_line_segments(self):
        data = pd.DataFrame({
            'Xstart': [0, 0, 0],
            'Ystart': [0, 1, 2],
            'Xend': [3, 3 , 3],
            'Yend': [0, 1, 2]
        })
        
        xc, yc = HT_center(data)
        
        assert xc == 1.5
        assert yc == 1.0

    # Test case 2: Test with a DataFrame containing vertical line segments
    def test_vertical_line_segments(self):
        data = pd.DataFrame({
            'Xstart': [0, 1, 2],
            'Ystart': [0, 0, 0],
            'Xend': [0, 1, 2],
            'Yend': [3, 3, 3]
        })
        
        xc, yc = HT_center(data)
        
        assert xc == 1.0
        assert yc == 1.5

    # Test case 3: Test with a DataFrame containing diagonal line segments
    def test_diagonal_line_segments(self):
        data = pd.DataFrame({
            'Xstart': [0, 1, 2],
            'Ystart': [0, 1, 2],
            'Xend': [2, 1, 0],
            'Yend': [2, 1, 0]
        })
        
        xc, yc = HT_center(data)
        
        assert xc == 1.0
        assert yc == 1.0


# Define a test class for Hough Transform functions
class TestHough:
    # Test case for an empty DataFrame input to HoughTransform function
    @staticmethod
    def test_emptyDataFrame_akhht():
        # Create an empty DataFrame with the required columns
        empty_df = pd.DataFrame(columns=['Xstart', 'Ystart', 'Xend', 'Yend'])
        
        # Call the HoughTransform function with the empty DataFrame
        with pytest.raises(ValueError):
            theta, rho, xc, yc = HoughTransform(empty_df)
    
    # Test case for a non-DataFrame input to HoughTransform function
    @staticmethod
    def test_notDataFrame_akhht():
        # Create an array (non-DataFrame)
        array = np.linspace(0, 100)
        
        # Call the HoughTransform function with the non-DataFrame input
        with pytest.raises(ValueError):
            theta, rho, xc, yc = HoughTransform(array)
            
    # Test case for points not lines
    @staticmethod
    def test_points_ht():
        data = pd.DataFrame({
            'Xstart': [0, 1, 2, 3],
            'Ystart': [0, 1, 2, 3],
            'Xend': [0, 1, 2, 3],
            'Yend': [0, 1, 2, 3]
        })

        # Call the HoughTransform function with the points input
        with pytest.raises(ValueError):
            theta, rho, xc, yc = HoughTransform(data)

    # Test case for HoughTransform with specified center coordinates
    @staticmethod
    def test_akht_with_center():
        # Create a sample dataframe with line segments
        data = pd.DataFrame({
            'Xstart': [0, 1, 2, 3],
            'Ystart': [0, 1, 2, 3],
            'Xend': [3, 2, 1, 0],
            'Yend': [3, 2, 1, 0]
        })
        
        # Define the expected center (xc, yc)
        xc, yc = 1.5, 1.5
        
        # Calculate the Hough Transform
        theta, rho, result_xc, result_yc = HoughTransform(data, xc=xc, yc=yc)
        
        # Perform assertions
        assert result_xc == xc
        assert result_yc == yc
        
        # Check calculated theta and rho values
        theta_true = np.array([-45., -45., -45., -45.])
        rho_true = np.array([0, 0, 0, 0])
        assert len(theta) == len(data)
        assert all([np.isclose(a, b) for a, b in zip(theta, theta_true)])
        assert len(rho) == len(data)
        assert all([np.isclose(a, b) for a, b in zip(rho, rho_true)])
    
    # Test case for HoughTransform without specifying center coordinates
    @staticmethod
    def test_akht_without_center():
        # Create a sample dataframe with line segments
        data = pd.DataFrame({
            'Xstart': [0, 1, 2, 3],
            'Ystart': [0, 1, 2, 3],
            'Xend': [3, 2, 1, 0],
            'Yend': [3, 2, 1, 0]
        })
        
        # Calculate the Hough Transform without specifying center
        theta, rho, xc, yc = HoughTransform(data)
        
        # Perform assertions
        assert isinstance(theta, np.ndarray)
        assert isinstance(rho, np.ndarray)
        theta_true = np.array([-45., -45., -45., -45.])
        rho_true = np.array([0, 0, 0, 0])
        assert len(theta) == len(data)
        assert all([np.isclose(a, b) for a, b in zip(theta, theta_true)])
        assert len(rho) == len(data)
        assert all([np.isclose(a, b) for a, b in zip(rho, rho_true)])
        assert xc == 1.5
        assert yc == 1.5
        assert isinstance(xc, float)
        assert isinstance(yc, float)
    
    # Test case for HoughTransform with vertical line segments
    @staticmethod
    def test_akht_vertical():
        # Create a sample dataframe with vertical line segments
        data = pd.DataFrame({
            'Xstart': [0, 1, 2, 3],
            'Ystart': [0, 0, 0, 0],
            'Xend': [0, 1 , 2, 3],
            'Yend': [3, 3, 3, 3]
        })
        
        # Calculate the Hough Transform with specified center (0, 0)
        theta, rho, xc, yc = HoughTransform(data, xc=0, yc=0)
        
        # Perform assertions
        assert isinstance(theta, np.ndarray)
        assert isinstance(rho, np.ndarray)
        theta_true = np.array([0., 0., 0., 0.])
        rho_true = np.array([0, 1, 2, 3])
        assert all([np.isclose(a, b) for a, b in zip(theta, theta_true)])
        assert all([np.isclose(a, b) for a, b in zip(rho, rho_true)])
        assert xc == 0
        assert yc == 0
    
    # Test case for HoughTransform with horizontal line segments
    @staticmethod
    def test_akht_horizontal():
        # Create a sample dataframe with horizontal line segments
        data = pd.DataFrame({
            'Xstart': [0, 0, 0, 0],
            'Ystart': [0, 1, 2, 3],
            'Xend': [3, 3 , 3, 3],
            'Yend': [0, 1, 2, 3]
        })
        
        # Calculate the Hough Transform with specified center (0, 0)
        theta, rho, xc, yc = HoughTransform(data, xc=0, yc=0)
        
        # Perform assertions
        assert isinstance(theta, np.ndarray)
        assert isinstance(rho, np.ndarray)
        theta_true = np.array([90., 90., 90., 90.])
        rho_true = np.array([0, 1, 2, 3])
        assert all([np.isclose(a, b) for a, b in zip(theta, theta_true)])
        assert all([np.isclose(a, b) for a, b in zip(rho, rho_true)])
        assert xc == 0
        assert yc == 0

# Define a test class for preprocessing functions
class TestPreProcessingFunctions():
    # Test case for CyclicAngleDist function
    @staticmethod
    def testCyclicAngleDist():
        #Test cases for Cyclic angle distance
        assert CyclicAngleDist([20.], [20.]) == 0
        assert CyclicAngleDist([90.], [-90.]) == 0
        assert CyclicAngleDist([0.], [90.]) == 90.0
        assert CyclicAngleDist([0.], [-90.]) == 90.0
        assert CyclicAngleDist([45.], [-45.]) == 90.0

    @staticmethod
    def test_rotateData():
        #Test rotation of vertical to horizontal 
        #Create dataframe
        data = pd.DataFrame({
            'Xstart': [0, 1, 2, 3],
            'Ystart': [0, 0, 0, 0],
            'Xend': [0, 1 , 2, 3],
            'Yend': [3, 3, 3, 3]
        })

        #Rotate by 90
        dataRotated=rotateData(data,90)
        
        
        assert len(data)==len(dataRotated)
        theta_true = np.array([-90., -90., -90., -90.])
        assert all([np.isclose(a, b) for a, b in zip(dataRotated['theta'].values, theta_true)])
        
        xc,yc=HT_center(data)
        
        xc_r, yc_r=HT_center(dataRotated)
        
        assert xc==xc_r
        assert yc==yc_r

class TestSynthetic: 
    
    @staticmethod
    def test_fromHT():
        
        # Test non equal arrays 
        theta_true = np.array([90., 90.])
        rho_true=np.array([1., 1., 1., 1.])
        l=5
        with pytest.raises(ValueError):
            df=fromHT(theta_true, rho_true, length=l, scale=1)
        
        # Test 90 degree angle
        theta_true = np.array([90., 90., 90., 90.])
        rho_true=np.array([1., 1., 1., 1.])
        l=5
        df=fromHT(theta_true, rho_true, length=l, scale=1)
        
        assert len(theta_true)==len(df)
        
        assert all([np.isclose(l,a) for a in df['seg_length'].values])
        assert all([np.isclose(a,b) for a,b in zip(df['theta'].values, theta_true)])
        assert all([np.isclose(a,b) for a,b in zip(df['rho'].values, rho_true)])
        
        # Test 0 degree angle 
        theta_true = np.array([0.,0., 0., 0.])
        rho_true=np.array([1., 1., 1., 1.])
        l=5
        df=fromHT(theta_true, rho_true, length=l, scale=1)
        
        assert len(theta_true)==len(df)
        
        assert all([np.isclose(l,a) for a in df['seg_length'].values])
        assert np.allclose(df['theta'].values, theta_true, atol=1e-7)
        assert np.allclose(df['rho'].values, rho_true, rtol=1e-2)
        
        
        