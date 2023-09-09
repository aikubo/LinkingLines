#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Sep  9 10:01:18 2023

@author: akh
"""

from htMOD import * 
from PrePostProcess import * 
from synthetic import * 

import pytest
import pandas as pd
import numpy as np


class TestHough:
    
    
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
        theta, rho, result_xc, result_yc = AKH_HT(data, xc=xc, yc=yc)
        
        # Perform assertions
        assert result_xc == xc
        assert result_yc == yc
        
        theta_true=np.array([-45., -45., -45., -45.])
        rho_true=np.array([0,0,0,0])
        assert len(theta) == len(data)
        assert all([np.isclose(a,b) for a, b in zip(theta, theta_true)])
        
        assert len(rho) == len(data)
        assert all([np.isclose(a,b) for a, b in zip(rho, rho_true)])

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
        theta, rho, xc, yc = AKH_HT(data)
        
        # Perform assertions
        # Add assertions to check the calculated center (xc, yc), theta, and rho
        assert isinstance(theta, np.ndarray)
        assert isinstance(rho, np.ndarray)
        theta_true=np.array([-45., -45., -45., -45.])
        rho_true=np.array([0,0,0,0])
        assert len(theta) == len(data)
        assert all([np.isclose(a,b) for a, b in zip(theta, theta_true)])
        
        assert len(rho) == len(data)
        assert all([np.isclose(a,b) for a, b in zip(rho, rho_true)])
        
        assert xc==1.5
        assert yc==1.5
        
        assert isinstance(xc, float)
        assert isinstance(yc, float)
        
    @staticmethod
    def test_emptyDataFrame_akhht():
       # Create an empty DataFrame with the required columns
        empty_df = pd.DataFrame(columns=['Xstart', 'Ystart', 'Xend', 'Yend'])
    
        # Call the AKH_HT function with the empty DataFrame
        theta, rho, xc, yc = AKH_HT(empty_df)
    
        # Add assertions to check the expected behavior with an empty DataFrame
        assert theta is None  # Theta should be None or a default value
        assert rho is None    # Rho should be None or a default value
        assert xc is None     # xc should be None or a default value
        assert yc is None     # yc should be None or a default value
        
        
class testRotate:
    
    def test_rotate_xy_shift():
        # Test case 1: Rotate 90 degrees counterclockwise about the origin
        x, y = 1, 0
        ang = np.pi / 2
        xp, yp = rotateXYShift(ang, x, y, 0, 0)
        assert np.isclose(xp, 0)
        assert np.isclose(yp, 1)
    
        # Test case 2: Rotate 45 degrees clockwise about (1, 1)
        x, y = 2, 1
        ang = -np.pi / 4
        h, k = 1, 1
        xp, yp = rotateXYShift(ang, x, y, h, k)
        assert np.isclose(xp, 1 + np.sqrt(2))
        assert np.isclose(yp, 1 - np.sqrt(2))

