#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Sep  9 10:01:04 2023

@author: akh
"""

import pytest
import pandas as pd
import numpy as np
import linkinglines as ll

# Create a test class for AggCluster
class TestAggCluster:

    # Create synthetic data for testing
    @pytest.fixture
    def synthetic_dikeset(self, autouse=True):
        # Create synthetic data here 
        theta=np.array([-75, -30, 5, 45, 60, 90])
        rho = np.random.uniform(low=-10, high=10, size=6)
        synthetic_dikeset = ll.fromHT(theta, rho, scale=1)
        return synthetic_dikeset

    # Test cases for AggCluster function without rotation
    def test_AggCluster_without_rotation(self, synthetic_dikeset):
        
        df = ll.fragmentDikes(synthetic_dikeset)
        
        dikeset, Z = ll.AggCluster(df, dtheta=1.0, drho=1.0, rotate=False)

        # Add assertions to check the results
        assert isinstance(dikeset, pd.DataFrame)
        assert isinstance(Z, np.ndarray)
        assert len(dikeset) == len(df)

    # Test cases for AggCluster function with rotation
    def test_AggCluster_with_rotation(self, synthetic_dikeset):
        dikeset, Z = ll.AggCluster(synthetic_dikeset, dtheta=1.0, drho=1.0, rotate=True)

        # Add assertions to check the results after rotation
        assert isinstance(dikeset, pd.DataFrame)
        assert isinstance(Z, np.ndarray)
        assert len(dikeset) == len(synthetic_dikeset)
        # Add more specific assertions based on your expected behavior

    # Test cases for AggCluster function with custom parameters
    def test_AggCluster_with_custom_parameters(self, synthetic_dikeset):
        dikeset, Z = ll.AggCluster(synthetic_dikeset, dtheta=0.5, drho=0.5)

        # Add assertions to check the results with custom parameters
        assert isinstance(dikeset, pd.DataFrame)
        assert isinstance(Z, np.ndarray)
        assert len(dikeset) == len(synthetic_dikeset)


