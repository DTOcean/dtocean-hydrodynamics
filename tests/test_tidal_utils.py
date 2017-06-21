# -*- coding: utf-8 -*-
"""
Created on Fri Jun 02 12:37:42 2017

@author: mtopper
"""

import pytest

import numpy as np

from dtocean_tidal.utils.distance_from_streamline import (
                                                distance_from_streamline)
from dtocean_tidal.utils.interpolation import interp_at_point


def test_distance_from_streamline():
    
    streamline = np.array([[0, 5, 10, 11],
                           [0, 0, 0, -5]])
    
    turbine_positions = np.array([[-5,  0],
                                  [ 9, -3],
                                  [12, -10]])
    
    turbine_distances = distance_from_streamline(streamline,
                                                 turbine_positions)
    
    assert len(turbine_distances) == 1
    
    dist_array = turbine_distances.values()[0]
        
    assert dist_array[0] > 10.
    assert dist_array[1] < 2.


@pytest.mark.parametrize("test_x, test_y, expected", [
    (4.5, 1, 0),
    (3.5, 1, -0.5),
    (3, 1, -1),
    (2.5, 1, -0.5),
    (1.5, 1, 0.5),
    (1, 1, 1),
    (0.5, 1, 0.5),
    (-0.5, 1, 0),
    (1.5, -1, 0.5),
    (1.5, 3, 0.5),
    ])
def test_interp_at_point(test_x, test_y, expected):
    
    X = np.array([0, 1, 2, 3, 4])
    Y = np.array([0, 1, 2])
    
    U = np.array([[np.nan, 1, 0, -1, np.nan],
                  [np.nan, 1, 0, -1, np.nan],
                  [np.nan, 1, 0, -1, np.nan]])
    Q = [U]
    
    result = interp_at_point(test_x, test_y, X, Y, Q)
    
    assert np.isclose(result[0], expected)
