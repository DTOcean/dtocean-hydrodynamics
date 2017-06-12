# -*- coding: utf-8 -*-
"""
Created on Fri Jun 02 12:37:42 2017

@author: mtopper
"""


import numpy as np

from dtocean_tidal.utils.distance_from_streamline import (
                                                distance_from_streamline)


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
