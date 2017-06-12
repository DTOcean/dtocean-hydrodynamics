# -*- coding: utf-8 -*-
"""
Created on Fri Jun 02 12:37:42 2017

@author: mtopper
"""


import numpy as np
import numpy.ma as ma

from dtocean_tidal.interface import get_indices


def test_get_indices():

    base = np.array([3, 5, 7, 1, 9, 8, 6, 6])
    search = np.array([2, 1, 5, 10, 100, 6])
    
    expected = ma.array([0, 3L, 1L, 0, 0, 6L],
                        mask=[1, 0, 0, 1, 1, 0])
    
    result = get_indices(search, base)
    
    assert (result == expected).all()
    

