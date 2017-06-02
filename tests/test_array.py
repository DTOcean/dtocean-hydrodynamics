# -*- coding: utf-8 -*-
"""
Created on Fri Jun 02 12:37:42 2017

@author: mtopper
"""

import pytest

import numpy as np
from shapely.geometry import Polygon

from dtocean_hydro.array import Array_pkg


def test_Array_pkg_init():
    
    lease = np.array([[-100,-100],
                      [100,-100],
                      [100,100],
                      [-100,100]],'f')
    
    Array_pkg(lease,
              lease*1.05,
              (10, 10),
              45.0/180*np.pi,
              Polygon(lease*0.2),
              True)
    
    assert True


@pytest.mark.parametrize("test_coords, expected",
    [(np.array([[0., 0.]]), False),
     (np.array([[0., 0.],[0., 11.]]), False),
     (np.array([[0., 0.],[11., 0.]]), False),
     (np.array([[0., 0.],[11., 11.]]), False),
     (np.array([[0., 0.],[0., 1.]]), True),
     (np.array([[0., 0.],[1., 0.]]), True),
     (np.array([[0., 0.],[1., 1.]]), True),
    ])
def test_checkMinDist(test_coords, expected):
    
    lease = np.array([[-100,-100],
                      [100,-100],
                      [100,100],
                      [-100,100]],'f')
    
    arr = Array_pkg(lease,
                    lease*1.05,
                    (10, 10),
                    45.0/180*np.pi,
                    Polygon(lease*0.2),
                    True)
    
    arr.coord = test_coords
    arr.checkMinDist()
    
    assert arr.minDist_constraint == expected
    

def test_checkMinDist_error():
    
    lease = np.array([[-100,-100],
                      [100,-100],
                      [100,100],
                      [-100,100]],'f')
    
    arr = Array_pkg(lease,
                    lease*1.05,
                    (10, 10),
                    45.0/180*np.pi,
                    Polygon(lease*0.2),
                    True)
    
    with pytest.raises(IOError):
        arr.checkMinDist()
