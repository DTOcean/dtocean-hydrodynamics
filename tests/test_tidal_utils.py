# -*- coding: utf-8 -*-

#    Copyright (C) 2017-2018 Mathew Topper
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.

"""
Created on Fri Jun 02 12:37:42 2017

.. moduleauthor:: Mathew Topper <mathew.topper@dataonlygreater.com>
"""

import pytest

import numpy as np

from dtocean_tidal.utils.distance_from_streamline import (
                                                distance_from_streamline)
from dtocean_tidal.utils.interpolation import interp_at_point
from dtocean_tidal.utils.misc import pi2pi


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


@pytest.mark.parametrize("A, B", [
    (3, 1),
    (2, 0),
    (3. / 2, -1. / 2),
    (1, 1),
    (1. / 2, 1. / 2),
    (0, 0),
    (-1. / 2, -1. / 2),
    (-1, 1),
    (-3. / 2, 1. / 2),
    (-2, 0),
    (-3, 1)
    ])
def test_pi2pi(A, B):
    
    result = pi2pi(A * np.pi)
    
    assert np.isclose(result, B * np.pi)
