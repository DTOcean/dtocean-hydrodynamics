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


def test_generator():
    
    lease = np.array([[-100,-100],
                      [100,-100],
                      [100,100],
                      [-100,100]],'f')
    
    arr = Array_pkg(lease,
                    lease*1.05,
                    (10, 10),
                    np.pi / 4)
    
    arr.generator(2, 2, 100, 100, np.pi / 2, 0)
    
    p1 = arr.coord[0, :]
    p2 = arr.coord[1, :]
    dist = np.sqrt((p2[0] - p1[0]) ** 2 + (p2[1] - p1[1]) ** 2)
    
    assert len(arr.coord) == 4
    assert np.isclose(dist, 100)


def test_show(mocker):
    
    lease = np.array([[-100,-100],
                      [100,-100],
                      [100,100],
                      [-100,100]],'f')
    
    arr = Array_pkg(lease,
                    lease*1.05,
                    (10, 10),
                    np.pi / 4)
    
    arr.generator(2, 2, 100, 100, np.pi / 2, 0)
    
    show_patch = mocker.patch("dtocean_hydro.array.plt.show")
    arr.show()
    
    assert show_patch.called
