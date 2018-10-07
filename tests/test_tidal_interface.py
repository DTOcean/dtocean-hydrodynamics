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
