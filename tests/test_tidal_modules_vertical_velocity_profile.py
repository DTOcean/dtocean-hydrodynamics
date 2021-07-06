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
Created on Wed Sep 26 13:17:52 2018

.. moduleauthor:: Mathew Topper <mathew.topper@dataonlygreater.com>
"""

import numpy as np
from dtocean_tidal.modules.vertical_velocity_profile import vvpw

def test_vvpw():
    test = vvpw(16, 0, 20, 0.2, 2)
    assert np.isclose(test, 2)

def test_vvpw_beta_neg():
    test = vvpw(16, 0, 20, -1, 2)
    assert test == 0

def test_vvpw_alpha_neg():
    test = vvpw(16, 0, 20, 0.2, -1)
    assert test == 0

def test_vvpw_debug():
    test = vvpw(16, 0, 20, 0.2, 2, True)
    assert np.isclose(test, 2)
