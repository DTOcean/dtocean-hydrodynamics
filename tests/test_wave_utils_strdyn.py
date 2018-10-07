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

from dtocean_wave.utils.StrDyn import EnergyProduction


def test_EnergyProduction():
    
    NBo = 1
    B = np.array([0, 180.]) / 180. * np.pi
    Hs = np.array([0.5, 1, 1.5, 2.])
    Tp = np.array([3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5])
    wdir = np.linspace(0, 360, 30, endpoint=False) / 180. * np.pi
    period = np.linspace(2, 15, 50)
    ScatDiag = (np.ones((7, 4, 2)) / (7 * 4 * 2), ('Jonswap', 3.3, 0))
    M = np.eye(3)
    Madd = np.array([np.eye(3)] * 50)
    Cpto = np.array([[[np.eye(3)] * 2] * 4] * 7)
    Crad = np.array([np.eye(3)] * 50)
    Khyd = np.eye(3)
    Fex = np.ones((50, 30, 3))
    Kfit = np.array([[[np.eye(3)] * 2] * 4] * 7)
    Cfit = np.array([[[np.eye(3)] * 2] * 4] * 7)
    RatedPower = 0.05

    Pyr, P_dev = EnergyProduction(NBo,
                                  B,
                                  Hs,
                                  Tp,
                                  wdir,
                                  period,
                                  ScatDiag,
                                  M,
                                  Madd,
                                  Cpto,
                                  Crad,
                                  Cpto,
                                  Khyd,
                                  Fex,
                                  Kfit,
                                  Cfit,
                                  RatedPower)

    assert Pyr.shape == (1L, 7L, 4L, 2L)
    assert P_dev.shape == (1L, 7L, 4L, 2L)
    assert P_dev.sum() > Pyr.sum()
