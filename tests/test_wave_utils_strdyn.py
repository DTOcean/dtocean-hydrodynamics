# -*- coding: utf-8 -*-
"""
Created on Fri Jun 02 12:37:42 2017

@author: mtopper
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
    RatedPower = 1e6

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
