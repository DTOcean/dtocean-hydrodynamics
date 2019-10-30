# -*- coding: utf-8 -*-
"""
Created on Thu Oct 30 13:06:34 2019
"""

import numpy as np
import numpy.testing as npt
import pytest

from dtocean_tidal.modules.array_yield import _get_turbine_power


@pytest.mark.parametrize("diameter, "
                         "cut_in, "
                         "cut_out, "
                         "rating, "
                         "angle_of_attack, "
                         "Cp, "
                         "u, "
                         "expected", [
    (20., 1, 3.2, 1.5e6, 0., 0.3, 3, 1304153.65),
    (20., 1, 3.2, 1.5e6, 30., 0.3, 3, 1129430.19),
    (20., 1, 3.2, 1.5e6, -30., 0.3, 3, 1129430.19),
    (20., 1, 3.2, 1.5e6, 0., 0.3, 3.15, 1.5e6),
    (20., 1, 3.2, 1.5e6, 0., 0.3, 0.9, 0.),
    (20., 1, 3.2, 1.5e6, 0., 0.3, 3.3, 0.),
    (20., 1, 3.2, 1.5e6, 0., -0.3, 3, 0.)
    ])
def test_get_turbine_power(diameter,
                           cut_in,
                           cut_out,
                           rating,
                           angle_of_attack,
                           Cp,
                           u,
                           expected):
    
    test = _get_turbine_power(diameter,
                              cut_in,
                              cut_out,
                              rating,
                              angle_of_attack,
                              Cp,
                              u)
    
    if expected == 0.:
        npt.assert_almost_equal(test, expected, decimal=10)
    else:
        assert np.isclose(test, expected)
