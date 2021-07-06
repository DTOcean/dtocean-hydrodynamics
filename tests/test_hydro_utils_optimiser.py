# -*- coding: utf-8 -*-

#    Copyright (C) 2017-2021 Mathew Topper
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
.. moduleauthor:: Mathew Topper <mathew.topper@dataonlygreater.com>
"""

import pytest

from copy import deepcopy

from dtocean_hydro.input import WP2input, WP2_MachineData, WP2_SiteData
from dtocean_hydro.main import WP2


@pytest.fixture
def searchoptimum(tidalsite, tidal, tidal_kwargs):
    
    tidal = deepcopy(tidal)
    tidal[-3]['Option'] = 1
    tidal[-3]['Value'] = "rectangular"
    
    site = WP2_SiteData(*tidalsite)
    machine = WP2_MachineData(*tidal, **tidal_kwargs)
    
    data = WP2input(machine, site)
    wp2 = WP2(data)
    hyd_obj = wp2._get_hyd_obj()
    
    return wp2._get_optim_obj(hyd_obj)


def test_SearchOptimum_estimate_start_point(searchoptimum):
    # This isn't a good test (or good code for testing)
    test = searchoptimum.estimate_start_point()
    assert len(test) == 2
    assert len(test[0]) == 2
