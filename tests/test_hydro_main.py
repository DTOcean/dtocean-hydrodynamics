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

from copy import deepcopy

import pytest
import numpy as np

from dtocean_hydro.input import WP2input, WP2_MachineData, WP2_SiteData
from dtocean_hydro.main import get_device_depths, WP2
from dtocean_hydro.utils.optimiser import SearchOptimum


class MockSearch(SearchOptimum):
    
    def eval_optimal_layout(self):
        
        xmap = [50, 25, 0, np.pi/2]
        
        if self._Opt == 1:
            NR, NC, IR, IC, beta, psi = self.param_conditioning(xmap)
            self._array.generator(NR, NC, IR, IC, beta, psi)
        else:
            self._array.coord = (self._Val*int(xmap[0]) /
                                                 self._normalisation_point)


def test_get_device_depths():

    bathy = np.array([[ 0,  0, -10],
                      [ 0, 10, -20],
                      [ 0, 20, -30],
                      [10,  0, -40],
                      [10, 10, -50],
                      [10, 20, -60],
                      [20,  0, -70],
                      [20, 10, -80],
                      [20, 20, -90]])
    
    layout = np.array([[ 0,   10],
                       [17.5, 17.5]])
    
    test_depths = get_device_depths(bathy, layout)
    
    assert np.isclose(test_depths[0], -20)
    assert np.isclose(test_depths[1], -90)


def test_get_device_depths_bad_layout1():
    
    bathy = np.array([[ 0,  0, -10],
                      [ 0, 10, -20]])
    
    layout = np.array([ 0,   10])
    
    with pytest.raises(IndexError):
        get_device_depths(bathy, layout)


def test_get_device_depths_bad_layout2():
    
    bathy = np.array([[ 0,  0, -10],
                      [ 0, 10, -20]])
    
    layout = np.array([[   0,   10, -10.],
                       [17.5, 17.5, -10.]])
    
    with pytest.raises(IndexError):
        get_device_depths(bathy, layout)


def test_WP2_init_wave(wavesite, wave, wave_data_folder):
    
    site = WP2_SiteData(*wavesite)
    machine = WP2_MachineData(*wave,
                              wave_data_folder=wave_data_folder)
    
    test = WP2input(machine, site)
    WP2(test)
    
    assert True
    
    
def test_WP2_init_tidal(tidalsite, tidal, tidal_kwargs):
    
    site = WP2_SiteData(*tidalsite)
    machine = WP2_MachineData(*tidal, **tidal_kwargs)
    
    test = WP2input(machine, site)
    WP2(test)
    
    assert True
    

def test_WP2_optimisationLoop_wave(wavesite, wave, wave_data_folder):
    
    site = WP2_SiteData(*wavesite)
    machine = WP2_MachineData(*wave,
                              wave_data_folder=wave_data_folder)
    
    data = WP2input(machine, site)
    test = WP2(data)
    result  = test.optimisationLoop()
    
    assert result


def test_WP2_optimisationLoop_tidal(tidalsite, tidal, tidal_kwargs):
    
    site = WP2_SiteData(*tidalsite)
    machine = WP2_MachineData(*tidal, **tidal_kwargs)
    
    data = WP2input(machine, site)
    test = WP2(data)
    result  = test.optimisationLoop()
    
    assert result


def test_WP2_optimisationLoop_one_outside(tidalsite, tidal, tidal_kwargs):
    
    # Append position
    tidal = deepcopy(tidal)
    tidal[-3]['Value'] = np.append(tidal[-3]['Value'], [[1100., 400.]],
                                   axis=0)
    
    site = WP2_SiteData(*tidalsite)
    machine = WP2_MachineData(*tidal, **tidal_kwargs)
    
    data = WP2input(machine, site)
    test = WP2(data)
    
    with pytest.raises(RuntimeError) as excinfo:
        test.optimisationLoop()
    
    assert "have been excluded." in str(excinfo.value)


def test_WP2_optimisationLoop_all_outside(tidalsite, tidal, tidal_kwargs):
    
    # Change all positions
    tidal = deepcopy(tidal)
    tidal[-3]['Value'] = np.array([[1100., 400.]])
    
    site = WP2_SiteData(*tidalsite)
    machine = WP2_MachineData(*tidal, **tidal_kwargs)
    
    data = WP2input(machine, site)
    test = WP2(data)
    
    assert test.optimisationLoop() == -1


def test_WP2_optimisationLoop_cma_rectangular(tidalsite, tidal, tidal_kwargs):
    
    # Change all positions
    tidal = deepcopy(tidal)
    tidal[-3]['Option'] = 1
    tidal[-3]['Value'] = "rectangular"
    
    site = WP2_SiteData(*tidalsite)
    machine = WP2_MachineData(*tidal, **tidal_kwargs)
    
    data = WP2input(machine, site)
    test = WP2(data, search_class=MockSearch)
    
    assert test.optimisationLoop()


def test_WP2_optimisationLoop_monte_staggered(tidalsite, tidal, tidal_kwargs):
    
    # Change all positions
    tidal = deepcopy(tidal)
    tidal[-3]['Option'] = 1
    tidal[-3]['Value'] = 'staggered'
    
    site = WP2_SiteData(*tidalsite)
    machine = WP2_MachineData(*tidal, **tidal_kwargs)
    
    data = WP2input(machine, site)
    test = WP2(data, search_class=MockSearch, optim_method=2)
    
    assert test.optimisationLoop()


def test_WP2_optimisationLoop_brute_full(tidalsite, tidal, tidal_kwargs):
    
    # Change all positions
    tidal = deepcopy(tidal)
    tidal[-3]['Option'] = 1
    tidal[-3]['Value'] = 'full'
    
    site = WP2_SiteData(*tidalsite)
    machine = WP2_MachineData(*tidal, **tidal_kwargs)
    
    data = WP2input(machine, site)
    test = WP2(data, search_class=MockSearch, optim_method=3)
    
    assert test.optimisationLoop()
