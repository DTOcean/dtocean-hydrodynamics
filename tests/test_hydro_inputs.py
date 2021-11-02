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
Created on Tue Sep 25 16:28:33 2018

.. moduleauthor:: Mathew Topper <mathew.topper@dataonlygreater.com>
"""

import numpy as np

from dtocean_hydro.input import WP2input, WP2_MachineData, WP2_SiteData


def test_WP2_MachineData_init(tidal):
    
    WP2_MachineData(*tidal)
    
    assert True
    
    
def test_WP2_SiteData_init(wavesite):
    
    WP2_SiteData(*wavesite)
    
    assert True


def test_WP2input_wave_init(wavesite, wave, wave_data_folder):
    
    site = WP2_SiteData(*wavesite)
    machine = WP2_MachineData(*wave,
                              wave_data_folder=wave_data_folder)
    
    WP2input(machine, site)
    
    assert True
    
    
def test_WP2input_wavebiggamma_init(wavesitebiggamma, wave, wave_data_folder):
    
    site = WP2_SiteData(*wavesitebiggamma)
    machine = WP2_MachineData(*wave,
                              wave_data_folder=wave_data_folder)
    
    WP2input(machine, site)
    
    assert True
    

def test_WP2input_tidal_init(tidalsite, tidal, tidal_kwargs):
    
    site = WP2_SiteData(*tidalsite)
    machine = WP2_MachineData(*tidal, **tidal_kwargs)
    
    WP2input(machine, site)
    
    assert True


def test_WP2input_tidal_mainAngle(tidalsite_simple, tidal, tidal_kwargs):
    
    site = WP2_SiteData(*tidalsite_simple)
    machine = WP2_MachineData(*tidal, **tidal_kwargs)
    test = WP2input(machine, site)
    
    assert np.isclose(test.S_data.mainAngle, np.pi/2)


def test_WP2input_tidal_MainDirection(tidalsite_simple, tidal, tidal_kwargs):
    
    tidalsite_simple[5] = np.array((-1, 0))
    
    site = WP2_SiteData(*tidalsite_simple)
    machine = WP2_MachineData(*tidal, **tidal_kwargs)
    test = WP2input(machine, site)
    
    assert np.isclose(test.S_data.mainAngle, np.pi)
