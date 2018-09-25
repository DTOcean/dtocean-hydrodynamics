# -*- coding: utf-8 -*-
"""
Created on Tue Sep 25 16:28:33 2018

@author: Mathew Topper
"""


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
