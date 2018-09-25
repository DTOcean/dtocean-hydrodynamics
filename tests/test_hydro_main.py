
import pytest
import numpy as np

from dtocean_hydro.input import WP2input, WP2_MachineData, WP2_SiteData
from dtocean_hydro.main import get_device_depths, WP2

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

