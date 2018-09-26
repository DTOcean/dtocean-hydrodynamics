# -*- coding: utf-8 -*-
"""
Created on Fri Jun 02 12:37:42 2017

@author: mtopper
"""

import numpy as np

from dtocean_wave.utils.spec_class import wave_spec


def test_wave_spec_Directional():
    
    freqs = 1. / np.linspace(1, 10, 10)
    t = np.linspace(0, 2 * np.pi, 4, endpoint=False)
    
    test = wave_spec(freqs, 0.5, 5, s=25, t=t)
    _, spec = test.Regular()
    result = test.Directional(spec)
    
    assert np.isclose(result[:,2], 0.).all()
