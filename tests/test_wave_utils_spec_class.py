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

from dtocean_wave.utils.spec_class import wave_spec


def test_wave_spec_Directional():
    
    freqs = 1. / np.linspace(1, 10, 10)
    t = np.linspace(0, 2 * np.pi, 4, endpoint=False)
    
    test = wave_spec(freqs, 0.5, 5, s=25, t=t)
    _, spec = test.Regular()
    result = test.Directional(spec)
    
    assert np.isclose(result[:,2], 0.).all()
