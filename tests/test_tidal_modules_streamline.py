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
Created on Wed Sep 26 13:17:52 2018

.. moduleauthor:: Mathew Topper <mathew.topper@dataonlygreater.com>
"""

import pytest
import numpy as np

from dtocean_tidal.main import Hydro
from dtocean_tidal.modules.streamline import Streamlines

YMAX = 350.0
BATHYGLOB = -31.5


@pytest.fixture
def data():

    ## Velocity field and site data
    xmax = 840.0
    ymax = YMAX
    lease = np.asarray([[0.0 , 0.0 ],
                        [0.0 , ymax],
                        [xmax, ymax],
                        [xmax, 0.0 ]])
    
    x = np.linspace(0.0, xmax, (xmax/10)+1) # dx = 10 m
    y = np.linspace(0.0, ymax, (ymax/10)+1) # dy = 10 m
    X, _ = np.meshgrid(x,y)
    BR = 1.0  # blockage ratio
    
    umax = 3.69 # maximum velocity in the X direction
    vmax = 0.0 # maximum velocity in the Y direction
    sshmax = 0.0 # maximum sea surface elevation
    timax= 0.1
    bathy = BATHYGLOB
    
    U = np.ones(X.shape) * umax
    V = np.ones(X.shape) * vmax
    SSH = np.ones(X.shape) * sshmax
    TI = np.ones(X.shape) * timax
    BATHY = np.ones(X.shape) * bathy
    beta = 0.4
    alpha = 7.
    
    data = {}
    data['TI'] = TI
    data['X'] = x  # save only 1D array as structured grid assumed
    data['Y'] = y  # save only 1D array as structured grid assumed
    data['U'] = U
    data['V'] = V
    data['SSH'] = SSH
    data['bathy'] = BATHY
    data['BR'] = BR
    data['lease'] = lease
    data['beta'] = beta
    data['alpha'] = alpha
    
    hydro = Hydro(data)
    
    data = {}
    xn = hydro.X.shape[0]
    yn = hydro.Y.shape[0]
    
    xmin = hydro.X.min()
    xmax = hydro.X.max()
    ymin = hydro.Y.min()
    ymax = hydro.Y.max()
    dx = (xmax - xmin) / (xn - 1)
    dy = (ymax - ymin) / (yn - 1)
    xi = np.arange(xmin, xmax + dx, dx)
    yi = np.arange(ymin, ymax + dy, dy)
    
    data['X'] = xi
    data['Y'] = yi
    data['U'] = hydro.U
    data['V'] = hydro.V
    data['interpU'] = hydro.interpU
    data['interpV'] = hydro.interpV
    data['lease'] = hydro.lease
    
    return data


def test_Streamlines_plot(mocker, data):
    
    positions = np.array([[420., 175.,  20.],
                          [514.5, 233.33333333,  20.]])
    
    SLs = Streamlines(data,
                      positions,
                      2)
    
    show_patch = mocker.patch("dtocean_tidal.modules.streamline.plt.show")
    SLs.plot()
    
    assert show_patch.called


def test_Streamlines_detect_loop(mocker, data):
    
    positions = np.array([[420., 175.,  20.],
                          [514.5, 233.33333333,  20.]])
    
    SLs = Streamlines(data,
                      positions,
                      2)
    
    x = y = [0, 0]
    
    assert SLs._detectLoop(x, y)


def test_Streamlines_dont_detect_loop(mocker, data):
    
    positions = np.array([[420., 175.,  20.],
                          [514.5, 233.33333333,  20.]])
    
    SLs = Streamlines(data,
                      positions,
                      2)
    
    x = [0, 0]
    y = [0, 20]
    
    assert not SLs._detectLoop(x, y)
    

def test_Streamlines_makeStreamline(mocker, data):
    
    positions = np.array([[420., 175.,  20.],
                          [514.5, 233.33333333,  20.]])
    
    SLs = Streamlines(data,
                      positions,
                      2)
    
    sx, sy = SLs._makeStreamline(1, YMAX / 2)
    
    assert sx[-1] > sx[0]
    assert np.isclose(sy, 175).all()
