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

from dtocean_tidal.main import (Array,
                                Hydro,
                                _is_within_yaw,
                                _get_angle_of_attack)

YMAX = 350.0
BATHYGLOB = -31.5
DIAM = 18.9
PIBY4 = np.pi / 4.
PIBY6 = np.pi / 6.
PIBY12 = np.pi / 12.


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
    X, Y = np.meshgrid(x,y)
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
    
    return data


@pytest.fixture
def turbines():
    
    # Turbines positions
    ymax = YMAX
    z = BATHYGLOB / 2.0 # hub height/depth
    coords = {}
    diam = DIAM
    first_row = 420.0 # x position of first row
    
    coords['turbine0'] = {}
    coords['turbine0']['position'] = np.asarray((first_row, (ymax/2.0), z))
    coords['turbine1'] = {}
    coords['turbine1']['position'] = np.asarray((first_row + 5.0 * diam,
                                                (ymax/1.5), z))
    
    return coords


@pytest.fixture
def features(turbines):
    
    ## Turbines features
    cut_in = 0.0 # cut-in speed
    cut_out = 10.0 # cut-out speed
    # actual turbine features
    speed = np.arange(0.0, 10.0, 0.2)
    
    CT = np.ones(speed.shape) * 0.76
    Ct = [speed,CT] # thrust curve
    
    CP = np.ones(speed.shape) * 0.3
    Cp = [speed, CP] # Power curve
    feat = {}
    
    for key in turbines.keys():
        
        feat[key] = {}
        feat[key]['OA'] = 315.0  # orientation angle (deg.),
                                 # turbines face North-West
        feat[key]['HAS'] = 0.0  # heading angle span (deg.),
                                # max = 180 deg. = full yaw
        feat[key]['Cp'] = Cp[:]  # Power curve
        feat[key]['Ct'] = Ct[:]  # thrust curve
        feat[key]['Diam'] = DIAM  # Diam = rotor diameter (m)
        feat[key]['cutIO'] = np.array([cut_in, cut_out])
        feat[key]['floating'] = False  # so hub height will be considered from
                                       # the seabottom upwards
        feat[key]['2way'] = False  # turbines work in one direction
    
    return feat


def test_Hydro_init(data):
    
    Hydro(data)
    
    assert True


def test_Array_init(data, turbines, features):
    
    hydro = Hydro(data)
    Array(hydro, turbines, features)
    
    return


@pytest.mark.parametrize("psi_turb, psi_yaw, psi_current, two_way, expected", [
    (5. * PIBY4, PIBY6, 0 * PIBY6, False, False),
    (5. * PIBY4, PIBY6, 1 * PIBY6, False, False),
    (5. * PIBY4, PIBY6, 2 * PIBY6, False, False),
    (5. * PIBY4, PIBY6, 3 * PIBY6, False, False),
    (5. * PIBY4, PIBY6, 4 * PIBY6, False, False),
    (5. * PIBY4, PIBY6, 5 * PIBY6, False, False),
    (5. * PIBY4, PIBY6, 6 * PIBY6, False, False),
    (5. * PIBY4, PIBY6, 7 * PIBY6, False, True),
    (5. * PIBY4, PIBY6, 8 * PIBY6, False, True),
    (5. * PIBY4, PIBY6, 9 * PIBY6, False, False),
    (5. * PIBY4, PIBY6, 10 * PIBY6, False, False),
    (5. * PIBY4, PIBY6, 11 * PIBY6, False, False),
    (5. * PIBY4, PIBY6, 0 * PIBY6, True, False),
    (5. * PIBY4, PIBY6, 1 * PIBY6, True, True),
    (5. * PIBY4, PIBY6, 2 * PIBY6, True, True),
    (5. * PIBY4, PIBY6, 3 * PIBY6, True, False),
    (5. * PIBY4, PIBY6, 4 * PIBY6, True, False),
    (5. * PIBY4, PIBY6, 5 * PIBY6, True, False),
    (5. * PIBY4, PIBY6, 6 * PIBY6, True, False),
    (5. * PIBY4, PIBY6, 7 * PIBY6, True, True),
    (5. * PIBY4, PIBY6, 8 * PIBY6, True, True),
    (5. * PIBY4, PIBY6, 9 * PIBY6, True, False),
    (5. * PIBY4, PIBY6, 10 * PIBY6, True, False),
    (5. * PIBY4, PIBY6, 11 * PIBY6, True, False),
    (0., np.pi, 0., False, True),
    (0., 0., 0., False, True),
    (np.pi, 0., np.pi, False, True),
    (2. * np.pi, 0., 0., False, True),
    (0., 0., -0.1, False, False),
    (np.pi, 0., -0.1, False, False),
    (2. * np.pi, 0., -0.1, False, False)
    ])
def test_is_within_yaw(psi_turb, psi_yaw, psi_current, two_way, expected):
    
    test = _is_within_yaw(psi_turb, psi_yaw, psi_current, two_way)
    
    assert test == expected


@pytest.mark.parametrize("psi_turb, psi_yaw, psi_current, two_way", [
    (0., -0.1, 0., False),
    (0., np.pi + 0.1, 0., False)
    ])
def test_is_within_yaw_error(psi_turb, psi_yaw, psi_current, two_way):
    
    with pytest.raises(ValueError):
        _is_within_yaw(psi_turb, psi_yaw, psi_current, two_way)


@pytest.mark.parametrize("psi_turb, psi_yaw, psi_current, two_way, expected", [
    (5. * PIBY4, PIBY6, 0 * PIBY6, False, 7 * PIBY12),
    (5. * PIBY4, PIBY6, 1 * PIBY6, False, 9 * PIBY12),
    (5. * PIBY4, PIBY6, 2 * PIBY6, False, 9 * PIBY12),
    (5. * PIBY4, PIBY6, 3 * PIBY6, False, 7 * PIBY12),
    (5. * PIBY4, PIBY6, 4 * PIBY6, False, 5 * PIBY12),
    (5. * PIBY4, PIBY6, 5 * PIBY6, False, 3 * PIBY12),
    (5. * PIBY4, PIBY6, 6 * PIBY6, False, PIBY12),
    (5. * PIBY4, PIBY6, 7 * PIBY6, False, 0.),
    (5. * PIBY4, PIBY6, 8 * PIBY6, False, 0.),
    (5. * PIBY4, PIBY6, 9 * PIBY6, False, PIBY12),
    (5. * PIBY4, PIBY6, 10 * PIBY6, False, 3 * PIBY12),
    (5. * PIBY4, PIBY6, 11 * PIBY6, False, 5 * PIBY12),
    (5. * PIBY4, PIBY6, 0 * PIBY6, True, PIBY12),
    (5. * PIBY4, PIBY6, 1 * PIBY6, True, 0.),
    (5. * PIBY4, PIBY6, 2 * PIBY6, True, 0.),
    (5. * PIBY4, PIBY6, 3 * PIBY6, True, PIBY12),
    (5. * PIBY4, PIBY6, 4 * PIBY6, True, 3 * PIBY12),
    (5. * PIBY4, PIBY6, 5 * PIBY6, True, 3 * PIBY12),
    (5. * PIBY4, PIBY6, 6 * PIBY6, True, PIBY12),
    (5. * PIBY4, PIBY6, 7 * PIBY6, True, 0.),
    (5. * PIBY4, PIBY6, 8 * PIBY6, True, 0.),
    (5. * PIBY4, PIBY6, 9 * PIBY6, True, PIBY12),
    (5. * PIBY4, PIBY6, 10 * PIBY6, True, 3 * PIBY12),
    (5. * PIBY4, PIBY6, 11 * PIBY6, True, 3 * PIBY12),
    (7. * PIBY4, 3 * PIBY4, 0., False, 0.),
    (7. * PIBY4, 3 * PIBY4, 3 * PIBY4, False, PIBY4),
    (7. * PIBY4, 3 * PIBY4, 3 * PIBY4, True, 0.),
    (3. * PIBY4, 3 * PIBY4, np.pi, False, 0.),
    (3. * PIBY4, 3 * PIBY4, 7 * PIBY4, False, PIBY4),
    (3. * PIBY4, 3 * PIBY4, 7 * PIBY4, True, 0.),
    (0., 0., 0., False, 0.),
    (0., 0., -0.1, False, 0.1),
    (2. * np.pi, 0., 0., False, 0.),
    (2. * np.pi, 0., -0.1, False, 0.1)
    ])
def test_get_angle_of_attack(psi_turb,
                             psi_yaw,
                             psi_current,
                             two_way,
                             expected):
    
    test = _get_angle_of_attack(psi_turb, psi_yaw, psi_current, two_way)
    
    assert np.isclose(test, expected)
