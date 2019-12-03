# -*- coding: utf-8 -*-

#    Copyright (C) 2016 Francesco Ferri
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
This module contains the class used to collect the parameters used in the
creation of the tidal or wave numerical models, and characterisaiton of the sea conditions

.. module:: hydro
   :platform: Windows
   :synopsis: Hydro module for the DTOcean WP2

.. moduleauthor:: Francesco Ferri <ff@civil.aau.dk>
"""

from math import pi
from dtocean_wave.utils.WatWaves import WNumber


class Hydro_pkg(object):
    """
    Args:
        the number of args is variable, depending on the machine type:
        Tidal: 10 inputs
            DeviceType (string): either a "T" or "W" for tidal or wave respectively
            V (numpy.ndarray): y-contribution (Northing) to the velocity field at each grid node
            U (numpy.ndarray): x-contribution (Easting) to the velocity field at each grid node
            p (numpy.ndarray): probability of occurency of the different sea states.
            TI (numpy.ndarray): turbulence intensity at each grid node for each sea state
            x (numpy.ndarray): x-coordinates (Easting) of the grid nodes
            y (numpy.ndarray): y-coordinates (Northing) of the grid nodes
            wdepth (numpy.ndarray): sea surface elevation at each grid node for each sea state
            bathy (numpy.ndarray): bathymetry of the lease area at each grid node
            beta (float): bed roughness
            alpha (float): power law exponent
        Wave: 9 inputs
            DeviceType (string): : either a "T" or "W" for tidal or wave respectively
            B (numpy.ndarray): Sea state wave directions.
            Hs (numpy.ndarray): Sea state wave heights.
            Tp (numpy.ndarray): Sea state wave periods.
            ScatDiag (numpy.ndarray): probability of occurency of the different sea states (B, Hs, Tp).
            depth (float): Flattened water depth of the site.
            cfrequency (numpy.ndarray): frequencies (rad/s) used to deiscretise the frequency domain model.
            dirs (numpy.ndarray): wave directions (rad) used to discretise the freqeucny domain model.
            specType (tuple): description of the wave spectrums, spectrum type, peak factor and spreding parameter

    Attributes:
        Same as the Args list
        for the wave case:
            .period: calculated as 2*pi/cfrequency
            .wnumber: calculated from the period and the water depth
    """
    def __init__(self, *arg):
        if arg[0]=='W':
            B = arg[1]
            Hs = arg[2] 
            Tp = arg[3]
            ScatDiag = arg[4]
            self.depth = arg[5]
            freqs = arg[6]
            wdirs = arg[7]
            self.specType = arg[8]
            self.B = B
            self.Hs = Hs
            self.Tp = Tp
            self.ScatDiag = ScatDiag
            self.cfrequency = freqs
            self.period = 2*pi/freqs 
            self.dir = wdirs
            self.wnumber= WNumber(self.period, self.depth)
        elif arg[0]=='T':
            self.V = arg[1]
            self.U = arg[2]
            self.p = arg[3]
            self.TI = arg[4]
            self.x = arg[5]
            self.y = arg[6]
            self.wdepth = arg[7]
            self.bathy = arg[8]
            self.beta = arg[9]
            self.alpha = arg[10]
        else:
            raise IOError('ERROR: Invalid input list to the Hydro class.')
