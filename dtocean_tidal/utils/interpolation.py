# -*- coding: utf-8 -*-

#    Copyright (C) 2016 Thomas Roc
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

from __future__ import division

import logging

import numpy as np
import scipy.spatial
from scipy.interpolate import LinearNDInterpolator

# Start logging
module_logger = logging.getLogger(__name__)


def interp_at_point(x, y, X, Y, Q):
    """
    Compute the quantity Q at point (x,y)

    Args:
      x (float): x axis coordinate, float
      y (float): y axis coordinate, float
      X (numpy.array): x axis coordinates, float array (N)
      Y (numpy.array): y axis coordinates, float array (M)
      Q (list): given list quantity matrix, list of float arrays (M,N)

    Returns:
      qi (list): interpolated values

    """
    # find first two nearest
    xDist = x-X
    i = np.argmin(np.abs(xDist))
    try:
        if np.sign(xDist[i]) < 0. and i - 1 < 0: raise IndexError
        if np.abs(xDist[i + 1]) < np.abs(xDist[i - 1]):
            i2 = i + 1
        else:
            i2 = i - 1
        # computes weight
        dist = np.abs(xDist[i]) + np.abs(xDist[i2])
        ai = np.abs(xDist[i2] / dist)
        ai2 = np.abs(xDist[i] / dist)
    except IndexError:  # the given point is outside the Q array
        i2 = 0
        ai = 1.0
        ai2 = 0.0
    yDist = y-Y
    j = np.argmin(np.abs(yDist))
    try:
        if np.sign(yDist[j]) < 0. and j - 1 < 0: raise IndexError
        if np.abs(yDist[j + 1]) < np.abs(yDist[j - 1]):
            j2 = j + 1
        else:
            j2 = j - 1
        # computes weight
        dist = np.abs(yDist[j]) + np.abs(yDist[j2])
        aj = np.abs(yDist[j2] / dist)
        aj2 = np.abs(yDist[j] / dist)
    except IndexError:  # the given point is outside the Q array
        j2 = 0
        aj = 1.0
        aj2 = 0.0

    # Bilinear interpolation
    qi = []
    
    for q in Q:
        
        # Don't allow NaN values
        bi = 0

        if not np.isnan(q[j, i]): bi += q[j, i] * aj * ai
        if not np.isnan(q[j, i2]): bi += q[j, i2] * aj * ai2
        if not np.isnan(q[j2, i]): bi += q[j2, i] * aj2 * ai
        if not np.isnan(q[j2, i2]): bi += q[j2, i2] * ai2 * aj2

        qi.append(bi)

    return qi


def volume_under_plane(x,y,z, debug=False):
    """
    Computes the volume under a given plane

    Args:
      x (numpy.array): x coordinates in m, 1d array
      y (numpy.array): y coordinates in m, 1d array
      z (numpy.array): z coordinates in m, 2d array, (y,x)

    Kwargs:
      debug (bool): debug flag

    Returns:
      volume (float): volume in m3, float
    """
    if debug: module_logger.info("Computing volume under plane...")
    X, Y = np.meshgrid(x, y)
    X = X.flatten()
    Y = Y.flatten()
    Z = z.flatten()
    xyz = np.vstack((X,Y,Z)).T

    d = scipy.spatial.Delaunay(xyz[:,:2])
    tri = xyz[d.vertices]

    a = tri[:,0,:2] - tri[:,1,:2]
    b = tri[:,0,:2] - tri[:,2,:2]
    proj_area = np.cross(a, b).sum(axis=-1)
    zavg = tri[:,:,2].sum(axis=1)
    volume = zavg * np.abs(proj_area) / 6.0

    return volume

def interpol_scatter2grid(x, y, z, xn, yn, debug=True):
    """
    Interpolate scattered points onto grid.

    Args:
      x (numpy.array): x coordinates, 1d array, n floats
      y (numpy.array): y coordinates, 1d array, n floats
      z (numpy.array): values, 1d array
      xn (numpy.array): x dimension of interpolated array, integer
      yn (numpy.array): y dimension of interpolated array, integer

    Kwargs:
      debug (bool): debug flag

    Returns:
      xi (numpy.array): gridded x coordinates, 1d array (xn)
      yi (numpy.array): gridded y coordinates, 1d array (yn)
      zi (numpy.array): gridded z values, 2d array (yn, xn)

    """
    if debug: module_logger.info("Interpolating onto grid...")
    xmin = x.min()
    xmax = x.max()
    ymin = y.min()
    ymax = y.max()
    dx = (xmax - xmin) / (xn - 1)
    dy = (ymax - ymin) / (yn - 1)
    xi = np.arange(xmin, xmax + dx, dx)
    yi = np.arange(ymin, ymax + dy, dy)

    X, Y = np.meshgrid(xi, yi)
    orig = np.vstack((x, y)).T
    asked = np.vstack((X.flatten(), Y.flatten())).T

    interpol = LinearNDInterpolator(orig, z.flatten())
    zi = interpol(asked)

    return xi, yi, zi.reshape(X.shape)