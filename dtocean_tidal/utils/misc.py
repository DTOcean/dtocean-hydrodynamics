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

import re
import math
import logging

import numpy as np
from numpy.linalg import LinAlgError
from shapely.geometry import LineString

# local imports
from .interpolation import interp_at_point

# Start logging
module_logger = logging.getLogger(__name__)


def closest_point( pt, pts, debug=False):
    """
    Finds the closest point in 'pts'
    and its associated distance
    to any given [x,y] 'pt'.

    Args:
      pt (list): [x,y]
      pts (list): points, [[x1,y1],...[xn,yn]]

    Kwargs:
      debug (bool): debug flag

    Returns:
      closest_point (list): [xi,yi]
      dist (float): distance, (m)

    """
    if debug: module_logger.info("Computing closest point index...")

    dist = np.inf
    ind = 0
    for test in pts:
        d = np.sqrt((test[0] - pt[0])**2.0 + (test[1] - pt[1])**2.0)
        if d < dist:
            closest_point = test
            dist = d
            index = ind
        ind += 1

    if debug: module_logger.info("index: {}".format(closest_point))
    if debug: module_logger.info("distance: {}".format(dist))

    return closest_point, index, dist


def line(p1, p2):
    """
    Defines line.

    Args:
      p1 (list): points, [x1,y1]
      p2 (list): points, [x2,y2]

    Returns:
      list: [A, B, C]
    """
    A = (p1[1] - p2[1])
    B = (p2[0] - p1[0])
    C = (p1[0]*p2[1] - p2[0]*p1[1])
    return [A, B, -C]


def intersection(L1, L2):
    """
    Computes the intersection point x,y between two lines

    Args:
      L1 (list): output of the dtocean_tidal.utils.line function
      L2 (list): output of the dtocean_tidal.utils.line function

    Returns:
      x (float): x coordinate of the intersection between L1 & L2
      y (float): y coordinate of the intersection between L1 & L2
    """
    D  = L1[0] * L2[1] - L1[1] * L2[0]
    Dx = L1[2] * L2[1] - L1[1] * L2[2]
    Dy = L1[0] * L2[2] - L1[2] * L2[0]
    if D != 0:
        x = Dx / D
        y = Dy / D
        return x,y
    else:
        return False


def natural_sort(l):
    """Sort string list based on number in string

    Args:
      l (list): list strings with numbers

    Returns:
      list: sorted list
    """
    convert = lambda text: int(text) if text.isdigit() else text.lower() 
    alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ] 
    return sorted(l, key = alphanum_key)


def radians_to_bearing(x):

    initial_bearing = 90 - math.degrees(x)
    compass_bearing = (initial_bearing + 360) % 360

    return compass_bearing


def vector_to_bearing(x, y):

    initial_bearing = math.atan2(y, x)
    compass_bearing = radians_to_bearing(initial_bearing)

    return compass_bearing


def pi2pi(angle):
    """ Ensures angle is in the interval (-pi, pi]

    Args:
      angle (float): angle in radians

    Returns:
      newDir (float): angle in radians in interval (-pi, pi]
    """
    
    while True:
        
        if angle > np.pi:
            
            angle -= 2 * np.pi
            continue
        
        elif angle < -np.pi or np.isclose(angle, -np.pi):
            
            angle += 2 * np.pi
            continue
        
        break

    return angle


def bearing_to_radians(bearing):
    """converts from 0-to-360 deg (North = 0, East = 90) to -pi-to-pi radian
    (North = pi/2; East = 0)

    Args:
      dir (float): angle in degree

    Returns:
      newDir (float): angle in radian

    """
    angle = np.radians(np.mod(90 - bearing, 360))
    angle = pi2pi(angle)
    
    return angle


def transec_surf(hydro, array, debug=False):
    """
    Computes transect surface (m2) and first row of the array
    
    Args:
      hydro (dtocean_tidal.main.Hydro): dtocean_tidal's Hydro object
      array (dtocean_tidal.main.Array): dtocean_tidal's Array object
    
    Kwargs:
      debug (bool): debug flag
    
    Returns:
      transect (float): lease's transect surface (m2), float
      first_row (list): first row of the array, list of ID numbers
      speed (float): averaged speed over transect (m/s), float
    
    """
    
    if debug: module_logger.info("Computing relative blockage ratio RBR...")
    
    Nturb = len(array.positions.keys())
    n_digits = len(str(Nturb))
    
    (xm, ym, xM, yM) = hydro.lease.bounds
    l = []
    ref = np.arange(Nturb)

    # In case there is only one turbine
    if Nturb == 1:
        l = array.positions.keys()[7:]  # removing 'turbine' from key
    else:
        for k1 in array.distances.keys():
            for k2 in array.distances[k1].keys():
                # find turbine away from 1 diam
                diam = array.features[k1]['Diam']
                # check distance along streamline
                if array.distances[k1][k2][0] > diam:
                    l.append(k2)
    
    # unique occurrence
    L = np.asarray(l).astype(int)
    un = np.unique(L)
    
    # first row composed of
    first_row = list(set(ref)-set(un))
    
    # Fix for odd cases/flows where no obvious first raw
    # First_row = to turbine with the least interactions (n.b.: first raw should not have any)
    if first_row == []:
        unique, counts = np.unique(L, return_counts=True)
        first_row = list(unique[np.where(counts == counts.min())])
    
    # check for missing turbine in the first row within a diam radius
    iterlist = first_row[:]
    
    for ii in iterlist:
        
        ii_turb_name = 'turbine{:0{width}d}'.format(ii, width=n_digits)
        
        diam = array.features[ii_turb_name]['Diam']
        ii_x = array.positions[ii_turb_name][0]
        ii_y = array.positions[ii_turb_name][1]
        
        for jj in range(Nturb):
            
            jj_turb_name = 'turbine{:0{width}d}'.format(jj, width=n_digits)
            
            jj_x = array.positions[jj_turb_name][0]
            jj_y = array.positions[jj_turb_name][1]
            xdist = np.abs(ii_x - jj_x)
            ydist = np.abs(ii_y - jj_y)
            
            if np.sqrt(xdist ** 2.0 + ydist ** 2.0) < diam:
                first_row.append(jj)
    
    first_row = np.unique(first_row)
    
    first_turb_name = 'turbine{:0{width}d}'.format(first_row[0],
                                                   width=n_digits)
    last_turb_name = 'turbine{:0{width}d}'.format(first_row[-1],
                                                  width=n_digits)
    
    # define linear function representative on the 1st row's ray line
    x1 = array.positions[first_turb_name][0]
    y1 = array.positions[first_turb_name][1]
    x2 = array.positions[last_turb_name][0]
    y2 = array.positions[last_turb_name][1]
    a = np.array([[x1,1], [x2,1]])
    b = np.array([y1,y2])
    
    try:
        coeffs = np.linalg.solve(a, b)
        ybm = coeffs[0] * xm + coeffs[1]
        ybM = coeffs[0] * xM + coeffs[1]
        xbm = xm
        xbM = xM
    except LinAlgError:  # error when row parallel of perpendicular to x axis
        if x1 == x2:
            ybm = ym
            ybM = yM
            xbm = x1
            xbM = x1
        if y1==y2:
            ybm = y1
            ybM = y1
            xbm = xm
            xbM = xM
    
    # define lineString within lease
    box = hydro.lease
    line = LineString([[xbm,ybm],[xbM,ybM]])
    inter_line = box.intersection(line)
    
    # calculate ratio between rotors surface and transect
    diam = array.features[last_turb_name]['Diam']
    n = inter_line.length//diam
    dl = inter_line.length / (n-1)
    
    # define points uniformly along line
    pts = map(inter_line.interpolate, np.arange(n) * dl)
    
    #  interpolate bathymetry at those points
    Q = [hydro.SSH, hydro.bathy, hydro.U, hydro.V]
    
    wh = 0.0  # water column height
    speed = 0.0  # speed over transect
    
    for pt in pts:
        [el, h, u, v] = interp_at_point(pt.x, pt.y, hydro.X, hydro.Y, Q)
        # Don't allow NaNs
        if np.isnan([el, h, u, v]).any(): continue
        # quantity formatting
        if h < 0.0: h *= -1.0
        wh += (h + el)
        speed += np.sqrt(u**2.0 + v**2.0)
    
    # relative blockage ratio
    transect = (wh / n) * inter_line.length  # average transect surface
    speed = speed / n

    return transect, first_row, speed
