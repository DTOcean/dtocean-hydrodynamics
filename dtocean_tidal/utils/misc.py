#!/usr/bin/python2.7
# encoding: utf-8
from __future__ import division

# Start logging
import logging
module_logger = logging.getLogger(__name__)

import numpy as np
import re
from shapely.geometry import LineString
from numpy.linalg import LinAlgError

# local imports
from dtocean_tidal.utils.interpolation import interp_at_point


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


def pi2pi(dir):
    """ keeps angle between -pi and pi

    Args:
      dir (float): angle in radian

    Returns:
      newDir (float): angle in radian between -pi and pi
    """
    if dir > np.pi:
        newDir = (dir-np.pi)-np.pi
    elif dir < -np.pi:
        newDir = (dir+np.pi)+np.pi
    else :
        newDir = dir
    return newDir


def deg360_to_radpi(dir):
    """converts from 0-to-360 deg (North = 0, East = 90) to -pi-to-pi radian (North = pi/2; East = 0)

    Args:
      dir (float): angle in degree

    Returns:
      newDir (float): angle in radian

    """
    dir = np.radians(np.mod(dir - 90, 360))
    newDir = pi2pi(dir)
    return newDir


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
                if array.distances[k1][k2][0] > diam:  # check distance along streamline
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
        for jj in range(Nturb):
            diam = array.features['turbine'+str(ii)]['Diam']
            xdist = np.abs(array.positions['turbine'+str(ii)][0] - array.positions['turbine'+str(jj)][0])
            ydist = np.abs(array.positions['turbine'+str(ii)][1] - array.positions['turbine'+str(jj)][1])
            if np.sqrt(xdist**2.0 + ydist**2.0) < diam:
                first_row.append(jj)
    first_row = np.unique(first_row)

    # define linear function representative on the 1st row's ray line
    x1 = array.positions['turbine'+str(first_row[0])][0]
    y1 = array.positions['turbine'+str(first_row[0])][1]
    x2 = array.positions['turbine'+str(first_row[-1])][0]
    y2 = array.positions['turbine'+str(first_row[-1])][1]
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
    box = hydro.bounding_box
    line = LineString([[xbm,ybm],[xbM,ybM]])
    inter_line = box.intersection(line)

    # calculate ratio between rotors surface and transect
    diam = array.features['turbine'+str(first_row[-1])]['Diam']
    n = inter_line.length//diam
    dl = inter_line.length / (n-1)
    pts = map(inter_line.interpolate, np.arange(n) * dl)  # define points uniformly along line
    #  interpolate bathymetry at those points
    Q = [hydro.SSH, hydro.bathy, hydro.U, hydro.V]
    wh = 0.0  # water column height
    speed = 0.0  # speed over transect
    for pt in pts:
        [el, h, u, v] = interp_at_point(pt.xy[0], pt.xy[1], hydro.X, hydro.Y, Q)
        # quantity formatting
        if h < 0.0: h *= -1.0
        wh += (h + el)
        speed += np.sqrt(u**2.0 + v**2.0)
    # relative blockage ratio
    transect = (wh / n) * inter_line.length  # average transect surface
    speed = speed / n

    return transect, first_row, speed