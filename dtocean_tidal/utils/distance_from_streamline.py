#!/usr/bin/python2.7
# encoding: utf-8
from __future__ import division

# Start logging
import logging
module_logger = logging.getLogger(__name__)

import numpy as np
import matplotlib.pyplot as plt

def distance_from_streamline(streamline, turbPos, NbTurb, tID, debug=False, debug_plot=True):
    """
    Compute the relative distance from a given streamline (SL) and its origin.

    Args:
      streamline (tuple): SL points, tuple of 2 1D arrays, i.e. Xs and Ys in meters
      turbPos (numpy.array): turbine x, y, z coordinates, 2D array of dimension (Nturb, 3)
      NbTurb (float): number of turbines in the array
      tID (integer): turbine's ID number

    Kwargs:
      debug (bool): debug flag
      debug_plot (bool): debug plot flag

    Returns:
      turbDist (numpy.array): distance from hub to SL, 2D array of dimension (Xturb, 3)
      where each vector is composed of: [turbine id, distance along SL, distance across SL]

    """
    if debug: module_logger.info("Computing relative distance from hub...")
    L = len(streamline[0])
    # Create vectors from streamline
    slVect = np.zeros((L-1,2))
    X = np.asarray(streamline[0])
    Y = np.asarray(streamline[1])
    slVect[:,0] = X[:-1] - X[1:]
    slVect[:,1] = Y[:-1] - Y[1:]
    # filtering Nans
    nanXind = np.where(np.isnan(slVect[:,0]))[0]
    nanYind = np.where(np.isnan(slVect[:,1]))[0]
    slVect[nanXind,0] = slVect[nanXind-1,0]
    slVect[nanYind,1] = slVect[nanYind-1,1]
    # Compute angle and length
    turbVect = np.zeros(slVect.shape)
    turbDist = {}

    for i in range(NbTurb):
        turbVect[:,0] = turbPos[i,0] - X[:-1]
        turbVect[:,1] = turbPos[i,1] - Y[:-1]
        # Scalar product of vectors
        dotProd = np.sum(slVect * turbVect, 1)
        # Search for scalar product = 0 == perpendicular point
        zeroCross = np.where(np.diff(np.sign(dotProd)))[0]

        if not zeroCross.shape[0]==0:
            # In case several distance possible:
            if zeroCross.shape[0]>1:
                zeroCross = zeroCross.max()
            # Check if it is not referring to itself
            if not i == tID:
                # Real coordinates of perpendicular point
                tot = np.abs(dotProd[zeroCross]) + np.abs(dotProd[zeroCross+1])
                weightA = np.abs(dotProd[zeroCross]) / tot
                # Real distance in between + and - scalar product step
                dAlongSL = np.sum(np.hypot(slVect[:zeroCross,0],
                                           slVect[:zeroCross,1]))\
                         + weightA*(np.hypot(slVect[zeroCross,0],
                                             slVect[zeroCross,1]))
                X1 = X[(zeroCross)] + weightA*(X[(zeroCross+1)]-X[(zeroCross)])
                Y1 = Y[(zeroCross)] + weightA*(Y[(zeroCross+1)]-Y[(zeroCross)])
                dAcrossSL =np.hypot(np.abs(turbPos[i,0]-X1),np.abs(turbPos[i,1]-Y1))
                if not dAlongSL == 0.0: # Sanity check.
                    try:
                        dist = np.asarray([dAlongSL[0], dAcrossSL[0]])
                    except IndexError:
                        dist = np.asarray([dAlongSL, dAcrossSL])
                    turbDist[str(i)] = dist
                    if debug_plot:
                        x = [turbPos[i,0],X1]
                        # print x
                        y = [turbPos[i,1],Y1]
                        # print y
                        plt.plot(x, y, 'k--', alpha=0.2)
    return turbDist