#!/usr/bin/python2.7
# encoding: utf-8
from __future__ import division

# Start logging
import logging
module_logger = logging.getLogger(__name__)

import numpy as np

# local imports
from dtocean_tidal.utils.misc import transec_surf


def blockage_ratio(hydro, array, debug=False):
    """
    Computes blockage ratio relative to the lease area

    Args:
      hydro (dtocean_tidal.main.Hydro): dtocean_tidal's Hydro object
      array (dtocean_tidal.main.Array): dtocean_tidal's Array object

    Kwargs:
      debug (bool): debug flag

    Returns:
      rbr (float): relative blockage ration (i.e. rotor's surface / lease's transect surface)

    """
    if debug: module_logger.info("Computing relative blockage ratio RBR...")
    transect, first_row, speed = transec_surf(hydro,array, debug=debug)
    rotor_surf = 0.0
    for i in first_row:
        diam = array.features['turbine'+str(i)]['Diam']
        ry = array.features['turbine'+str(i)]['RY']  # relative yawing angle
        surf = np.pi * ((diam/2.0)**2.0) * np.cos(np.radians(ry))  # ellipse area
        rotor_surf += surf

    # final check
    rbr = rotor_surf/transect
    if rbr > 1.0:
        rbr = 1.0
    if debug:
        module_logger.info("...transect surface = " + str(transect) + " m2...")
        module_logger.info("...rotors surface = " + str(rotor_surf) + " m2...")
        module_logger.info("...RBR = " + str(rbr))

    return rbr