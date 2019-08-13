# -*- coding: utf-8 -*-

#    Copyright (C) 2016 Thomas Roc
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
    
    transect, first_row, speed = transec_surf(hydro, array, debug=debug)
    rotor_surf = 0.0
    
    n_digits = len(str(array.turbine_count))
    
    for i in first_row:
        turb_name = 'turbine{:0{width}d}'.format(i, width=n_digits)
        diam = array.features[turb_name]['Diam']
        ry = array.features[turb_name]['RY']  # relative yawing angle
        surf = np.pi * ((diam/2.0)**2.0) * abs(np.cos(np.radians(ry)))  # ellipse area
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
