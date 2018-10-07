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

# Start logging
import logging
module_logger = logging.getLogger(__name__)

import numpy as np

def vvpw(u, v, z, el, h, n, debug=False):
    """
    Computes the Vertical Velocity Profile Weight

    Args:
      u (float): velocity component in the x direction (m/s)
      v (float): velocity component in the y direction (m/s)
      z (float): hub height (m), positive from sea bottom upwards
      el (float): sea surface elevation (m)
      h (float): depth (m), must be positive
      n (float): Manning coefficient

    Kwargs:
      debug (bool): debug flag

    Returns:
      w (float): Vertical Velocity Profile Weight

    References:
      - M. Reiner and R. Schoenfeld-Reiner, Viskosimetrische Untersuchungen an Lösungen hochmolekularer Naturstoffe.
        I. Mitteilung. Kautschuk in Toluol., Kolloid-Zeitschrift, vol. 65, no. 1, pp. 44-62, 1933.
      - H.-E. Lee, C. Lee, Y.-J. Kim, J.-S. Kim and W. Kim,
        Power Law Exponents for Vertical Velocity Distributions in Natural Rivers,
        Engineering, vol. 5, pp. 933-942, 2013.

    """
    if debug: module_logger.info("Computing the Vertical Velocity Profile Weight...")

    # Constants
    k = 0.41  # Von Karman constant
    g = 9.81  # gravitational constant
    nu = 1.0e-4  # kinematic viscosity (m2/s)

    # Variables
    U = np.sqrt(u**2.0 + v**2.0) # flow speed
    H = (el + h)
    Re = (U * H) / nu
         
    # If the Reynolds number is zero return 1.
    if np.isclose(Re, 0.): return 1.
         
    m = k * np.sqrt((Re**(1.0/3.0)) / ((n**2.0) * g))  # exponent from Lee et al. 2013
    C = (H**((1.0/m) - 1.0)) / ((1.0/m) + 1.0)

    # Vertical Velocity Profile Weight
    w = (z / (C * H))**(1.0/m)
    if debug: module_logger.info("....Vertical Velocity Profile Weight = " + str(w) + "...")
    if debug: module_logger.info("...done.")
    
    if np.isnan(w):
        errStr = "Vertical velocity profile weight is NaN"
        raise ArithmeticError(errStr)

    return w

def vvpw_soulsby(z, el, h, z0, ple, debug=False):
    """
    Computes the Vertical Velocity Profile Weight

    Args:
      z (float): hub height (m), positive from sea bottom upwards
      el (float): sea surface elevation (m)
      h (float): depth (m), must be positive
      z0 (float): bed roughness, usually between 0 and 0.32
      ple (float): power law exponent, usually = 7.0

    Kwargs:
      debug (bool): debug flag

    Returns:
      w (float): Vertical Velocity Profile Weight

    References:
      - R. Soulsby, Dynamics  of  Marine  Sands:  A  Manual  for  Practical Applications
        London: Telford, 249p. 1997.
    """
    if debug: module_logger.info("Computing the Vertical Velocity Profile Weight...")

    # Water column height
    H = (el + h)

    # Vertical Velocity Profile Weight
    w = (z/(z0*H))**(1.0/ple)
    if debug: module_logger.info("....Vertical Velocity Profile Weight = " + str(w) + "...")
    if debug: module_logger.info("...done.")

    return w