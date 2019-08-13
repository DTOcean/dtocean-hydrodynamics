# -*- coding: utf-8 -*-

#    Copyright (C) 2016 Thomas Roc
#    Copyright (C) 2017-2019 Mathew Topper
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

# Start logging
module_logger = logging.getLogger(__name__)


def vvpw(z, el, h, beta, alpha, debug=False):
    """
    Computes the vertical velocity profile weight, with a fixed power
    law expression.
    
    Args:
      z (float): hub height (m), positive from sea bottom upwards
      el (float): sea surface elevation (m)
      h (float): depth (m), must be positive
      beta (float): bed roughness, usually between 0 and 0.32
      alpha (float): power law exponent, usually = 7.0
    
    Kwargs:
      debug (bool): debug flag
    
    Returns:
      w (float): Vertical Velocity Profile Weight
    
    References:
      - R. Soulsby, Dynamics  of  Marine  Sands:  A  Manual  for  Practical 
        Applications,
        London: Telford, 249p. 1997.
      - M. Lewis, S. P. Neill, P. Robins, M. R. Hashemi, and S. Ward, 
        'Characteristics of the velocity profile at tidal-stream energy sites', 
        Renewable Energy, vol. 114, pp. 258–272, Dec. 2017.

    """
    
    total_height = (el + h)
    
    if alpha <= 0. or beta <= 0. or total_height <= 0.:
        return 0.
    
    w = (z / (beta * total_height)) ** (1.0 / alpha)
    
    if debug:
        log_msg = "Vertical velocity profile weight = {}".format(w)
        module_logger.debug(log_msg)
    
    return w
