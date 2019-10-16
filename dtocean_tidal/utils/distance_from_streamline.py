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

# Start logging
module_logger = logging.getLogger(__name__)


def distance_from_streamline(streamline,
                             turbine_positions,
                             debug=False,
                             debug_plot=True):
    """
    Compute the relative distance from a given streamline (SL) and its origin.
    
    Args:
      streamline (tuple): SL points, tuple of 2 1D arrays, i.e. Xs and Ys in
          meters
      turbine_positions (numpy.array): turbine x, y, z coordinates, 2D array of
          dimension (Nturb, 3)
      n_turbines (float): number of turbines in the array
      turbine_id (integer): turbine's ID number
    
    Kwargs:
      debug (bool): debug flag
      debug_plot (bool): debug plot flag
    
    Returns:
      turbine_distances (dict): distance from hub to SL, 2D array of dimension
           (Xturb, 3) where each item is composed of:
           turbine id: np.array(distance along SL, distance across SL)
    
    """
    if debug: module_logger.info("Computing relative distance from hub...")
    
    # Initialise output
    turbine_distances = {}
    
    # Split streamline into x and y coordinates
    X = np.asarray(streamline[0])
    Y = np.asarray(streamline[1])
    
    # Edge cases: No streamline or single point
    if len(X) < 2: return turbine_distances
    
    # Create vectors from streamline point n to point n+1
    streamline_vector = np.zeros((len(X) - 1, 2))
    streamline_vector[:, 0] = X[1:] - X[:-1]
    streamline_vector[:, 1] = Y[1:] - Y[:-1]
    
    # filtering Nans
    streamline_vector[~np.isnan(streamline_vector).any(axis=1)]
    
    # Compute angle and length
    turbine_vector = np.zeros(streamline_vector.shape)
    
    n_turbines = turbine_positions.shape[0]
    
    for i in range(n_turbines):
        
        # Calculate vector from turbine to each streamline point
        turbine_vector[:, 0] = turbine_positions[i, 0] - X[:-1]
        turbine_vector[:, 1] = turbine_positions[i, 1] - Y[:-1]
        
        # Another edge case when the turbine is at the beginning of the
        # streamline
        if np.isclose(np.sum(turbine_vector[0, :]), 0.): continue
        
        # Inner product of vectors
        inner_product = np.sum(streamline_vector * turbine_vector, 1)
        segment_lengths = np.sum(np.square(streamline_vector), 1)
        
        # The point lines perpendicular if there is some point such that the
        # inner product is greater than zero or less than than the segment
        # length squared
        check_zero = 0. < inner_product
        check_segment = inner_product <= segment_lengths
        check_both = np.logical_and(check_zero, check_segment)
        
        # If the turbine doesn't project onto any of the segments then loop.
        if not check_both.any(): continue
        
        # Filter the valid segments
        valid_segments = streamline_vector[check_both, :]
        valid_product = inner_product[check_both]
        valid_lengths = segment_lengths[check_both]
        valid_startx = X[:-1][check_both]
        valid_starty = Y[:-1][check_both]
        
        # Calculate perpendicular distances
        segment_weights = valid_product / valid_lengths
        
        perpendicular_x = valid_startx - turbine_positions[i, 0] + \
                                     segment_weights * valid_segments[:, 0]
        perpendicular_y = valid_starty - turbine_positions[i, 1] + \
                                     segment_weights * valid_segments[:, 1]
        perdendicual_distance = np.hypot(perpendicular_x, perpendicular_y)
        
        # Find nearest perpendicular point on the streamline
        min_dist = np.min(perdendicual_distance)
        min_idx = np.argmin(perdendicual_distance)
        min_segment = valid_segments[min_idx, :]
        min_weight = segment_weights[min_idx]
        
        # Collect segments leading up to nearest
        where_segments = np.where(streamline_vector == min_segment)
        upto_segments = streamline_vector[:where_segments[0][0], :]
        
        # Find distance along streamline to the nearest point
        upto_distance = np.sum(np.sqrt(np.sum(np.square(upto_segments), 1)))
        part_distance = np.sqrt(np.sum(np.square(min_segment))) * min_weight
        streamline_distance = upto_distance + part_distance
        
        turbine_distances[i] = np.array([streamline_distance,
                                         min_dist])
    
    if debug: module_logger.info(turbine_distances)
    
    return turbine_distances
