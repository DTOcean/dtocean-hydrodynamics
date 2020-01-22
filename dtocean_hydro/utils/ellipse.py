# -*- coding: utf-8 -*-

#    Copyright (C) 2020 Mathew Topper
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

import numpy as np


def distance(x1, x2):
    return np.sqrt((x2[0] - x1[0]) ** 2 + (x2[1] - x1[1]) ** 2)


def in_ellipse(x, a, b):
    return (x[0] ** 2 / a ** 2 + x[1] ** 2 / b ** 2) <= 1


def ellipse_intersection(x, a, b, pos=True):
    
    if pos:
        sign = 1
    else:
        sign = -1
    
    multiplier = a * b / np.sqrt(a ** 2 * x[1] ** 2 + b ** 2 * x[0] ** 2)
    
    return sign * multiplier * x[0], sign * multiplier * x[1]


def nearest_ellipse_intersection(x, a, b):
    
    pos_intersect = ellipse_intersection(x, a, b)
    neg_intersect = ellipse_intersection(x, a, b, False)
    
    if distance(pos_intersect, x) < distance(neg_intersect, x):
        result = pos_intersect
    else:
        result = neg_intersect
    
    return result


def inside_ellipse_percent(x, a, b):
    
    p1 = nearest_ellipse_intersection(x, a, b)
    percent = distance(p1, x) / distance(p1, [0, 0])
    
    return percent


def get_grid_inside_ellipse_percent(grid,
                                    angle,
                                    max_tran,
                                    max_long,
                                    atol=1e-08):
    
    """Determine whether any points in the grid lie within an ellipse placed
    around each point with semi-major axis max_tran and semi-minor axis
    max_long, rotated by angle. For each point the percentage distance along
    the transect of the closest point and the ellipse is returned."""
    
    rotator = np.array([[np.cos(angle), -np.sin(angle)],
                        [np.sin(angle),  np.cos(angle)]])
    
    f = lambda x: np.dot(rotator, x)
    g = lambda x: in_ellipse(x, max_tran, max_long)
    h = lambda x: inside_ellipse_percent(x, max_tran, max_long)
    
    result = []
    max_dist = max(max_long, max_tran)
    
    for p0 in grid:
        
        centred = grid - p0
        rotated = np.apply_along_axis(f, 1, centred)
        all_dist = np.hypot(rotated[:, 0], rotated[:, 1])
        
        logic = np.less_equal(all_dist, max_dist) & np.greater(all_dist, 0)
        suspect = rotated[logic, :]
        if suspect.shape[0] < 1: continue
        
        inside = suspect[np.apply_along_axis(g, 1, suspect), :]
        if inside.shape[0] < 1: continue
            
        percents = np.apply_along_axis(h, 1, inside)
        max_percents = max(percents)
        if max_percents < atol: continue
        
        result.append(max_percents)
    
    return np.array(result)
