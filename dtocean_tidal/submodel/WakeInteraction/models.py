# -*- coding: utf-8 -*-

#    Copyright (C) 2019 Mathew Topper
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

import numpy as np


class DominantWake(object):
    
    # Returns the velocity coefficient and indexes of the dominant wake 
    # affecting each turbine
    
    def __init__(self, turb_velocity,
                       wake_matrix):
        
        turb_speed = np.sqrt(turb_velocity[0, :] ** 2 +
                                                 turb_velocity[1, :] ** 2)
        
        self._coefficients = get_wake_coefficients(turb_speed,
                                                   wake_matrix)
        self.indexes = np.argmin(self._coefficients, axis=0)
        
        return
    
    @property
    def coefficients(self):
        
        coefficients = self.get_dominant(self._coefficients)
        
        return coefficients
    
    def get_dominant(self, matrix):
        
        if matrix.shape != self._coefficients.shape:
            
            err_str = ("Supplied matrix must have shape "
                       "{0}x{0}").format(self._coefficienst.shape)
            raise ValueError(err_str)
        
        dominant = matrix[self.indexes,
                          range(self._coefficients.shape[1])]
        
        return dominant


def get_wake_coefficients(baseline,
                          interaction_matrix):
    
    coefficient = interaction_matrix / baseline[:, None]
    
    return coefficient
