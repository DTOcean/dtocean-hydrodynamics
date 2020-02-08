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

import time
import logging

import numpy as np
import matplotlib.pyplot as plt
from numpy.linalg import norm
from descartes import PolygonPatch

# Local import
from .models import DominantWake, get_wake_coefficients
from ..ParametricWake import WakeShape, Wake
from ...modules.blockage_ratio import blockage_ratio

# Start logging
module_logger = logging.getLogger(__name__)


class WakeInteraction:
    """
    -WakeInteraction class-
    
    Computes wake interaction and their impacts of flow field
    
    """
    
    def __init__(self, hydro,
                       array,
                       cfd_data,
                       U_dict,
                       V_dict,
                       TKE_dict,
                       criterior=1e-8,
                       max_loop=50,
                       debug=False,
                       debug_plot=False):
        
        self._turbine_count = len(array.positions.keys())
        self._debug = debug
        self._debug_plot = debug_plot
        self._hydro = hydro
        self._bounding_box = hydro.bounding_box
        self._array = array
        self._df = cfd_data
        self._U_dict = U_dict
        self._V_dict = V_dict
        self._TKE_dict = TKE_dict
        
        nan_array = np.empty(self._turbine_count)
        nan_array[:] = np.nan
        
        self.coefficient = np.ones(self._turbine_count)
        self.inducedTKE = nan_array
        self._wake = {}
        self._wakeShape = {}
        
        rbr = blockage_ratio(hydro, array, debug=debug)
        self._blockage = rbr * hydro.BR
        
        self._criterior = criterior
        self._max_loop = max_loop
        
        return
    
    def solve_flow(self, debug=False):
        
        """
        Find velocities and TIs at each turbine by iterating the velocity and
        TI coefficients based on superposition of wakes
        
        Kwargs:
          criterior (float): convergence criterior
          max_loop (int): maximum number of loops
          debug (bool): debug flag
        """
        
        debug = debug or self._debug
        
        iniSpeed = np.empty(self._turbine_count)
        iniVel = np.empty((2, self._turbine_count))
        iniTI = np.empty(self._turbine_count)
        iniTKE = np.empty(self._turbine_count)
        
        n_digits = len(str(self._turbine_count))
        
        for i in range(self._turbine_count):
            
            turb = 'turbine{:0{width}d}'.format(i, width=n_digits)
            
            p = self._array.velHub[turb]
            
            iniSpeed[i] = np.sqrt(p[0]**2.0 + p[1]**2.0)
            iniVel[:, i] = self._array.velHub[turb][:]
            iniTI[i] = self._array.features[turb]['TIH']
            iniTKE[i] = _get_tke(iniTI[i], iniSpeed[i])
            
            self._wake[turb] = Wake(self._df,
                                    self._U_dict,
                                    self._V_dict,
                                    self._TKE_dict,
                                    self._array.features[turb],
                                    self._blockage,
                                    debug=self._debug)
        
        if debug:
            module_logger.info("Querying parametric wake model...")
            start = time.time()
        
        ind_err = np.inf
        old_coefficient = np.ones(self._turbine_count)
        newVel = iniVel.copy()
        newTI = iniTI.copy()
        newTKE = iniTKE.copy()
        loop_counter = 0
        
        while (ind_err > self._criterior and
               loop_counter < self._max_loop):
            
            (newVel,
             newSpeed,
             newTI,
             newTKE) = _solve_flow(self._turbine_count,
                                   self._array.distances,
                                   self._wake,
                                   newVel,
                                   newTI,
                                   newTKE,
                                   iniVel,
                                   iniTKE,
                                   debug)
            
            new_coefficient = newSpeed / iniSpeed
            ind_err = norm(abs(old_coefficient - new_coefficient))
            old_coefficient = new_coefficient
            
            loop_counter += 1
        
        if debug:
            end = time.time()
            module_logger.info("...done after " + str(end - start) +
                                                               " seconds.")
        
        if ind_err <= self._criterior:
            
            log_msg = ("Device interaction solver converged in {} "
                       "iteration(s)").format(loop_counter)
            module_logger.debug(log_msg)
        
        else:
            
            log_msg = ("Device interaction solver failed to converge after "
                       "the maximum {} iterations. Residual at last "
                       "interation was {}").format(self._max_loop, ind_err)
            module_logger.warning(log_msg)
        
        for i in range(self._turbine_count):
            
            turb = 'turbine{:0{width}d}'.format(i, width=n_digits)
            
            self._array.velHub[turb][:] = newVel[:, i]
            self._array.features[turb]['TIH'] = newTI[i]
        
        self.coefficient = newSpeed / iniSpeed
        self.inducedTKE = newTKE
        
        return
    
    # Simple formula for wake expansion while waiting for further development...
    # see dtocean_tidal/submodel/ParametricWake/wakeClass.py
    # This feature is not needed for DTOcean as is.
    def solv_wake(self, debug=False, debug_plot=False):
        """
        Computes turbine wakes and overlapping wake polygons

        Kwargs:
          debug (bool): debug flag
          debug_plot (bool): debug plot flag

        """ 
        debug = self._debug or debug 
        debug_plot = debug_plot or self._debug_plot
        
        n_digits = len(str(self._turbine_count))
        
        # Compute wake shape
        for i in range(self._turbine_count):
            turb = 'turbine{:0{width}d}'.format(i, width=n_digits)
            self._wakeShape[turb] = WakeShape(self._array.velHub[turb][:],
                                              self._array.streamlines[turb][:],
                                              self._wake[turb],
                                              self._bounding_box,
                                              debug=debug)#,
                                             #debug_plot=debug_plot)
        
        if debug_plot:
            
            fig = plt.figure(figsize=(18,10))
            ax = fig.add_subplot(111)
            
            for i in range(self._turbine_count):
                turb = 'turbine{:0{width}d}'.format(i, width=n_digits)
                x, y = self._wakeShape[turb].polygon.exterior.xy
                patch = PolygonPatch(self._wakeShape[turb].polygon,
                        alpha=0.1, zorder=2)
                ax.plot(x, y, color='#999999', alpha=0.1, zorder=1)
                ax.add_patch(patch)
                ax.set_aspect('equal')
                ax.set_ylabel('Distance (m)', fontsize = 12)
                ax.set_xlabel('Distance (m)', fontsize = 12)
            
            plt.show()
        
        return


def _solve_flow(turbine_count,
                turb_distances,
                turb_wakes,
                turb_velocity,
                turb_TI,
                turb_TKE,
                base_velocity,
                base_TKE,
                debug=False):
    
    new_vel = np.empty((2, turbine_count))
    new_TKE = np.empty(turbine_count)
    new_TI = np.empty(turbine_count)
    new_speed = np.empty(turbine_count)
    wake_mat = np.empty((turbine_count, turbine_count))
    tke_mat = np.empty((turbine_count, turbine_count))
    
    n_digits = len(str(turbine_count))
    
    for i in range(turbine_count):
        
        turb = 'turbine{:0{width}d}'.format(i, width=n_digits)
        
        for j in range(turbine_count):
            
            if j in turb_distances[turb].keys():
                
                (wake_mat[i, j],
                 tke_mat[i, j]) = turb_wakes[turb].get_velocity_TKE(
                                                 turb_distances[turb][j][:],
                                                 turb_velocity[:, i],
                                                 turb_TI[i],
                                                 debug=debug)
            
            else:
                
                wake_mat[i, j] = np.sqrt(turb_velocity[0, i] ** 2 +
                                                     turb_velocity[1, i] ** 2)
                tke_mat[i, j] = np.nan
    
    superposition_model = DominantWake(turb_velocity,
                                       wake_mat)
    coefficients = superposition_model.coefficients
    
    TKE_coeff_mat = get_wake_coefficients(turb_TKE,
                                          tke_mat)
    TKE_coeff = superposition_model.get_dominant(TKE_coeff_mat)
    new_TKE = TKE_coeff * base_TKE
    
    # Replace any nan values with TKE of turbines
    if np.isnan(new_TKE).any():
        new_TKE = np.where(np.isnan(new_TKE), turb_TKE, new_TKE)
    
    for i in range(turbine_count):
        
        new_vel[:,i] = base_velocity[:, i] * coefficients[i]
        new_speed[i] = np.sqrt(new_vel[0, i] ** 2 + new_vel[1, i] ** 2)
        new_TI[i] = _get_ti(new_TKE[i], new_speed[i])
    
    return new_vel, new_speed, new_TI, new_TKE


def _get_tke(I, U):
    
    return 1.5 * (I * U) ** 2.0


def _get_ti(k, U):
    
    two_thirds = 2. / 3
    
    return np.sqrt(two_thirds * k) / U
