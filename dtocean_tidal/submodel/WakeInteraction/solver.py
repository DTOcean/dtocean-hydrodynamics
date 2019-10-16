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
from descartes import PolygonPatch

# Local import
from .models import DominatingWake
from ..ParametricWake import WakeShape, Wake
from ...modules.blockage_ratio import blockage_ratio
from ...utils.misc import MovingAverage

# Start logging
module_logger = logging.getLogger(__name__)


class WakeInteraction:
    """
    -WakeInteraction class-
    
    Computes wake interaction and their impacts of flow field
    
    Args:
      hydro (dtocean_tidal.main.Hydro): Hydro class/object
      array (dtocean_tidal.main.Array): Array class/object
    
    Kwargs:
      debug (bool): debug flag
      debug_plot (bool): debug plot flag
    
    Attributes:
      indMat (numpy.array): matrix of induction factors,
                            float, (nb. device, nb. device)
      tiMat (numpy.array): matrix of turbulence intensity,
                           float, (nb. device, nb. device)
      induction (numpy.array): turbines' resultant induction factors,
                               float, (nb. device)
      inducedTI (numpy.array): turbines' resultant turbulence intensities,
                               float, (nb. device)
      wake (dict): collection of dtocean_tidal.submodel.ParametricWake.Wake
                   class objects for each device
      wakeShape (dict):  collection of
                         dtocean_tidal.submodel.ParametricWake.WakeShape class
                         objects for each device
    
    """
    
    def __init__(self, hydro,
                       array,
                       cfd_data,
                       criterior=1e-4,
                       max_loop=6,
                       moving_average_length=3,
                       debug=False,
                       debug_plot=False):
        
        self._turbine_count = len(array.positions.keys())
        self._debug = debug
        self._debug_plot = debug_plot
        self._hydro = hydro
        self._bounding_box = hydro.bounding_box
        self._array = array
        self._df = cfd_data
        
        nan_array = np.empty(self._turbine_count)
        nan_array[:] = np.nan
        
        self.induction = np.ones(self._turbine_count)
        self.inducedTKE = nan_array
        self._wake = {}
        self._wakeShape = {}
        
        rbr = blockage_ratio(hydro, array, debug=debug)
        self._blockage = rbr * hydro.BR
        
        self._criterior = criterior
        self._max_loop = max_loop
        self._average_induction = MovingAverage(moving_average_length)
        self._average_TKE = MovingAverage(moving_average_length)
        
        return
    
    def solv_induction(self, debug=False):
        
        """
        Compute the induction factor at each turbine hub
        based on each others interactions
        
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
        
        for i in range(self._turbine_count):
            
            turb = 'turbine' + str(i)
            
            p = self._array.velHub['turbine' + str(i)]
            
            iniSpeed[i] = np.sqrt(p[0]**2.0 + p[1]**2.0)
            iniVel[:, i] = self._array.velHub[turb][:]
            iniTI[i] = self._array.features[turb]['TIH']
            iniTKE[i] = _get_tke(iniTI[i], iniSpeed[i])
            
            self._wake[turb] = Wake(self._df,
                                    self._array.features[turb],
                                    self._blockage,
                                    debug=self._debug)
        
        if debug:
            module_logger.info("Querying parametric wake model...")
            start = time.time()
        
        ind_err = np.inf
        old_induction = 2
        newVel = iniVel.copy()
        newTI = iniTI.copy()
        newTKE = iniTKE.copy()
        loop_counter = 0
        
        while (abs(ind_err) > self._criterior and
               loop_counter < self._max_loop):
            
            (newVel,
             newSpeed,
             newTI,
             newTKE,
             last_induction) = _solve_induction(self._turbine_count,
                                                self._array.distances,
                                                self._wake,
                                                newVel,
                                                newTI,
                                                newTKE,
                                                iniVel,
                                                iniTKE,
                                                self._average_induction,
                                                self._average_TKE,
                                                debug)
            
            new_induction = np.mean(newSpeed / iniSpeed)
            ind_err = old_induction - new_induction
            old_induction = new_induction
            meanTI = np.mean(newTI)
            
            debug_msg = ("loop: {}; induction: {}; error: {}; TI: "
                         "{}").format(loop_counter,
                                      old_induction,
                                      abs(ind_err),
                                      meanTI)
            module_logger.debug(debug_msg)
            
            loop_counter += 1
        
        if debug:
            end = time.time()
            module_logger.info("...done after " + str(end - start) +
                                                               " seconds.")
        
        for i in range(self._turbine_count):
            
            turb = 'turbine' + str(i)
            
            self._array.velHub[turb][:] = newVel[:, i]
            self._array.features[turb]['TIH'] = newTI[i]
        
        self.induction = newSpeed / iniSpeed
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

        # Compute wake shape
        for i in range(self._turbine_count):
            turb = 'turbine' + str(i)        
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
                turb = 'turbine' + str(i)
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


def _solve_induction(turbine_count,
                     turb_distances,
                     turb_wakes,
                     turb_velocity,
                     turb_TI,
                     turb_TKE,
                     base_velocity,
                     base_TKE,
                     average_induction=None,
                     average_TKE=None,
                     debug=False):
    
    new_vel = np.empty((2, turbine_count))
    new_TKE = np.empty(turbine_count)
    new_TI = np.empty(turbine_count)
    new_speed = np.empty(turbine_count)
    wake_mat = np.empty((turbine_count, turbine_count))
    tke_mat = np.empty((turbine_count, turbine_count))
    
    for i in range(turbine_count):
        
        turb = 'turbine' + str(i)
        
        for j in range(turbine_count):
            
            if j in turb_distances[turb].keys():
                        
                (wake_mat[i, j],
                 tke_mat[i, j]) = turb_wakes[turb].induction(
                                                 turb_distances[turb][j][:],
                                                 turb_velocity[:, j],
                                                 turb_TI[j],
                                                 debug=debug)
                
            else:
                
                wake_mat[i, j] = np.sqrt(turb_velocity[0, j] ** 2 +
                                                     turb_velocity[1, j] ** 2)
                tke_mat[i, j] = turb_TI[j]
    
    superposition_model = DominatingWake(turb_velocity,
                                         wake_mat)
    induction = superposition_model.induction
    
    if average_induction is not None:
        induction = average_induction(induction)
    
    wake_TKE = tke_mat[superposition_model.indexes, range(tke_mat.shape[1])]
    TKE_perturbation = wake_TKE / turb_TKE[superposition_model.indexes]
    induced_TKE = base_TKE * TKE_perturbation
    
    # Replace any nan values with original TKE
    check_nan = np.isnan(induced_TKE)
    
    if check_nan.any():
        nan_idx = np.argwhere(check_nan)
        induced_TKE[nan_idx] = base_TKE[nan_idx]
    
    if average_TKE is not None:
        induced_TKE = average_TKE(induced_TKE)
    
    for i in range(turbine_count):
        
        new_vel[:,i] = base_velocity[:, i] * induction[i]
        new_speed[i] = np.sqrt(new_vel[0, i] ** 2 + new_vel[1, i] ** 2)
        
        new_TKE[i] = induced_TKE[i]
        new_TI[i] = _get_ti(new_TKE[i], new_speed[i])
    
    return new_vel, new_speed, new_TI, new_TKE, induction


def _get_tke(I, U):
    
    return 1.5 * (I * U) ** 2.0


def _get_ti(k, U):
    
    two_thirds = 2. / 3
    
    return np.sqrt(two_thirds * k) / U
