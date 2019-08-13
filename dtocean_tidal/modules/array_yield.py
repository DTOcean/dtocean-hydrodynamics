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

import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import interp1d

from dtocean_tidal.utils.misc import natural_sort

# Start logging
module_logger = logging.getLogger(__name__)


class ArrayYield:
    """
    Array yield class
    
    Args:
      Array (dtocean_tidal.main.Array): array class
    
    Kwargs:
      debug (bool): debug flag
      debug_plot (bool): debug plot flag
    
    Attributes:
      turbine_capacity (dict):
          Gathers wake induced turbines' capacities, Watts. Dictionary with
          turbine's ID as key
      turbine_capacity_no_interaction (dict):
          Gathers turbines' capacities, Watts. Dictionary with turbine's ID as
          key
      array_capacity (float): array's capacity, Watts
      array_capacity_no_interaction (float): array's capacity, Watts
    
    Note:
       wake interactions are not accounted here
    """
    def __init__(self, array, debug=False, debug_plot=False):
        
        self._debug = debug
        self._debug_plot = debug_plot
        self._array = array
        self.turbine_capacity = {}
        self.turbine_capacity_no_interaction = {}
        self.array_capacity = 0.0
        self.array_capacity_no_interaction = 0.0
        self._turbine_count = len(array.positions.keys())
        
        return
    
    def performance(self, debug=False, debug_plot=False):
        """
        Method computes both array and device performance power capacity,
        Watts.
        
        Kwargs:
          debug (bool): debug flag
          debug_plot (bool): debug plot flag
        """
        
        debug = debug or self._debug
        debug_plot = debug_plot or self._debug_plot
        if debug: module_logger.info("...Computing capacity & performance")
        
        # Water density
        rho = 1025.0
        
        l = natural_sort(self._array.velHub.keys())
        turbNorm = np.zeros(self._turbine_count)
        turbNormIni = np.zeros(self._turbine_count)
        turbGene = np.zeros(self._turbine_count)
        turbGeneIni = np.zeros(self._turbine_count)
        turbID = []
        
        for i, key in enumerate(l):
            
            # Turbine features
            diam = self._array.features[key]['Diam']
            cutIn = self._array.features[key]['cutIO'][0]
            cutOut = self._array.features[key]['cutIO'][1]
            ry = self._array.features[key]['RY']  # angle of attack
            rating = self._array.features[key]['Rating']
            
            
            # Flow speed at hub 
            u = self._array.velHub[key][0]
            v = self._array.velHub[key][1]
            uIni = self._array.velHubIni[key][0]
            vIni = self._array.velHubIni[key][1]
            norm = np.sqrt((u**2.0) + (v**2.0))
            normIni = np.sqrt((uIni**2.0) + (vIni**2.0))
            
            # Cp curve
            cp = interp1d(self._array.features[key]['Cp'][0],
                          self._array.features[key]['Cp'][1],
                          bounds_error=False, fill_value=0.0)
            
            # Capacity & perf
            Cp = cp(norm)
            Cpini = cp(normIni)
            
            power = _get_turbine_power(diam,
                                       cutIn,
                                       cutOut,
                                       rating,
                                       ry,
                                       Cp,
                                       norm,
                                       rho=rho)
            
            powerIni = _get_turbine_power(diam,
                                          cutIn,
                                          cutOut,
                                          rating,
                                          ry,
                                          Cpini,
                                          normIni,
                                          rho=rho)
            
            turbNorm[i] = norm
            turbNormIni[i] = normIni
            turbGene[i] = power
            turbGeneIni[i] = powerIni
            turbID.append(key)
            
            self.turbine_capacity[key] = power
            self.turbine_capacity_no_interaction[key] = powerIni
        
        totGene = np.nansum(turbGene)
        totGeneIni = np.nansum(turbGeneIni)
        
        self.array_capacity = totGene
        self.array_capacity_no_interaction = totGeneIni
        
        if debug_plot:
            
            fig1 = plt.figure(figsize=(18,10))
            ax1 = fig1.add_subplot(1,1,1)
            x = np.arange(self._turbine_count)
            y = turbGene / 1e6
            ax1.scatter(x,y, s=20)
            plt.xticks(x, turbID, rotation=90)
            # ax.set_xlabel("Turbine's ID", fontsize = 12)
            ax1.set_ylabel('Power generation (MWatts)', fontsize = 12)
            ax1.grid(True)
            text = "{:,.2f}".format(totGene/1e6) + "MWatts"
            text2 = "{:,.2f}".format(totGeneIni/1e6) + "MWatts"
            fig1.text(0.1,
                      0.95,
                      "Array production: " + text,
                      color='k',
                      fontsize=14)
            fig1.text(0.1,
                      0.925,
                      "Array production without interaction: " + text2,
                      color='k',
                      fontsize=14)
            plt.show()
            
            fig2 = plt.figure(figsize=(18,10))
            ax2 = fig2.add_subplot(1,1,1)
            x = np.arange(self._turbine_count)
            ax2.scatter(x, turbNorm, s=20)
            ax2.scatter(x, turbNormIni, s=20)
            plt.xticks(x, turbID, rotation=90)
            ax1.set_ylabel('Normal Velocity (m/s)', fontsize = 12)
            ax1.grid(True)
            plt.show()
        
        return
    
    def dissipated_mass_flow_rate(self, debug=False):
        """
        Computes the equivalent mass flow rate (kg/s) dissipated by turbine.
        
        Kwargs:
          debug (bool): debug flag
          
        Returns:
          float: dissipated mass flow rate
        """
        
        debug = debug or self._debug
        if debug: module_logger.info("Computing dissipated power...")
        #Constants
        rho = 1025.0 #water density
        #Initialize some stuffs
        l=natural_sort(self._array.velHub.keys())
        turbDiss = np.zeros(self._turbine_count)
        
        for i, key in enumerate(l):
            # turbine features
            diam = self._array.features[key]['Diam']
            cutIn = self._array.features[key]['cutIO'][0]
            cutOut = self._array.features[key]['cutIO'][1]
            ry = self._array.features[key]['RY']  # angle of attack
            A = np.pi * ((diam/2.0)**2.0) * np.cos(np.radians(ry))  # ellipse
            # Flow speed at hub
            u = self._array.velHub[key][0]
            v = self._array.velHub[key][1]
            norm = np.sqrt((u**2.0) + (v**2.0))
            #Ct curve
            ct = interp1d(self._array.features[key]['Ct'][0],
                          self._array.features[key]['Ct'][1],
                          bounds_error=False, fill_value=0.0)
            #capacity & perf.
            Ct = ct(norm)
            if norm < cutIn:
                disspow = 0.0
            elif norm > cutOut:
                disspow = 0.0
            else:
                disspow = rho*Ct*A*norm
            
            turbDiss[i] = disspow
        
        return np.nansum(turbDiss)


def _get_turbine_power(diameter,
                       cut_in,
                       cut_out,
                       rating,
                       angle_of_attack,
                       Cp,
                       u,
                       rho=1025.):
    
    if (np.abs(angle_of_attack) >= 90 or
        u < cut_in or
        u > cut_out): return 0.
    
    aoa_rad = np.radians(np.abs(angle_of_attack))
    power = np.pi / 8. * rho * diameter ** 2 * np.cos(aoa_rad) * Cp * u ** 3.0
    
    if power > rating:
        power = rating
    elif power < 0.:
        power = 0.
    
    return power
