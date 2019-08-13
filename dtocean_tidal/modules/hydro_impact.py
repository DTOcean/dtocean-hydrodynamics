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
import matplotlib.pyplot as plt

# local imports
from dtocean_tidal.utils.misc import transec_surf

class HydroImpact:
    """
    Hydrodynamic impact class

    Args:
      array (dtocean_tidal.main.Array): Array class/object
      hydro (dtocean_tidal.main.Hydro): Hydro class/object
      interaction (dtocean_tidal.main.WakeInteraction): WakeInteraction class/object
      arrayyield (dtocean_tidal.modules.ArrayYield): ArrayYield class/object

    Attributes:
      diss_avai_mass_flow_rate (float): the ratio between dissipated and available mass flow rate
      from the flow for the considered bin, [0 to 1]
      u_reduced (numpy.array): turbine induced u velocity component/field, 2D numpy array
      v_reduced (numpy.array): turbine induced v velocity component/field, 2D numpy array

    """
    def __init__(self, array, hydro, interaction, arrayyield, debug=False, debug_plot=False):
        self._debug = debug
        self._debug_plot = debug_plot
        self._array = array
        self._hydro = hydro
        self._arrayyield = arrayyield
        self._interaction = interaction
        self.diss_avai_mass_flow_rate = 0.0
        self.u_reduced = np.zeros(self._hydro.U.shape)
        self.v_reduced = np.zeros(self._hydro.V.shape)
        self._turbine_count = len(array.positions.keys())

    # This function is still under development as it is not explicitly need
    # for the DTOcean suite of tools
    def reduced_velocity_field(self, debug=False, debug_plot=False):
        """
        Computes reduced velocity field

        Kwargs:
          debug (bool): debug flag
          debug_plot (bool): debug plot flag

        """
        debug = debug or self._debug
        debug_plot = debug_plot or self._debug_plot
        #TR comment: simplied version for information flow sake, alpha-version
        #Intial velocity field
        u = self._hydro.U[:]
        v = self._hydro.V[:]
        x = self._hydro.X[:]
        y = self._hydro.Y[:]

        if debug:
            module_logger.info('...find mesh elements within 5 diam circle..')
        
        n_digits = len(str(self._turbine_count))
        
        for i in range(self._turbine_count):
            ID = 'turbine{:0{width}d}'.format(i, width=n_digits)
            xt = self._array.positions[ID][0]
            yt = self._array.positions[ID][1]
            diam = self._array.features[ID]['Diam']
            # in case only one turbine
            if self._interaction == None:
                coefficient = 0.0
            else:
                coefficient = self._interaction.coefficient[i]
            #find mesh elements within 5 diam circle
            xx,yy=np.meshgrid(x,y)
            distance=np.sqrt((xx-xt)**2.0 + (yy-yt)**2.0)
            #build mask proportional to distance
            mask = distance/(2.5*diam)#reduces flow only around turbine
            index = np.where(mask<1.0)
            mask = 1.0-mask#)**2.0#cubic distance mask
            #apply mask
            #TR quick fix for area with surging
            if (1-coefficient) < 0.0: coefficient = 1.0
            u[index]=u[index]-(mask[index]*(1-coefficient)*u[index])
            v[index]=v[index]-(mask[index]*(1-coefficient)*v[index])

        #load into attributs
        self.u_reduced[:] = u[:]
        self.v_reduced[:] = v[:]

        #plotting
        if debug_plot:
            fig = plt.figure(figsize=(18,10))
            ax = fig.add_subplot(1,1,1)
            norm = np.sqrt((self.u_reduced)**2.0 + (self.v_reduced)**2.0)
            mappy = ax.contourf(xx,yy,norm,100)
            fig.colorbar(mappy)
            ax.set_ylabel('distance (m)', fontsize = 12)
            ax.set_xlabel('distance (m)', fontsize = 12)
            fig.suptitle('Velocity field impact (m/s)')
            ax.axis('tight')
            ax.set_aspect('equal')
            ax.grid(True)

    def ini_mass_flow_rate(self, debug=False):
        """
        Computes initial mass flow rate (kg/s)

        Kwargs:
          debug (bool): debug flag

        Returns:
          ini_mfr (float): initial mass flow rate

        """
        if debug or self._debug: module_logger.info("Computes initial momentum flux...")
        transect, first_row, speed = transec_surf(self._hydro, self._array, debug=debug)
        rho = 1025.0
        coeff = 1.0 / self._hydro.BR  # assuming the mass flow rate of the entire site
        ini_mfr = rho * transect * speed * coeff

        return ini_mfr
            
    def dissipated_over_available_flux_ratio(self, debug=False):
        """
        Computes the ratio between dissipated and available mass flow rate
        from the flow for the considered bin, float, [0 to 1]

        Kwargs:
          debug (bool): debug flag

        Returns:
          float: ratio between dissipated and available mass flow rate

        """
        if debug or self._debug: module_logger.info("Computes extracted-over-available energy ratio...")
        ini_mfr = self.ini_mass_flow_rate(debug=debug)
        diss_mfr = self._arrayyield.dissipated_mass_flow_rate(debug=debug)

        # load to attribute
        self.diss_avai_mass_flow_rate = diss_mfr / ini_mfr