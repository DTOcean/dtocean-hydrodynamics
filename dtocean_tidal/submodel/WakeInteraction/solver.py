#!/usr/bin/python2.7
# encoding: utf-8
from __future__ import division

# Start logging
import logging
module_logger = logging.getLogger(__name__)

import time
import numpy as np
import matplotlib.pyplot as plt
import os
from descartes import PolygonPatch

# Local import
from dtocean_tidal.submodel.ParametricWake import WakeShape, Wake
from dtocean_tidal.modules.blockage_ratio import blockage_ratio


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
      indMat (numpy.array): matrix of induction factors, float, (nb. device, nb. device)
      tiMat (numpy.array): matrix of turbulence intensity, float, (nb. device, nb. device)
      induction (numpy.array): turbines' resultant induction factors, float, (nb. device)
      inducedTI (numpy.array): turbines' resultant turbulence intensities, float, (nb. device)
      wake (dict): collection of dtocean_tidal.submodel.ParametricWake.Wake class objects for each device
      wakeShape (dict):  collection of dtocean_tidal.submodel.ParametricWake.WakeShape class objects for each device

    """

    def __init__(self, hydro, array, cfd_data, debug=False, debug_plot=False):

        self._turbine_count = len(array.positions.keys())
        self._debug = debug
        self._debug_plot = debug_plot
        self._hydro = hydro
        self._bounding_box = hydro.bounding_box
        self._array = array
        self.indMat = np.ones((self._turbine_count, self._turbine_count))
        self.tkeMat = np.zeros((self._turbine_count, self._turbine_count))
        self.induction = np.zeros(self._turbine_count)
        self.inducedTKE = np.zeros(self._turbine_count)
        self.wake = {}
        self.wakeShape = {}
        # compute global blockage ratio
        rbr = blockage_ratio(hydro, array, debug=debug)
        self.blockage = rbr * hydro.BR
        # Load the dataframe.
        self._df = cfd_data

        return

    def solv_induction(self, speed_superimpo='renkema', tke_superimpo='max', debug=False):
        """
        Compute the induction factor at each turbine hub
        based on each others interactions

        Kwargs:
          speed_superimpo (str): wake speed superimposition method
                                 it can be 'linear', 'geometric', 'rss', 'linrss', 'max', 'mean', 'renkema', 'prod', 'geosum'
          tke_superimpo (str): wake tke superimposition method
                                 it can be 'max', 'min', 'mean', 'sum', 'prod', 'rss'
          debug (bool): debug flag
          debug_plot (bool): debug plot flag
        """
        # Compute velocities at each turbine hub position
        debug = debug or self._debug

        # Initial speed and T.I. into two big matrix
        for i in range(self._turbine_count):
            p = self._array.velHub['turbine' + str(i)]
            speed = np.sqrt(p[0]**2.0 + p[1]**2.0)
            u = np.abs(p[0])
            v = np.abs(p[1])
            tke = 1.5 * (self._array.features['turbine' + str(i)]['TIH'] * speed)**2.0
            if i == 0:
                iniSpeed = speed
                iniTKE = tke
            else:
                iniSpeed = np.hstack((iniSpeed, speed))
                iniTKE = np.hstack((iniTKE, tke))

        # initialization Convergence criteria and other matrices
        newVel = np.zeros((2, self._turbine_count))
        newTKE = np.zeros(self._turbine_count)
        self.tkeMat[:] = np.resize(iniTKE, (self._turbine_count, self._turbine_count))

        if debug:
           module_logger.info("Querying parametric wake model...")
           start = time.time()
        for i in range(self._turbine_count):
            turb = 'turbine' + str(i)
            self.wake[turb] = Wake(self._df,
                                   self._array.features[turb],
                                   self.blockage,
                                   debug=self._debug)

            for turbNb in self._array.distances[turb].keys():
                self.indMat[i, int(turbNb)], self.tkeMat[i, int(turbNb)] = \
                self.wake[turb].ind_fac(self._array.distances[turb][turbNb][:],
                                        self._array.velHub[turb][:],
                                        self._array.features[turb]['TIH'],
                                        debug=self._debug)
        if debug:
            end = time.time()
            module_logger.info("...done after " + str(end - start) + " seconds.")
        # Cumulative induction: Sigma (1-U/Uinf)
        # TODO: perform validation against CFD to determine which method to choose
        speed_superimpo = speed_superimpo.lower()
        if not speed_superimpo in ['prod', 'linear', 'geometric', 'rss', 'linrss', 'max', 'renkema', 'geosum', 'mean']:
            module_logger.warning("Wrong superimposition method for momentum. Default value used")
            speed_superimpo = 'renkema'  # default method
        # ## Wake Interaction methods for intersection wakes cf. Palm 2011
        if speed_superimpo == 'linear':
            # Wake Inter. Method 1: linear superposition
            self.induction = np.sum(self.indMat,0)
            # filtering induction = 1.0 for summation
            indMat = self.indMat[:]
            indMat[self.indMat == 1.0] = 0.0
            # summation
            self.induction = np.sum(indMat, 0)
            # filtering back
            self.induction[np.where(self.induction == 0.0)] = 1.0
        elif speed_superimpo == 'geometric':
            # Wake Inter. Method 2: geometric superposition
            ratio = 1.0 - self.indMat
            # Filtering value = 0.0
            ratio[ratio == 0.0] = 1.0
            self.induction = 1.0 - np.prod(ratio,0)
            # filtering back
            self.induction[np.where(self.induction == 0.0)] = 1.0
        elif speed_superimpo == 'rss':
            # Wake Inter. Method 3: RSS of velocity deficits
            # filtering induction = 1.0 for summation
            indMat = self.indMat[:]
            indMat[self.indMat == 1.0] = 0.0
            # summation
            self.induction = np.sqrt(np.sum(indMat**2.0, 0))
            # filtering back
            self.induction[np.where(self.induction == 0.0)] = 1.0
        elif speed_superimpo == 'linrss':
            # Wake Inter. Method 4: average of RSS & linear superposition
            # filtering induction = 1.0 for summation
            indMat = self.indMat[:]
            indMat[self.indMat == 1.0] = 0.0
            # summation
            urss = np.sqrt(np.sum(indMat**2.0,0))
            uls = np.sum(indMat,0)
            self.induction = (urss + uls) / 2.0
            # filtering back
            self.induction[np.where(self.induction == 0.0)] = 1.0
        elif speed_superimpo == 'max':
            # Wake Inter. Method 5: Maximum wake deficit
            self.induction = self.indMat.max(axis=0)
        # ## Wake Superposition method cf. Renkema 2007
        elif speed_superimpo == 'renkema':
            # Wake Inter. Method 3: RSS of velocity deficits
            c = 1.0 - self.indMat
            self.induction = 1.0 - np.sqrt(np.sum(c**2.0, 0))
        # ## Wake Interaction methods for wake superimposition Thomas Roc
        elif speed_superimpo == 'prod':
            # Wake Inter. Method 1: simple multiplication of induction
            self.induction = np.prod(self.indMat, 0)
        elif speed_superimpo == 'geosum':
            # Wake Inter. Method 2: geometric superposition type
            ratio = 1.0 - self.indMat
            self.induction = 1.0 - np.sum(ratio, 0)
        elif speed_superimpo == 'mean':
            # Wake Inter. Method 6: Mean wake deficit# filtering induction = 1.0 for summation
            indMat = self.indMat[:]
            indMat[indMat == 1.0] = np.nan
            # summation
            self.induction = np.nanmean(self.indMat, axis=0)
            # filtering back
            self.induction[np.where(np.isnan(self.induction))] = 1.0

        tke_superimpo = tke_superimpo.lower()
        if not tke_superimpo in ['prod', 'min', 'rss', 'max', 'sum', 'mean']:
            module_logger.warning("Wrong superimposition method for tke. Default value used")
            tke_superimpo = 'max'  # default method
        # ##Wake Interaction methods for TKE superimposition Thomas Roc
        if tke_superimpo == 'max':
            # Wake Inter. Method 1: simple maximum of TKE
            self.inducedTKE = np.nanmax(self.tkeMat, axis=0)
        elif tke_superimpo == 'min':
            # Wake Inter. Method 2: simple average of TKE
            self.inducedTKE = np.nanmin(self.tkeMat, axis=0)
        elif tke_superimpo == 'mean':
            # Wake Inter. Method 2: simple average of TKE
            self.inducedTKE = np.nanmean(self.tkeMat, axis=0)
        elif tke_superimpo == 'sum':
            # Wake Inter. Method 3: simple sum of TKE
            self.inducedTKE = np.nansum(self.tkeMat, axis=0)
            self.inducedTKE[np.where(self.inducedTKE == 0.0)] = np.nan  # empty slice case
        elif tke_superimpo == 'prod':
            # Wake Inter. Method 4: simple prod of TKE
            self.inducedTKE = 1.0 - np.nanprod(1.0 + self.tkeMat, axis=0)
            self.inducedTKE[np.where(self.inducedTKE == 1.0)] = np.nan  # empty slice case
        elif tke_superimpo == 'rss':
            # Wake Inter. Method 3: RSS of tke
            self.inducedTKE = np.sqrt(np.nansum(self.tkeMat**2.0, 0))
            self.inducedTKE[np.where(self.inducedTKE == 0.0)] = np.nan  # empty slice case

        # iteratively compute and re-assign velocity and T.I. at hub
        for i in range(self._turbine_count):
            turb = 'turbine' + str(i)
            newVel[:,i] = self._array.velHub[turb][:] * self.induction[i]
            self._array.velHub[turb][:] = newVel[:, i]
            speed = np.sqrt(newVel[0, i]**2.0 + newVel[1, i]**2.0)
            if not np.isnan(self.inducedTKE[i]):
                newTKE[i] = self.inducedTKE[i]
            else:
                newTKE[i] = 1.5 * (self._array.features[turb]['TIH'] * speed)**2.0
            self._array.features[turb]['TIH'] = np.sqrt((2.0/3.0) * newTKE[i]) / speed

        # Recompute induction factor based on result speed / initial speed
        for i in range(self._turbine_count):
            p = self._array.velHub['turbine' + str(i)]
            speed = np.sqrt(p[0]**2.0 + p[1]**2.0)
            if i == 0:
                resSpeed = speed
            else:
                resSpeed = np.hstack((resSpeed, speed))
        self.induction = resSpeed / iniSpeed
         
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
            self.wakeShape[turb] = WakeShape(self._array.velHub[turb][:],
                                             self._array.streamlines[turb][:],
                                             self.wake[turb],
                                             self._bounding_box,
                                             debug=debug)#,
                                             #debug_plot=debug_plot)
        if debug_plot:
            fig = plt.figure(figsize=(18,10))
            ax = fig.add_subplot(111)
            for i in range(self._turbine_count):
                turb = 'turbine' + str(i)
                x, y = self.wakeShape[turb].polygon.exterior.xy
                patch = PolygonPatch(self.wakeShape[turb].polygon,
                        alpha=0.1, zorder=2)
                ax.plot(x, y, color='#999999', alpha=0.1, zorder=1)
                ax.add_patch(patch)
                ax.set_aspect('equal')
                ax.set_ylabel('Distance (m)', fontsize = 12)
                ax.set_xlabel('Distance (m)', fontsize = 12)
            plt.show()