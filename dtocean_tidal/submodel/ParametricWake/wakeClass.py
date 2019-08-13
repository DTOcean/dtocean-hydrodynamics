# -*- coding: utf-8 -*-

#    Copyright (C) 2016 Thomas Roc, Mathew Topper
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
import matplotlib.pyplot as plt
from shapely.geometry import Polygon, Point, LineString
from shapely.geos import TopologicalError
from descartes import PolygonPatch
from scipy.interpolate import interp1d

# Local import
from dtocean_tidal.utils.misc import line, intersection, closest_point
from dtocean_tidal.utils.interpolation import interp_at_point

# Start logging
module_logger = logging.getLogger(__name__)


class Wake:
    """
    Wake class computes velocities and TKEs for any distance (x,y) and any set
    of turbine parameters
    
    Args:
      data_reader (fortran): database reader
      turbParams (dict): turbine's parameters
      BR (float): blockage ratio, float
    
    """
    
    def __init__(self, dataframe,
                       U_dict,
                       V_dict,
                       TKE_dict,
                       turbParams,
                       BR,
                       debug=False):
        
        self._debug = debug
        
        # ...blockage ratio
        if BR < 0.0:
            self.BR = 0.0
        elif BR > 0.7:
            self.BR = 0.7
            module_logger.debug("Blockage ratio > 70%... reaching model's "
                                "theoretical limit.")
        else:
            self.BR = BR
        
        # ...relative yawing angle
        ry = np.radians(turbParams['RY'])
        
        if ry < 0.0:
            self.RY = ry * -1.0
        elif ry > np.radians(89.0):
            module_logger.debug("Relative angle with inflow > 90.0 deg... "
                                "reaching model's theoretical limit.")
            self.RY = ry
        else:
            self.RY = ry
            
        self.TIH = turbParams['TIH']
        self.Ct = interp1d(turbParams['Ct'][0],
                           turbParams['Ct'][1],
                           bounds_error=False,
                           fill_value=0.0)
        self.Diam = turbParams['Diam']
        self.cutIn = turbParams['cutIO'][0]
        self.cutOut = turbParams['cutIO'][1]
        
        # database values
        self._dfU = U_dict
        self._dfV = V_dict
        self._dfTKE = TKE_dict
        self._cts = np.copy(dataframe['cts'][:])
        self._tis = np.copy(dataframe['tis'][:])
        self._cts.sort()
        self._tis.sort()
        
        # Real distance Empirical relationship for vertical blockage ratio and
        # yawing
        ratio = 0.489*self.BR**2.0 - 1.472*self.BR + 1.0
        self._dfx = dataframe['dfX'] * ratio * self.Diam
        self._dfy = dataframe['dfY'] * abs(np.cos(ry)) * self.Diam
        
        return
    
    def read_at_point(self, x, y, ct, ti, debug=False):
        """
        Interpolates CFD datasets based on Ct and TI values
        Args:
          x (float): relative distance along streamline, m
          y (float): relative distance across streamline, m
          ct (float): turbine's thrust coefficient, dimensionless
          ti (float): turbulence intensity at hub's location,
                      dimensionless [0 ; 1]
        
        Returns:
          u (numpy array): 2D array, flow field, m/s
          tke (numpy array): 2D array, TKE filed, m2/s2
        """
        
        debug = debug or self._debug
        if debug: module_logger.info("Picking CFD datasets...")
        
        # Looking through available Ct's and TI's values
        ctbounds = []
        
        for ss, ee in zip(self._cts[:-1], self._cts[1:]):
            if ss <= ct <= ee:
                ctbounds = [ss, ee]
                break
        
        tibounds = []
        
        for ss, ee in zip(self._tis[:-1], self._tis[1:]):
            if ss <= ti <= ee:
                tibounds = [ss, ee]
                break
        
        # If requested values outside of range
        if ctbounds == []:
            
            module_logger.debug(("Requested Ct value {} is not covered by the "
                                 "current CFD database.").format(ct))
            
            idx = (np.abs(self._cts - ct)).argmin()
            
            try:
                ctbounds = [self._cts[idx-1], self._cts[idx]]
            except IndexError:
                ctbounds = [self._cts[idx], self._cts[idx + 1]]
            
            module_logger.debug(("Ct value of {} will be "
                                 "used.").format(self._cts[idx]))
        
        if tibounds == []:
            
            module_logger.debug(("Requested TI value {} is not covered by "
                                 "the current CFD database.").format(ti))
            
            idx = (np.abs(self._tis - ti)).argmin()
            
            try:
                tibounds = [self._tis[idx-1], self._tis[idx]]
            except IndexError:
                tibounds = [self._tis[idx], self._tis[idx+1]]
            
            module_logger.debug(("TI value of {} will be "
                                 "used.").format(self._tis[0]))
        
        # compute distance
        dist = np.zeros(4)
        dist[0] = (ct - ctbounds[0])**2.0 + (ti - tibounds[0])**2.0
        dist[1] = (ct - ctbounds[0])**2.0 + (ti - tibounds[1])**2.0
        dist[2] = (ct - ctbounds[1])**2.0 + (ti - tibounds[0])**2.0
        dist[3] = (ct - ctbounds[1])**2.0 + (ti - tibounds[1])**2.0
        mloc = dist.argmax()
        
        # computing weights
        nn = 0
        atmp = np.zeros(3)
        btmp = np.zeros(3)
        
        for n in range(4):
            
            if not n == mloc:
                
                if n==0:
                    atmp[nn] = ctbounds[0]
                    btmp[nn] = tibounds[0]
                elif n==1:
                    atmp[nn] = ctbounds[0]
                    btmp[nn] = tibounds[1]
                elif n==2:
                    atmp[nn] = ctbounds[1]
                    btmp[nn] = tibounds[0]
                elif n==3:
                    atmp[nn] = ctbounds[1]
                    btmp[nn] = tibounds[1]
                
                nn += 1
        
        wght = np.ones(3)
        
        wght[0] = (btmp[1] - btmp[2]) * (ct - atmp[2]) + \
                                      (atmp[2] - atmp[1]) * (ti - btmp[2])
        wght[1] = (btmp[2] - btmp[0]) * (ct - atmp[2]) + \
                                      (atmp[0] - atmp[2]) * (ti - btmp[2])
        wght[0:1] = wght[0:1] / ((btmp[1] - btmp[2]) * (atmp[0] - atmp[2]) +
                                 (atmp[2] - atmp[1]) * (btmp[0] - btmp[2]))
        wght[2] = 1.0 - wght[0] - wght[1]
        
        X = self._dfx
        Y = self._dfy
        
        if debug: module_logger.info("...interpolation...")
        u = 0.0
        v = 0.0
        tke = 0.0
        
        for n in range(3):
            
            if wght[n] == 0: continue
            
            tmpCt = atmp[n]
            tmpTi = btmp[n]
            
            umat = self._dfU['ti' + str(tmpTi)]['ct' + str(tmpCt)] * wght[n]
            vmat = self._dfV['ti' + str(tmpTi)]['ct' + str(tmpCt)] * wght[n]
            tkemat = self._dfTKE['ti' + str(tmpTi)][
                                                'ct' + str(tmpCt)] * wght[n]
            
            tmp = interp_at_point(x, y, X, Y, [umat.T, vmat.T, tkemat.T])
            
            u += tmp[0]
            v += tmp[1]
            tke += tmp[2]
        
        return u, v, tke
    
    def get_velocity_TKE(self, distance, velHub, tiHub, debug=False):
        """Return velocity and T.K.E at location behind turbine for input
        velocity and T.I. conditions.
        
        Args:
          distance (numpy.array or list): along and across distances to hub
                                          axis (m)
          or 2D array (Nturb, 2) or list [x,y]
          velHub (list): velocity components at hub, [u, v], float list
          tiHub (float): turbulence intensity at hub (%)
        
        Kwargs:
          debug (bool): debug flag
        
        Returns:
          wake_speed (float): wake velocity magnitude
          newTI (float): new turbulence intensity
        
        """
        
        debug = debug or self._debug
        
        x = np.asarray(distance[0])
        y = np.asarray(distance[1])
        
        norm = np.sqrt((velHub[0]**2.0) + (velHub[1]**2.0))

        try:
            
            Ct = self.Ct(norm)
            
        except ValueError:
            
            if debug: module_logger.info("Flow speed value outside of Ct's "
                                         "curve range. Ct will be set to 0.0")
            
            Ct = 0.0
        
        ry = np.radians(self.RY)
        
        # Checking bounds
        if Ct > self._cts[-1]:
            
            Ct = self._cts[-1]
            module_logger.debug(("Thrust coefficient > {}... reaching "
                                 "parametric database's "
                                 "limit.").format(self._cts[-1]))
        
        if Ct < self._cts[0]:
            
            module_logger.debug(("Thrust coefficient < {}... reaching "
                                 "parametric database's "
                                 "limit.").format(self._cts[0]))
        
        if tiHub < self._tis[0]:
            
            tiHub = self._tis[0]
            module_logger.debug(("T.I. < {}... reaching "
                                 "parametric database's "
                                 "limit.").format(self._tis[0]))
        
        if tiHub > self._tis[-1]:
            
            tiHub = self._tis[-1]
            module_logger.debug(("T.I. > ... reaching "
                                 "parametric database's "
                                 "limit.").format(self._tis[-1]))
        
        X = self._dfx
        Y = self._dfy
        
        if (Ct > np.min(self._cts) and ry < np.radians(89.0)) and \
           (X.min() < x < X.max()) and (Y.min() < y < Y.max()):
            
            u, v, tke = self.read_at_point(x, y, Ct, tiHub, debug=debug)
            
            indFac = np.sqrt(u**2.0 + v**2.0)
            newTKE = norm * tke
        
        else:
            
            indFac = 1.0
            newTKE = np.nan
        
        wake_speed = norm * indFac
        
        return wake_speed, newTKE


# Simple formula for wake expansion while waiting for further development...
# This feature is not needed for DTOcean as is.
class WakeShape:
    """
    Wake class

    Computes the wake shape for any set of turbine parameters, ie:
      - VBR = vertical blockage ratio
      - TIH = turbulence intensity at hub (%)
      - YA = yaw angle (deg.)
      - Ct = thrust coefficient
      - A = rotor diameter (m2)
      - velHub = velocity components at hub, [u, v], float list

    Args:
      velHub (list): velocity components at hub, [u, v], float list
      streamline (dtocean_tidal.modules.streamline): streamline object
      wake (dtocean_tidal.modules.Wake): Wake class/object
      bounding_box (shapely.geometry.polygon.Polygon): domain's limits

    Kwargs:
      debug (bool): debug flag
      debug_plot (bool): debug plot flag

    Attributes:
      streamtop (list): top streamline's coordinates
      streambot (list): bottom streamline's coordinates

    """

    def __init__(self, velHub, streamline, wake, bounding_box,
                 debug=False, debug_plot=False):

        #Bounding box, shapely polygon
        self._bounding_box = bounding_box
        #initialise top and bottom wake stream boundaries
        self.streamtop = [[0]*len(streamline[0]),[0]*len(streamline[1])]
        self.streambot = [[0]*len(streamline[0]),[0]*len(streamline[1])]
        L = len(streamline[0])
        #Compute hypothenus
        hypo = np.zeros(L-1)
        X = np.asarray(streamline[0])
        Y = np.asarray(streamline[1])
        hypo = np.sqrt(((X[:-1] - X[1:])**2.0) + ((Y[:-1] - Y[1:])**2.0)) 
        #Compute the outer lines of the wake
        ##Flow signs
        deltaY = Y[1]-Y[0]
        sy = np.sign(deltaY)
        deltaX = X[1]-X[0]
        sx = np.sign(deltaX)
        ##Tilt first outer points
        l = wake.Diam/2.0 #wake expansion
        gamma = np.arctan2(deltaY,deltaX)           
        self.streamtop[0][0] = X[0] - (l*np.sin(gamma))
        self.streamtop[1][0] = Y[0] + (l*np.cos(gamma))
        self.streambot[0][0] = X[0] + (l*np.sin(gamma))
        self.streambot[1][0] = Y[0] - (l*np.cos(gamma))

        
        distance = 0.0
        for i in range(L-1):
            ##Flow signs
            deltaY = Y[i+1]-Y[i]
            sy = np.sign(deltaY)
            deltaX = X[i+1]-X[i]
            sx = np.sign(deltaX)

            distance += hypo[i] #cumulative distance along streamline
            l = (wake.Diam/2.0) + (wake.l * distance) #wake expansion
            AB = np.sqrt((deltaY**2.0) + (deltaX**2.0))
            AC = np.sqrt((AB**2.0) + (l**2.0))
            gamma = np.arctan((Y[i+1]-Y[i])/(X[i+1]-X[i]))
            alpha1 = np.arctan(l/AB)
            alpha2 = np.arctan(-l/AB)           
            beta = gamma + alpha1
            delta = gamma + alpha2 
            self.streamtop[0][i+1] = X[i] + ((AC*np.cos(beta))*sx)
            self.streamtop[1][i+1] = Y[i] + ((AC*np.sin(beta))*sx) 
            self.streambot[0][i+1] = X[i] + ((AC*np.cos(delta))*sx)
            self.streambot[1][i+1] = Y[i] + ((AC*np.sin(delta))*sx)

        #Make sure that the last points are outside or on the bounding box
        pt1 = Point(self.streambot[0][-1], self.streambot[1][-1])
        pt2 = Point(self.streamtop[0][-1], self.streamtop[1][-1])
        ptEnd = [streamline[0][-1], streamline[1][-1]]
        pts = list(self._bounding_box.exterior.coords)[:-1]
        L0 = line(pts[0], pts[1])
        L1 = line(pts[1], pts[2])
        L2 = line(pts[2], pts[3])
        L3 = line(pts[3], pts[0])
        Ls = [L0, L1, L2, L3]
        #Need to find if one or the two last points outside of box
        if (self._bounding_box.contains(pt1) or self._bounding_box.contains(pt2)):
            if self._bounding_box.contains(pt1):
                Lref = line([self.streambot[0][-2], self.streambot[1][-2]],
                            [self.streambot[0][-1], self.streambot[1][-1]])
                Lbis = line([self.streamtop[0][-2], self.streamtop[1][-2]],
                            [self.streamtop[0][-1], self.streamtop[1][-1]])
            if self._bounding_box.contains(pt2):
                Lref = line([self.streamtop[0][-2], self.streamtop[1][-2]],
                            [self.streamtop[0][-1], self.streamtop[1][-1]])
                Lbis = line([self.streambot[0][-2], self.streambot[1][-2]],
                            [self.streambot[0][-1], self.streambot[1][-1]])
            #if one point: find intersection with bounding box
            Rs = []
            Rbiss = []
            for L in Ls:
                R = intersection(L,Lref)
                Rbis = intersection(L,Lbis)
                if not R == False:
                    Rs.append(R)
                if not R == False:
                    Rbiss.append(Rbis)
            clPt, ind, dist = closest_point(ptEnd, Rs, debug=debug)
            if self._bounding_box.contains(pt1):
                self.streambot[0].append(clPt[0])
                self.streambot[1].append(clPt[1])
                if not Rbiss[ind] == False:
                    self.streamtop[0].append(Rbiss[ind][0])
                    self.streamtop[1].append(Rbiss[ind][1])                    
            if self._bounding_box.contains(pt2):
                self.streamtop[0].append(clPt[0])
                self.streamtop[1].append(clPt[1])
                if not Rbiss[ind] == False:
                    self.streambot[0].append(Rbiss[ind][0])
                    self.streambot[1].append(Rbiss[ind][1]) 
        else:
            #if two points: add closest boundingbox corner to wake
            c1 = [self.streambot[0][-1], self.streambot[1][-1]]
            c2 = [self.streamtop[0][-1], self.streamtop[1][-1]] 
            c3 = [self.streamtop[0][0], self.streamtop[1][0]]
            c4 = [self.streambot[0][0], self.streambot[1][0]]
            test_poly = Polygon([c1, c2, c3, c4])
            test_line = LineString([c1,c2])
            cPt = []
            for pt in pts:
                if test_line.crosses(self._bounding_box):
                    trgl = Polygon([c1, c2, pt])
                    test = trgl.intersection(test_poly)
                    if (test.area==0.0):                
                        self.streambot[0].append(pt[0])
                        self.streambot[1].append(pt[1])

        if debug_plot:
            fig = plt.figure(figsize=(18,10))
            ax1 = fig.add_subplot(121)
            ax1.plot(streamline[0], streamline[1])
            ax1.plot(self.streamtop[0], self.streamtop[1])
            ax1.plot(self.streambot[0], self.streambot[1])
            ax1.set_aspect('equal')
            #plt.show()

        #Define wake shapes with shapely
        extPt = []
        for i, j in zip(self.streambot[0][:], self.streambot[1][:]):
            extPt.append((i,j))
        for i, j in zip(self.streamtop[0][::-1], self.streamtop[1][::-1]):
            extPt.append((i,j))
        self.polygon = Polygon(extPt)

        #Cut-out wake bits outside of bounding box
        try:
            self.polygon = self.polygon.intersection(self._bounding_box)
        #In case added points mess up the shapely polygon
        except TopologicalError:
            extPt.pop(len(extPt)//2)
            self.polygon = Polygon(extPt)
            try:
                self.polygon = self.polygon.intersection(self._bounding_box)
            except TopologicalError:
                extPt.pop(len(extPt)//2)
                self.polygon = Polygon(extPt)
                self.polygon = self.polygon.intersection(self._bounding_box)

        #Control plots 
        if debug_plot:
            ax = fig.add_subplot(122)
            x, y = self.polygon.exterior.xy
            ax.plot(x, y, 'o', color='#999999', zorder=1)
            patch = PolygonPatch(self.polygon, alpha=0.5, zorder=2)
            ax.add_patch(patch)
            ax.set_aspect('equal')
            #plt.show()