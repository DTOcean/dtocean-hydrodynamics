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
from shapely.geometry import Polygon
from scipy.interpolate import RectBivariateSpline

# Local import
from .utils import distance_from_streamline
from .utils.interpolation import interp_at_point
from .utils.misc import bearing_to_radians, vector_to_bearing
from .modules.vertical_velocity_profile import vvpw
from .modules import Streamlines, ArrayYield, HydroImpact
from .submodel.WakeInteraction import WakeInteraction
from .submodel.ParametricWake import read_database

# Start logging
module_logger = logging.getLogger(__name__)


class Hydro:
    """
    Site data & Hydrodynamic conditions

    Args:
      data (dict): dictionary gathering all site's information

    kwargs:
      debug (bool): output additional info. during runtime
      debug_plot (bool): plot additional info. during runtime

    Attributes:
      U (numpy.array): depth-averaged u velocity component (West-East), 2d numpy array, m/s
      V (numpy.array): depth-averaged v velocity component (South-North), 2d numpy array, m/s
      TI (numpy.arra): turbulence intensity, 2d numpy array
      SSH (numpy.array): sea surface elevation, 2d numpy array, m
      bathy (numpy.array): bathymetry, 2d numpy array, m
      geophy (numpy.array): geophysic info, 2d numpy array, Manning coefficient or bed roughness if using Soulsby formulation
      BR (float): blockage ratio, defined as lease surface / site area, user input, float
      X (numpy.array): x data positions (West-East), assuming regular grid, 1d array, m
      Y (numpy.array): y data positions (South-North), assuming regular grid, 1d array, m
      dx (float): x spatial step, assuming regular grid, float, m
      dy (float): y spatial step, assuming regular grid, float, m
      beta (float): bed roughness
      alpha (float): power law exponent

    """
    def __init__(self, data, debug=False, debug_plot=False):
        if debug: module_logger.info("Loading Hydrodynamic conditions...")
        self._debug = debug
        self.data = data
        self.U = self.data['U']
        self.V = self.data['V']
        self.TI = self.data['TI']
        self.SSH = self.data['SSH']
        self.bathy = self.data['bathy']
        self.X = self.data['X']
        self.Y = self.data['Y']
        self.lease = Polygon(self.data['lease'])
        (xm, ym, xM, yM) = self.lease.bounds
        self.bounding_box = Polygon([[xm, ym], [xM, ym], [xM, yM], [xm, yM]])
        self.BR = self.data['BR'] # lease surface / site area, user input
        self.beta = self.data['beta']
        self.alpha = self.data['alpha']
        
        # Set up interpolators for U and V
        U0 = np.nan_to_num(self.U)
        V0 = np.nan_to_num(self.V)
        
        self.interpU = RectBivariateSpline(self.X, self.Y, U0.T)
        self.interpV = RectBivariateSpline(self.X, self.Y, V0.T)

        if debug_plot:
            plt.figure(figsize=(18,10))
            norm = np.sqrt(self.U**2.0 + self.V**2.0)
            #plt.quiver( self.X, self.Y, self.U, self.V, norm, alpha=.5)
            plt.contourf(self.X, self.Y, norm, alpha=.5)
            plt.colorbar()
            plt.quiver(self.X, self.Y, self.U, self.V, alpha=.5,
                        edgecolor='k', facecolor='None', linewidth=.5)
            plt.show()
        return


class Array:
    """
    Gathers site data & Hydrodynamic conditions, turbine coordinates & features
    Computes streamlines, vertical profiles & yaw angles at hub heights

    Args:
      hydro(dtocean_tidal.main.Hydro) hydro class object
      turbines (dict): dictionnary gathering all turbines' positions
      features (dict): dictionnary gathering all turbines' features

    kwargs:
      debug (bool): output additional info. during runtime
      debug_plot (bool): plot additional info. during runtime

    Attributes:
      positions (dict): dictionary with turbine's ID as key, gathers turbine's coordinates
        ex: Array.positions[ID] = [x,y,z], 1d numpy array
      features (dict): dictionary gathering turbine's coordinates with the following keys:
          . HAS = heading angle span (deg.), float
          . RY = yaw angle relative to incoming flow (deg.), float
          . OA = orientation angle (deg.), float
          . Cp = power curve, velocity vs Cp, (2,n) array, [vel., Cp]
          . Ct = thrust curve, velocity vs Ct, (2,n) array, [vel., Ct]
          . Diam = turbine diameter, m, float
          . cutIO = cut-in and cut-out speed, m/s, [float,float]
          . floating = floating or not, boolean
          . 2way = working both ways or not, boolean
      streamlines (dict): dictionary with turbine's ID as key, gathers turbine's streamline
        ex: Array.streamlines[ID] = [X,Y] where X and Y are lists of coordinates
      distances (dict): dictionary with turbine's ID as key, gathers turbine's Relative
        distances from streamline to one another
      velHub (dict): dictionary with turbine's ID as key, gathers turbine's flow speed at hub

    """
    def __init__(self, hydro,
                       turbines,
                       features,
                       debug=False,
                       debug_plot=False):
        
        turbine_count = len(turbines.keys())
        n_digits = len(str(turbine_count))
        
        self.turbine_count = turbine_count
        
        #Turn turbine positions into one big matrix
        self.positions = {}
        
        for i in range(turbine_count):
            
            turb_name = 'turbine{:0{width}d}'.format(i, width=n_digits)
            
            p = turbines[turb_name]['position']
            
            if i == 0:
                positions = p
            else:
                positions = np.vstack((positions, p))
            
            self.positions[turb_name] = p
        
        self.features = features
        
        # Streamlines
        data = {}
        xn = hydro.X.shape[0]
        yn = hydro.Y.shape[0]
        
        xmin = hydro.X.min()
        xmax = hydro.X.max()
        ymin = hydro.Y.min()
        ymax = hydro.Y.max()
        dx = (xmax - xmin) / (xn - 1)
        dy = (ymax - ymin) / (yn - 1)
        xi = np.arange(xmin, xmax + dx, dx)
        yi = np.arange(ymin, ymax + dy, dy)
        
        data['X'] = xi
        data['Y'] = yi
        data['U'] = hydro.U
        data['V'] = hydro.V
        data['interpU'] = hydro.interpU
        data['interpV'] = hydro.interpV
        data['lease'] = hydro.lease
        
        # In case there is only one turbine in the array
        if not turbine_count == 1:
            
            #  computes streamlines
            turb_zero = 'turbine{:0{width}d}'.format(0, width=n_digits)
            
            diam = features[turb_zero]['Diam']
            max_len = 20 * diam
            
            SLs = Streamlines(data,
                              positions,
                              turbine_count,
                              maxlen=max_len,
                              debug=debug)
            
            if debug_plot: SLs.plot(turbine_count)
            
            #Relative distance from streamlines
            self.streamlines={}
            self.distances={}
            
            for i in range(turbine_count):
                
                if debug:
                    module_logger.info("Streamline nb: {}".format(i))
                
                turb_name = 'turbine{:0{width}d}'.format(i, width=n_digits)
                streamline = SLs.streamlines[i]
                
                self.streamlines[turb_name] = streamline
                self.distances[turb_name] = \
                            distance_from_streamline(streamline,
                                                     positions,
                                                     debug=debug,
                                                     debug_plot=debug_plot)
        
        else:
            
            if debug: module_logger.info("Single turbine...no streamlines")
        
        # Velocity and TI at hub
        self.velHub = {}
        self.velHubIni = {}  # save initial velocity for computing Q factor
        
        if debug:
            module_logger.info("Computing Hub velocities...")
        
        for i in range(turbine_count):
            
            if debug:
                module_logger.info("Interpolating quantities...")
            
            turb_name = 'turbine{:0{width}d}'.format(i, width=n_digits)
            
            Q = [hydro.U,
                 hydro.V,
                 hydro.SSH,
                 hydro.bathy,
                 hydro.TI]
            
            (u, v, el, h, ti) = interp_at_point(self.positions[turb_name][0],
                                                self.positions[turb_name][1],
                                                hydro.X,
                                                hydro.Y,
                                                Q)
            
            # Quantity formatting
            if h < 0.0: h *= -1.0
            
            # Get hub height accounting for floating options
            hh = self.positions[turb_name][2]
            if hh > 0.0: hh *= -1.0
            
            if self.features[turb_name]['floating']:
                z = (el + h) + hh
            else:
                z = -hh
            
            # Computing velocity vertical profile weight
            radius = self.features[turb_name]['Diam'] / 2.0
            z_top = z + radius
            z_bottom = z - radius
            args = [el,
                    h,
                    hydro.beta,
                    hydro.alpha,
                    debug]
            
            w_top = vvpw(z_top, *args)
            w_bottom = vvpw(z_bottom, *args)
            
            # Linear interpolation along the vertical
            w = (w_top + w_bottom) / 2.0
            
            urootw = u * np.sqrt(w)
            vrootw = v * np.sqrt(w)
            
            self.velHub[turb_name] = np.array([urootw, vrootw])
            self.velHubIni[turb_name] = np.array([urootw, vrootw])
            
            # Compute TIH, tubulence intensity at hub (%)
            # TR: assuming TKE uniform throughout the water column height and
            # TI = sqrt(2/3 * k) / U
            if w == 0.:
                wTI = 0.
            else:
                wTI = 1. / w # TR: due to above assumption
            
            self.features[turb_name]['TIH'] = wTI * ti
            
            if debug: module_logger.info("Computing yawing angles...")
            
            # current angle of attack
            bearing = vector_to_bearing(u, v)
            psi_current = bearing_to_radians(bearing) + np.pi
            
            # turbine direction
            psi_turb = bearing_to_radians(self.features[turb_name]['OA'])
            
            # maximum yaw
            psi_yaw = np.abs(np.radians(self.features[turb_name]['HAS'])) / 2.0
            
            ry = _get_angle_of_attack(psi_turb,
                                      psi_yaw,
                                      psi_current,
                                      two_way=self.features[
                                                          turb_name]['2way'])
            
            self.features[turb_name]['RY'] = np.degrees(ry)
            
            if debug:
                
                logMsg = ("Relative yawing angle for turbine{} = "
                          "{}").format(i, self.features[turb_name]['RY'])
                module_logger.info(logMsg)
        
        return


###Functions definitions###
def wp2_tidal(data,
              turbines,
              features,
              cfd_data,
              U_dict,
              V_dict,
              TKE_dict,
              debug=False,
              debug_plot=False):
    """
    DTOcean tidal model

    Args:
      data (dict) dictionary gathering site information with the following entries:
        . 'bathy' = bathymetry (m), 2D numpy array, dimension: (ny,nx)
        . 'SSH' = sea surface elevation (m), 2D numpy array, dimension: (ny,nx)
        . 'TI' = depth averaged turbulence intensity (m), 2D numpy array, dimension: (ny,nx)
        . 'U' = depth averaged West-East velocity component (m/s),
              2D numpy array, dimension: (ny,nx)
        . 'V' = depth averaged South-North velocity component (m/s),
              2D numpy array, dimension: (ny,nx)
        . 'X' = West-East corrdinates (m), 1D numpy array, dimension: (nx)
        . 'Y' = South-North corrdinates (m), 1D numpy array, dimension: (ny)
          'beta' = bed roughness
          'alpha' = power law exponent
      turbines (dict): dictionary gathering turbines' locations (m), 1D numpy array [x,y,z].
      features (dict): dictionary gathering individual turbine's features with the following entries:
        . 'Cp' = Power curve, Cp vs flow speed, list of two 1D numpy arrays [speed, Cp] where
               speed's and Cp's dimensions are identical
        . 'Cp' = Thrust curve, Ct vs flow speed, list of two 1D numpy arrays [speed, Ct] where
               speed's and Ct's dimensions are identical
        . 'cutIO' = cut-in and cut-out speeds (m/s), 1D numpy array [cut-in,cut-out]
        . 'Diam' = rotor diameter (m), float
        . 'floating' = floating-turbine-or-not info., boolean
        . 'HAS' = heading angle span (deg.), float

    Kwargs:
      debug (bool): debug flag
      debug_plot (bool): debug plot flag

    Returns:
      pow_perf_dev_no_int (dict): dictionary gathering turbines' power performance (MW) without interaction, float.
      pow_perf_dev )dict): dictionary gathering turbines' power performance (MW), float.
      pow_perf_array_no_int (float): Array power performance (MW) without interaction, float
      pow_perf_array (float): Array power performance (MW), float
      impact (float): ratio between dissipated and available mass flow rate from the flow
      for the considered statistical bin, [0 to 1]
      ti (dict) turbulence intensities at hub height for each turbine

    Notes: the 'turbines', 'features' and output dictionaries have the following structure:
      turbines[turbine's ID]['position'] where turbine's ID = 'turbine{number}'
      where the numbers should be padded so all ids are the same length


    """
    #performance benchmark
    if debug: start = time.time()
    
    # Check that the keys of turbines are the same length
    key_lengths_set = set([len(k) for k in turbines.keys()])
    
    if not len(key_lengths_set) == 1:
        
        err_msg = ("Key lengths in turbines argument differ. Ensure turbine "
                   "numbers are padded with zeros.")
        raise ValueError(err_msg)
    
    # Check that the keys of turbines and features are equal
    turbines_keys_set = set(turbines.keys())
    features_keys_set = set(features.keys())
    
    if turbines_keys_set != features_keys_set:
        
        err_msg = ("The arguments 'turbines' and 'features' have "
                   "non-matching keys")
        raise ValueError(err_msg)
    
    if debug: module_logger.info("initiating classes...")
    
    hydro = Hydro(data, debug=debug, debug_plot=debug_plot)
    array = Array(hydro,
                  turbines,
                  features,
                  debug=debug,
                  debug_plot=debug_plot)
    
    # In case there is only one turbine in the array
    NbTurb = len(array.positions.keys())
    if not NbTurb == 1:
        if debug: module_logger.info("Starting solver...")
        interaction = WakeInteraction(hydro,
                                      array,
                                      cfd_data,
                                      U_dict,
                                      V_dict,
                                      TKE_dict,
                                      debug=debug,
                                      debug_plot=debug_plot)
        interaction.solve_flow(debug=debug)
    else:
        if debug: module_logger.info("Single turbine...no interactions")
        
    if debug: module_logger.info("Computing performances...")
    arrayYield = ArrayYield(array, debug=debug, debug_plot=debug_plot)
    arrayYield.performance()

    if debug: module_logger.info("Computing hydrodynamic impacts...")
    # In case there is only one turbine in the array
    if not NbTurb == 1:
        impacts = HydroImpact(array,
                              hydro,
                              interaction,
                              arrayYield,
                              debug=debug,
                              debug_plot=debug_plot)
        impacts.dissipated_over_available_flux_ratio(debug=debug)
    else:
        impacts = HydroImpact(array,
                              hydro,
                              None,
                              arrayYield,
                              debug=debug,
                              debug_plot=debug_plot)
        impacts.dissipated_over_available_flux_ratio(debug=debug)

    #Outputs
    #  performances
    pow_perf_dev = arrayYield.turbine_capacity
    pow_perf_dev_no_int = arrayYield.turbine_capacity_no_interaction
    for key in pow_perf_dev.keys():
        pow_perf_dev[key] = pow_perf_dev[key] / 1e6
        pow_perf_dev_no_int[key] = pow_perf_dev_no_int[key] / 1e6
    pow_perf_array = arrayYield.array_capacity / 1e6
    pow_perf_array_no_int = arrayYield.array_capacity_no_interaction / 1e6
    #  impacts
    ratio = impacts.diss_avai_mass_flow_rate
    ti = {}
    
    n_digits = len(str(NbTurb))
    
    for i in range(NbTurb):
        turb = 'turbine{:0{width}d}'.format(i, width=n_digits)
        ti[turb] = array.features[turb]['TIH']

    #exit function
    if debug:
        end = time.time()
        logMsg = "Overall computation time: {} seconds.".format(end - start)
        module_logger.info(logMsg)

    return (pow_perf_dev_no_int,
            pow_perf_dev,
            pow_perf_array_no_int,
            pow_perf_array,
            ratio,
            ti)


def _get_angle_of_attack(psi_turb, psi_yaw, psi_current, two_way=False):
    
    if _is_within_yaw(psi_turb, psi_yaw, psi_current, two_way=two_way):
        return 0.0
    
    check_angles = []
                
    for a in [-1, 1]:
        
        x = (psi_current - psi_turb
                                 + np.pi) % (2 * np.pi) - np.pi + a * psi_yaw
        
        check_angles.append(abs(x))
    
    if two_way:
        
        for a in [-1, 1]:
        
            x = (psi_current - psi_turb) % (2 * np.pi) - np.pi + a * psi_yaw
            
            check_angles.append(abs(x))
    
    return min(check_angles)


def _is_within_yaw(psi_turb, psi_yaw, psi_current, two_way=False):
    
    if psi_yaw < 0. or psi_yaw > np.pi:
        
        psi_yaw = ("Given psi_yaw={} exceeds the valid range of "
                   "[0, pi]").format(psi_yaw)
        raise ValueError(psi_yaw)
    
    if psi_yaw == np.pi: return True
    
    check_uni = (psi_current - psi_turb + np.pi) % (2 * np.pi) - np.pi
    
    if abs(check_uni) <= psi_yaw: return True
    
    if not two_way: return False
    
    check_bi = (psi_current - psi_turb) % (2 * np.pi) - np.pi
    
    if abs(check_bi) <= psi_yaw: return True
    
    return False


def test():
    """
    Test function.
    Reproduce a two-device case in a uniform flow.
    """

    module_logger.info("Starting computation...")
    
    # Read and load CFD database into a dataframe at import level
    data_path = '/home/thomas/Desktop/Bitbucket/tidal_wake_data/'
    cfd_data = read_database(data_path)

    # Upload inputs
    ## Velocity field and site data
    xmax = 840.0
    ymax = 350.0
    lease = np.asarray([[0.0 , 0.0 ],
                    [0.0 , ymax],
                    [xmax, ymax],
                    [xmax, 0.0 ]])

    x = np.linspace(0.0, xmax, (xmax/10)+1) # dx = 10 m
    y = np.linspace(0.0, ymax, (ymax/10)+1) # dy = 10 m
    X, Y = np.meshgrid(x,y)
    BR = 1.0  # blockage ratio

    umax = 3.69 # maximum velocity in the X direction
    vmax = 0.0 # maximum velocity in the Y direction
    sshmax = 0.0 # maximum sea surface elevation
    timax= 0.1
    bathy = -31.5

    U = np.ones(X.shape) * umax
    V = np.ones(X.shape) * vmax
    SSH = np.ones(X.shape) * sshmax
    TI = np.ones(X.shape) * timax
    BATHY = np.ones(X.shape) * bathy

    data = {}
    data['TI'] = TI
    data['X'] = x  # save only 1D array as structured grid assumed
    data['Y'] = y  # save only 1D array as structured grid assumed
    data['U'] = U
    data['V'] = V
    data['SSH'] = SSH
    data['bathy'] = BATHY
    data['BR'] = BR
    data['lease'] = lease
    data['beta'] = 0.4
    data['alpha'] = 7.

    ## Turbines positions
    z = bathy/2.0 # hub height/depth
    coords = {}
    diam = 18.9
    first_row = 420.0 # x position of first row
    # turbine position
    coords['turbine0'] = {}
    coords['turbine0']['position'] = np.asarray((first_row, (ymax/2.0), z))
    coords['turbine1'] = {}
    coords['turbine1']['position'] = np.asarray((first_row + 5.0 * diam,
                                                (ymax/1.5), z))

    ## Turbines features
    cut_in = 0.0 # cut-in speed
    cut_out = 10.0 # cut-out speed
    # actual turbine features
    speed = np.arange(0.0, 10.0, 0.2)

    CT = np.ones(speed.shape) * 0.76
    Ct = [speed,CT] # thrust curve

    CP = np.ones(speed.shape) * 0.3
    Cp = [speed, CP] # Power curve
    feat = {}

    for key in coords.keys():
        feat[key] = {}
        feat[key]['OA'] = 270.0  # orientation angle (deg.), turbines face the West
        feat[key]['HAS'] = 0.0  # heading angle span (deg.), max = 180 deg. = full yaw
        feat[key]['Cp'] = Cp[:]  # Power curve
        feat[key]['Ct'] = Ct[:]  # thrust curve
        feat[key]['Diam'] = diam  # Diam = rotor diameter (m)
        feat[key]['cutIO'] = np.array([cut_in, cut_out])
        feat[key]['floating'] = False  # so hub height will be considered from the seabottom upwards
        feat[key]['2way'] = False  # turbines work in both direction.


    (pow_perf_dev_no_int,
     pow_perf_dev,
     pow_perf_array_no_int,
     pow_perf_array,
     ratio,
     ti) = wp2_tidal(data,
                     coords,
                     feat,
                     cfd_data,
                     debug=True,
                     debug_plot=False)
    module_logger.info("Array performance: " + str(pow_perf_array) + " MWatts")
    module_logger.info("Ratio of dissipated-over-available energy: "
                        + str(ratio*100) + " %")

    return
