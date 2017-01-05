#!/usr/bin/python2.7
# encoding: utf-8
from __future__ import division

# Start logging
import logging
module_logger = logging.getLogger(__name__)

import time
from shapely.geometry import Polygon
import numpy as np
import matplotlib.pyplot as plt

# Local import
from .utils import distance_from_streamline
from .utils.interpolation import interpol_scatter2grid, interp_at_point
from .utils.misc import pi2pi, deg360_to_radpi
#from .modules.vertical_velocity_profile import vvpw
# TR: alternative using Soulsby formulation
from .modules.vertical_velocity_profile import vvpw_soulsby
from .modules import Streamlines, ArrayYield, HydroImpact
from .submodel.WakeInteraction import WakeInteraction
from .submodel.ParametricWake import read_database

###Classes definition###
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

    """
    def __init__(self, data, debug=False, debug_plot=False):
        if debug: module_logger.info("Loading Hydrodynamic conditions...")
        self._debug = debug
        self.data = data
        self.U = self.data['U']
        self.V = self.data['V']
        self.TI = self.data['TI']
        self.PLE = self.data['PLE']
        self.SSH = self.data['SSH']
        self.bathy = self.data['bathy']
        self.geophy = self.data['geophy']
        self.X = self.data['X']
        self.Y = self.data['Y']
        self.lease = Polygon(self.data['lease'])
        (xm, ym, xM, yM) = self.lease.bounds
        self.bounding_box = Polygon([[xm, ym], [xM, ym], [xM, yM], [xm, yM]])
        self.BR = self.data['BR'] # lease surface / site area, user input

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
        """ Initialise Array class"""
        self._turbine_count = len(turbines.keys())
        #Turn turbine positions into one big matrix
        self.positions = {}
        for i in range(self._turbine_count):
            p = turbines['turbine' + str(i)]['position']
            if i == 0:
                positions = p
            else:
                positions = np.vstack((positions, p))
            self.positions['turbine' + str(i)] = p            
 
        #Turn turbine features into one big matrix
        self.features = {}
        for i in range(self._turbine_count):
            self.features['turbine' + str(i)] = features['turbine' + str(i)]

        # Streamlines
        #  generate regular structured grid U, V, X, Y to pass on to Streamlines
        data = {}
        x, y = np.meshgrid(hydro.X, hydro.Y)
        xn = hydro.X.shape[0]
        yn = hydro.Y.shape[0]
        data['X'], data['Y'], data['U'] = interpol_scatter2grid(
                                                            x.flatten(),
                                                            y.flatten(),
                                                            hydro.U.flatten(),
                                                            xn,
                                                            yn, 
                                                            debug=debug)
        data['X'], data['Y'], data['V'] = interpol_scatter2grid(
                                                            x.flatten(),
                                                            y.flatten(),
                                                            hydro.V.flatten(),
                                                            xn,
                                                            yn,
                                                            debug=debug)
        # In case there is only one turbine in the array
        if not self._turbine_count == 1:
            #  computes streamlines
            SLs = Streamlines(data,
                              positions,
                              self._turbine_count,
                              debug=debug)
            if debug_plot: SLs.plot(self._turbine_count)

            #Relative distance from streamlines
            self.streamlines={}
            self.distances={}
            if debug_plot:
                fig = plt.figure(figsize=(18,10))
                ax = fig.add_subplot(1,1,1)

            for i in range(self._turbine_count):
                if debug:
                    module_logger.info("Streamline nb: {}".format(i))
                streamline = SLs.streamlines[i]
                self.streamlines['turbine' + str(i)] = streamline
                self.distances['turbine' + str(i)] = \
                            distance_from_streamline(streamline,
                                                     positions,
                                                     self._turbine_count,
                                                     i,
                                                     debug=debug,
                                                     debug_plot=debug_plot)
            if debug_plot:
                SLs.plot(self._turbine_count)
                ax.set_ylabel('Distance (m)', fontsize = 12)
                ax.set_xlabel('Distance (m)', fontsize = 12)
                plt.show()
        else:
            if debug: module_logger.info("Single turbine...no streamlines")

        # Velocity and TI at hub
        self.velHub = {}
        self.velHubIni = {}  # save initial velocity for computing Q factor
        if debug: module_logger.info("Computing Hub velocities...")
        for i in range(self._turbine_count):
            if debug: module_logger.info("Interpolating quantities...")
            Q = [hydro.U,
                 hydro.V,
                 hydro.SSH,
                 hydro.bathy,
                 hydro.geophy,
                 hydro.PLE,
                 hydro.TI]
            [u, v, el, h, n, ple, ti] = interp_at_point(
                                        self.positions['turbine' + str(i)][0],
                                        self.positions['turbine' + str(i)][1],
                                        hydro.X, hydro.Y, Q)
            # quantity formatting
            if h < 0.0: h *= -1.0
            #  Hub height account for floating options
            hh = self.positions['turbine' + str(i)][2]  # hub height, negative by convention
            if hh > 0.0: hh *= -1.0
            if self.features['turbine' + str(i)]['floating']:
                z = (el + h) + hh
            else:
                z = -hh
            # Computing velocity vertical profile weight
            radius = self.features['turbine' + str(i)]['Diam'] / 2.0
            #wTop = vvpw(u, v, z+radius, el, h, n, debug=debug)
            #wBottom = vvpw(u, v, z-radius, el, h, n, debug=debug)
            # TR: alternative using Soulsby formulation
            wTop = vvpw_soulsby(z+radius, el, h, n, ple, debug=debug)
            wBottom = vvpw_soulsby(z-radius, el, h, n, ple, debug=debug)
            w = (wTop + wBottom) / 2.0  # Linear interpolation along the vertical
            self.velHub['turbine' + str(i)] = np.array([u * w, v * w])
            self.velHubIni['turbine' + str(i)] = np.array([u * w, v * w])
            # Compute TIH, tubulence intensity at hub (%)
            # TR: assuming TKE uniform throughout the water column height and TI = sqrt(2/3 * k) / U
            wTI = 1/w # TR: due to above assumption
            self.features['turbine' + str(i)]['TIH'] = wTI * ti

        # Computes actual yawing angles
        # Convention: angle between -pi and pi, 0 of trigonometric circle coincides with East
        if debug: module_logger.info("Computing yawing angles...")
        for i in range(self._turbine_count):
            # absolute angle
            absAngle = pi2pi(np.arctan2(v,u) - np.pi)
            # -pi in order to test turbine aligned with flow rather than facing flow
            turbDir = deg360_to_radpi(self.features['turbine' + str(i)]['OA'])
            turbSpan = np.abs(
                    np.radians(self.features['turbine' + str(i)]['HAS']))/2.0
            if not turbSpan == 180.0:  # case where turbine rotates 360.0 degrees
                # angle between -pi and pi
                turbDir = pi2pi(turbDir)
                # define intervals
                inter1 = self._interval(turbDir, turbSpan)
                # check if 2way turbine and define additional interval if needed
                if self.features['turbine' + str(i)]['2way']:
                    turbDir2 = pi2pi(turbDir + np.pi)
                    inter2 = self._interval(turbDir2, turbSpan)
                # check if flow direction with interval
                inFlag = False
                if type(inter1)==list:
                    if inter1[0][0]<=absAngle<=inter1[0][1]: inFlag=True
                    if inter1[1][0]<=absAngle<=inter1[1][1]: inFlag=True
                else:
                    if inter1[0]<=absAngle<=inter1[1]: inFlag=True
                if self.features['turbine' + str(i)]['2way']:
                    if type(inter2)==list:
                        if inter2[0][0]<=absAngle<=inter2[0][1]: inFlag=True
                        if inter2[1][0]<=absAngle<=inter2[1][1]: inFlag=True
                    else:
                        if inter2[0]<=absAngle<=inter2[1]: inFlag=True
            else:
                inFlag=True
            # compute relative yaw, RY
            if inFlag:
                self.features['turbine' + str(i)]['RY'] = 0.0
            else:
                bounds = np.array(inter1).flatten()
                if self.features['turbine' + str(i)]['2way']:
                    bounds = np.hstack((bounds, np.array(inter2).flatten()))
                #  find closest bound
                ry = bounds - absAngle
                self.features['turbine' + str(i)]['RY'] = np.degrees(ry.min())  # relative yaw angle in deg., RY
            if debug:
                logMsg = ("Relative yawing angle for turbine{} = "
                          "{}").format(i, 
                                       self.features['turbine' + str(i)]['RY'])                                       
                module_logger.info(logMsg)

    def _interval(self, turbDir, turbSpan):
        """
        define angle interval

        Args:
          turbDir (float): turbine's orientation
          turbSpan (float): turbine's yawing span

        Returns:
          inter1 (list): yawing angle's interval

        """
        inter1 = np.array([turbDir - turbSpan, turbDir + turbSpan])
        # angle between -pi and pi
        if np.any(inter1 >= np.pi):
            inter1 = [np.array([turbDir - turbSpan, np.pi]),
                     np.array([-np.pi, (turbDir + turbSpan) - 2.0*np.pi])]
        elif np.any(inter1 <= -np.pi):
            inter1 = [np.array([-np.pi, turbDir + turbSpan]),
                     np.array([(turbDir - turbSpan) + 2.0*np.pi, np.pi])]
        return inter1

###Functions definitions###
def wp2_tidal(data,
              turbines,
              features,
              cfd_data,
              debug=False,
              debug_plot=False):
    """
    DTOcean tidal model

    Args:
      data (dict) dictionary gathering site information with the following entries:
        . 'bathy' = bathymetry (m), 2D numpy array, dimension: (ny,nx)
        . 'geophy' = Manning's roughness coeff., 2D numpy array, dimension: (ny,nx)
        . 'PLE' = Power law exponent for velocity vertical profile,
                2D numpy array, dimension: (ny,nx)
        . 'SSH' = sea surface elevation (m), 2D numpy array, dimension: (ny,nx)
        . 'TI' = depth averaged turbulence intensity (m), 2D numpy array, dimension: (ny,nx)
        . 'U' = depth averaged West-East velocity component (m/s),
              2D numpy array, dimension: (ny,nx)
        . 'V' = depth averaged South-North velocity component (m/s),
              2D numpy array, dimension: (ny,nx)
        . 'X' = West-East corrdinates (m), 1D numpy array, dimension: (nx)
        . 'Y' = South-North corrdinates (m), 1D numpy array, dimension: (ny)
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
      turbines[turbine's ID]['position'] where turbine's ID = 'turbine'+str(integer)
      for integers going from 0 to infinity


    """
    #performance benchmark
    if debug: start = time.time()

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
                                      debug=debug,
                                      debug_plot=debug_plot)
        interaction.solv_induction(debug=debug)
    else:
        if debug: module_logger.info("Single turbine...no interactions")

    if debug: module_logger.info("Computing performances...")
    arrayYield = ArrayYield(array, debug=debug, debug_plot=debug_plot)
    arrayYield.performance()

    if debug: module_logger.info("Computing hydrodynamic impacts...")
    # In case there is only one turbine in the array
    if not NbTurb == 1:
        impacts = HydroImpact(array, hydro,
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
    for i in range(NbTurb):
        turb = 'turbine' + str(i)
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
    nx = len(x)
    ny = len(y)
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
    PLE = 12.0 * np.ones(X.shape)
    # bed roughness coefficient here
    manning = np.ones(X.shape) * 0.3

    data = {}
    data['TI'] = TI
    data['PLE'] = PLE
    data['geophy'] = manning
    data['X'] = x  # save only 1D array as structured grid assumed
    data['Y'] = y  # save only 1D array as structured grid assumed
    data['U'] = U
    data['V'] = V
    data['SSH'] = SSH
    data['bathy'] = BATHY
    data['BR'] = BR
    data['lease'] = lease

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
