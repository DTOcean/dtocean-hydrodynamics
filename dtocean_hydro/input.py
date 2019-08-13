# -*- coding: utf-8 -*-

#    Copyright (C) 2016 Francesco Ferri
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

"""
This module contains the classes used to collect and modify the WP2 inputs

.. module:: input
   :platform: Windows
   :synopsis: Input module to DTOcean WP2
   
.. moduleauthor:: Francesco Ferri <ff@civil.aau.dk>
.. moduleauthor:: Mathew Topper <mathew.topper@dataonlygreater.com>
"""

#Built in module inport
from __future__ import division

# Start logging
import logging
module_logger = logging.getLogger(__name__)

import os
import copy
from pprint import pformat
from shapely.geometry import Polygon, Point

# External module inport
import numpy as np
from scipy import interpolate


# Local module import
from dtocean_hydro.utils.Visualise_polygons import *
from dtocean_hydro.utils.bathymetry_utility import get_unfeasible_regions
from dtocean_hydro.utils.set_wdirs_multibody import anglewrap, convertangle

import dtocean_wave.utils.read_bem_solution as reader


class WP2_SiteData:
    """
    #------------------------------------------------------------------------------  
    #------------------------------------------------------------------------------
    #------------------ WP2 site data class
    #------------------------------------------------------------------------------
    #------------------------------------------------------------------------------
    WP2_SiteData class: The class contains all the information relative to the area of deployment 
     of the array.
    
     Args:
        LeaseArea (numpy.ndarray) [m]: UTM coordinates of the lease area poligon expressed as [X(Northing),Y(Easting)]. 
                                   X is the column vector containing the easting coordinates
                                   Y is the column vector containing the northing coordinates
        NogoAreas (list) [m]: list containing the UTM coordinates of the nogo areas poligons expressed as [X(Northing),Y(Easting)].
        MeteoceanConditions (dict): dictionary gathering all the information related to the metocean conditions 
                                    of the site. The dictionary is different for wave and tidal cases:
                                    Wave keys:
                                        'Te' (numpy.ndarray)[s]: Vector containing the wave energy periods 
                                        'Hs' (numpy.ndarray)[m]: Vector containing the significant wave height
                                        'dir' (numpy.ndarray)[rad]: Vector containing the wave direction
                                        'p' (numpy.ndarray)[-]: Probability of occurence of the sea state
                                        'specType' (tup): description of the spectral shape recorded at the site
                                                                    (spectrum name, gamma factor, spreading parameter)
                                        'SSH' (float)[m]: Sea Surface Height wrt the bathymetry datum at a single point
                                    Tidal keys:
                                       'V' (numpy.ndarray)[m/s]: northing-direction of the velocity field
                                        'U' (numpy.ndarray)[m/s]: easting-direction of the velocity field
                                        'p' (numpy.ndarray)[-]: probability of occurency of the state
                                        'TI' (numpy.ndarray)[-]: Turbulence intensity. It can be a single float or a matrix where the 
                                                  TI is specified at each grid node.
                                        'x' (numpy.ndarray)[m]: Vector containing the easting coordinate of the grid nodes
                                        'y' (numpy.ndarray)[m]: Vector containing the northing coordinate of the grid nodes
                                        'SSH' (numpy.ndarray)[m]: Sea Surface Height wrt the bathymetry datum
        
        Beta (float, optional if None): TIDAL ONLY bed roughness (default = 0.4)
        Alpha (float, optional if None): TIDAL ONLY power law exponent (default = 7.)
        Main_Direction (numpy.ndarray, optional) [m]: xy vector defining the main orientation of the array. If not provided it will be 
                                                        assessed from the Metocean conditions, expressed as [X(Northing),Y(Easting)].
        Bathymetry (numpy.ndarray) [m]: Describes the vertical profile of the sea bottom at each (given) UTM coordinate. 
                                        Expressed as [X(Northing),Y(Easting),Z(Down)]
        BR (float) [-]: describes the ratio between the lease area surface over the site area surface enclosed in a channel. 
                        1. - closed channel
                        0. - open sea
        electrical_connection_point (numpy.ndarray) [m]: UTM coordinates of the electrical connection point at the shore line, expressed as [X(Northing),Y(Easting)].
        boundary_padding (float, optional) [m]: Padding added to inside of the lease area in which devices may not be placed 

    Attributes:
        Same as the Args list
        mainAngle (numpy.ndarray)[rad]: defines the orientation angle of the array.
                                     The angle is calculated from the north in a clockwise direction
                                     0 - north
                                     pi/2 - east
                                     pi - south
                                     3/2 pi - west
    
    Returns:
         None
    """  
    
    def __init__(self,LeaseArea,
                      NogoAreas,
                      MeteoceanConditions,
                      Alpha,
                      Beta,
                      Main_Direction,
                      Bathymetry,
                      BR,
                      electrical_connection_point,
                      boundary_padding=None):
                          
        self.LeaseArea = LeaseArea
        self.NogoAreas = NogoAreas
        self.MeteoceanConditions = MeteoceanConditions
        self.Alpha = Alpha
        self.Beta = Beta
        self.Main_Direction = Main_Direction
        self.Bathymetry = Bathymetry
        self.BR = BR
        self.electrical_connection_point = electrical_connection_point
        self.boundary_padding = boundary_padding
            
    def printInput(self, indent=4):
        
        """print the Site class input arguments
        
        Args:
            indent(int,Optional): define the indent level of the printed string
        """
        vars_string = pformat(vars(self), indent=indent)
        logMsg = '--> SITE INPUT SUMMARY:\n\n{}'.format(vars_string)
        module_logger.info(logMsg)
        
        return
        

        
class WP2_MachineData:
    """
    #------------------------------------------------------------------------------  
    #------------------------------------------------------------------------------
    #------------------ WP2 machine data class
    #------------------------------------------------------------------------------
    #------------------------------------------------------------------------------  
    
    MachineData class: The class contains all the information relative to the machine 
    deployed in the array.
    
    Args:
        Type (str)[-]: defines the type of device either 'tidal' or 'wave'. No other strings are accepted
        lCS (numpy.ndarray)[m]: position vector of the local coordinate system from the given reference point.
                                Wave: represents the position vector of the body CS wrt the mesh CS
                                Tidal: represents the position of the hub from the machine (reference) CS.
        Clen (numpy.ndarray)[m]: characteristic lenght of the device:
                         Wave: unused
                         Tidal: turbine diameter and (optional) distance of the hub from the center line
                                used in case the machine is composed by two parallel turbines
        YawAngle (float)[rad]: Yaw angle span, wrt the main direction of the array. 
                            The total yawing range is two-fold the span. -Yaw/+Yaw
        Float_flag (bool)[-]: defines whether the machine is floating (True) or not (False)
        InstalDepth (list)[m]: defines the min and maximum water depth at which the device can be installed
        MinDist (tuple)[m]: defines the minimum allowed distance between devices in the array configuration
                            the first element is the distance in the x axis
                            the second element is the distance in the y axis
        OpThreshold (float)[-]: defines the minimum allowed q-factor
        UserArray (dict): dictionary containing the description of the array layout to be optimise. Keys:
                        'Option' (int): 1-optimisation over the internal parametric array layouts
                                        2-fixed array layout specified by the user not subject to optimisation
                                        3-array layout specified by the user subject to optimisation via expantion of the array
                        'Value' options:
                                (str) 'rectangular'
                                (str) 'staggered'
                                (str) 'full'
                                (numpy.ndarray) [m]: [X(Northing),Y(Easting)] coordiantes of the device
        RatedPowerArray (float)[W]: Rated power of the array.
        RatedPowerDevice (float)[W]: Rated power of the single isolated device.
        UserOutputTable (dict, optional): dictionary of dictionaries where all the array layouts inputed and analysed by the user are
                                  collected. Using this option the internal WP2 calculation is skipped, and the optimisaton
                                  is performed in the given data. The dictionaies keys are the arguments of the WP2 Output class.
        wave_data_folder (string, optional): path name of the hydrodynamic results generate by the wave external module
        tidal_power_curve (numpy.ndarray, optional)[-]: Non-dimensional power curve coefficients function of the stream velocity
        tidal_trust_curve (numpy.ndarray, optional)[-]: Non-dimensional trust curve coefficients function of the stream velocity
        tidal_velocity_curve (numpy.ndarray, optional)[m/s]: Vector containing the stream velocity
        tidal_cutinout (numpy.ndarray, optional)[m/s]: contain the cut_in and cut_out velocity of the turbine.
                                                        Outside the cut IN/OUT velocity range the machine will not produce
                                                        power. The generator is shut down, but the machine will still interact
                                                        with the others.
        tidal_bidirectional (bool, optional)[-]: bidirectional working principle of the turbine
        tidal_data_folder (string, optional)[-]: Path to tidal device CFD data files

        Attributes:
        Same as the Args list plus:
            MaxNumDevices (float) [-]: max number of devices allowed in the array estimated from the array rated power
                                        and the machine rated power
    """
    
    def __init__(self, Type,
                       lCS,
                       Clen,
                       YawAngle,
                       Float_flag,
                       InstalDepth,
                       MinDist,
                       OpThreshold,
                       UserArray,
                       RatedPowerArray,
                       RatedPowerDevice,
                       UserOutputTable=None,
                       wave_data_folder=None,
                       tidal_power_curve=None,
                       tidal_thrust_curve=None,
                       tidal_velocity_curve=None,
                       tidal_cutinout=None,
                       tidal_bidirectional=None,
                       tidal_data_folder=None):
        
        self.Type = Type
        self.lCS = lCS
        self.Clen = Clen
        self.YawAngle = YawAngle
        self.Floatflag = Float_flag
        self.InstalDepth = InstalDepth
        self.MinDist = MinDist
        self.OptThreshold = OpThreshold
        self.wave_data_folder = wave_data_folder
        self.tidal_data_folder = tidal_data_folder
        self.tidal_power_curve = tidal_power_curve
        self.tidal_thrust_curve = tidal_thrust_curve
        self.tidal_bidirectional = tidal_bidirectional
        self.tidal_cutinout = tidal_cutinout
        self.tidal_velocity_curve = tidal_velocity_curve
        self.UserArray = UserArray
        self.RatedPowerArray = RatedPowerArray
        self.RatedPowerDevice = RatedPowerDevice
        self.MaxNumDevices = int(RatedPowerArray / RatedPowerDevice)
        self.UserOutputTable = UserOutputTable
                    
    def printInput(self, indent=4):
        """print the Machine class input arguments

        Args:
            indent(int,Optional): define the indent level of the printed string
        """
        
        vars_string = pformat(vars(self), indent=indent)
        logMsg = '--> MACHINE INPUT SUMMARY:\n\n{}'.format(vars_string)
        module_logger.info(logMsg)
        
        return


class WP2input:
    """
    #------------------------------------------------------------------------------
    #------------------------------------------------------------------------------
    #------------------ WP2 input class
    #------------------------------------------------------------------------------
    #------------------------------------------------------------------------------

    WP2input class: The class collect the machine and site input and perform both data modification
    and checks.

    Args:
        MachineData (class): class containing the machine specification
        SiteData (class): class containing the site specification
        debug (bool): Debug flag
        debug_plot (bool): Plot debug flag


    Attributes:
        M_data (class): copy of the input class
        S_data (class): copy of the input class
        debug (bool): copy of the input debug flag
        debug_plot (bool): copy of the input debug flag
        internalOptim (bool): flag identifying whether the user provided a result table or not
        stopWP2run (bool): flag to trigger the interruption of the WP2 run

    Note:
        The M_data class attribute(s) are modified as follow:
            .InstalDepth:   set to [-inf,0] if None
                            switch the numbers if the first element is bigger than the second

        The M_data class is updated with the following attribute(s):
            .tidal_flag (bool): switch the call between the tidal or wave modules
    """
    
    def __init__(self, MachineData, SiteData, debug=False, debug_plot=False):
        # The input class is composed by the Machinedata and Sitedata classes
        self.stopWP2run = False
        self.M_data = copy.deepcopy(MachineData)
        self.debug = debug
        self.debug_plot = debug_plot
        self.S_data = copy.deepcopy(SiteData)
        self.internalOptim = True

        if self.M_data.UserOutputTable is not None:
            self.internalOptim = False

        if len(self.M_data.MinDist) == 1:
                self.M_data.MinDist = (self.M_data.MinDist[0], self.M_data.MinDist[0])

        if self.S_data.NogoAreas is None:
            self.S_data.NogoAreas = []

        if 't' in self.M_data.Type or 'T' in self.M_data.Type:
            self.M_data.tidalFlag = True
        else:
            self.M_data.tidalFlag = False

        # check the inputs before the WP2 run
        self.compress_lease_area()
        (status, errStr) = self.checkInput()
        if status < 0:
            raise ValueError(errStr)

        if not self.M_data.tidalFlag:
            
            # convert the scatter diagram angle convention to fit the
            # East-North convention
            self.change_angle_convention()  
            
            # convert the Te given in the scatter diagram into Tp for later use
            self.S_data.MeteoceanConditions['Tp'] = convert_te2tp(
                            self.S_data.MeteoceanConditions['Te'],
                            self.S_data.MeteoceanConditions['specType'][0],
                            self.S_data.MeteoceanConditions['specType'][1])
            
        self.getMainAngle()  # identify the main angle (reference) of the array

        # InstallationDepth
        if self.M_data.InstalDepth is None:
            self.M_data.InstalDepth = [-np.inf,0]
        else:
            if self.M_data.InstalDepth[0] > self.M_data.InstalDepth[1]:
                self.M_data.InstalDepth = [self.M_data.InstalDepth[1],
                                           self.M_data.InstalDepth[0]]
        
        self.getInstallationAreaConstraints()  # identifies possible areas where the devices cannot be installed
        self.rated_power_consistency_check()  # check the rated power consistency
        self.wec_power_matrix_check()  # check if the WEC power matrix represents the given site
        
        self.printSummary()
            
    def checkInput(self):
        """
        checkInput: used to assess the validity of the input given to the WP2 prior to perform any calculation

        Args:
            self (class)

        Returns:
            status (int): identify whether an error is found (status<0) or not (status=0)
            errStr (list): error strings appended during the error occurence.
        """
        errStr = []
        status = 0
        
        if self.M_data.MaxNumDevices<1:
            status -= 1
            errStr.append('The ratio between array and machine rated powers '
                           'is below 1.\n'
                           'Verify the rated power inputs.')
        
        if not self.S_data.NogoAreas is None:
            if not isinstance(self.S_data.NogoAreas,(list, tuple)):
                status -= 1
                errStr.append('The input format of the nogo area is '
                              'incorrect.\n'
                              'Accepted format are: None or list.')
            else:
                if not len(self.S_data.NogoAreas):
                    self.S_data.NogoAreas = None
                else:
                    if not isinstance(self.S_data.NogoAreas[-1],(np.ndarray)):
                        status -= 1
                        errStr.append('The input format of the nogo area is '
                                      'correct, but the list element format '
                                      'is incorrect.\n'
                                      'Accepted format is: numpy.ndarray.')
            
        if not self.S_data.Main_Direction is None:
            if not isinstance(self.S_data.Main_Direction, (np.ndarray)):
                status -= 1
                errStr.append('The input format of the main direction is '
                              'incorrect.\n'
                              'Accepted format are: None or numpy.ndarray')
            elif not len(self.S_data.Main_Direction)==2:
                status -= 1
                errStr.append('The input format of the main direction is '
                              'incorrect.\n'
                              'Too many input values. Accepted dimension: 2.')
                                    
        if self.M_data.tidalFlag:
            
            if not isinstance(self.S_data.MeteoceanConditions['SSH'],
                              np.ndarray):
                status -= 1
                errStr.append('The SSH value should be numpy.ndarray for '
                              'tidal devices.')
                              
        else:
            
            try:
                sshfloat = float(self.S_data.MeteoceanConditions['SSH'])
                self.S_data.MeteoceanConditions['SSH'] = sshfloat
            except (ValueError, TypeError):
                status -= 1
                errStr.append('The SSH value should be float for wave '
                              'devices.')
        
        if not self.M_data.tidalFlag and self.M_data.wave_data_folder is None:
                status -= 1
                errStr.append('The wave data folder must be given for a wave '
                              'simulation.')
                                    
        if not self.M_data.tidalFlag and self.M_data.wave_data_folder is not None:
            fname = os.path.join(self.M_data.wave_data_folder,'wec_solution.h5')
            
            if not os.path.isfile(fname):
                 status -= 1
                 errStr.append('The wec_solution.h5 file is missing!\n'
                               'The wec_solution.h5 needs to be contained into '
                               'the DataFolder.')
         
        if self.M_data.tidalFlag and self.M_data.tidal_data_folder is not None:
            module_logger.info('Need to check on the user database format!!')

        if not self.M_data.InstalDepth is None:
            if not isinstance(self.M_data.InstalDepth,(list,np.ndarray,tuple)):
                status -= 1
                errStr.append('The installation depth constraints format is '
                              'incorrect.\n'
                              'Accepted format are: None, list, tuple, '
                              'numpy.ndarray.')
            else:
                if len(self.M_data.InstalDepth) == 1:
                    status-=1
                    errStr.append('The installation depth requires two '
                                  'water depth values but only one is given.')
                elif len(self.M_data.InstalDepth) == 2:
                    if self.M_data.InstalDepth[0] == self.M_data.InstalDepth[1]:
                        status-=1
                        errStr.append('The two installation depth values need '
                                      'to be dinstinct.')
                else:
                    status-=1
                    errStr.append('Too many values in the installation depth '
                                  'input.')
        if self.M_data.tidalFlag:
            # the compress lease area has been moved above to allow the polygon buffering required by the mooring footprint
            # self.compress_lease_area()
            self.check_site_data_boundary()
            self.irrotational_flow()
        
        Opt = self.M_data.UserArray['Option']
        Value = self.M_data.UserArray['Value']
        
        if Opt==2:
            if Value.shape == (2L,):
                Value = Value.reshape((1,2))
                self.M_data.UserArray['Value'] = Value
            
            if not Value.shape[1]==2:
                if Value.shape[0]==2:
                    Value = Value.T
                    self.M_data.UserArray['Value'] = Value
                
                else:
                    raise ValueError("The fixed array layout need to have",
                                     "a shape of nbodyx2.",
                                     "The given array shape is: ",
                                     "{}x{}".format(Value.shape[0], Value.shape[1]))
        
        if self.M_data.tidalFlag:
            if self.M_data.lCS[-1] < self.M_data.Clen[0]/2:
                raise ValueError("The turbine radius is larger than the ",
                                 "distance between the hub height and",
                                 " the reference point (either the mwl or sea",
                                 " bottom). ")
                pass
            
                                  
        errStr = "\n".join(errStr)

        return status, errStr

    def printSummary(self):
        """
        printSummary: used to report a summary of the user inputs

        Args:
            self (class)
        """
        self.M_data.printInput()
        self.S_data.printInput()
        
        
    def check_site_data_boundary(self):
        lease = self.S_data.LeaseArea
        lx_max = lease[:,0].max()
        lx_min = lease[:,0].min()
        ly_max = lease[:,1].max()
        ly_min = lease[:,1].min()
        
        char_len = self.M_data.Clen
        offset = 0.
        if not np.size(char_len) == 1:
            if char_len[1] > 0:  # MCT case, two turbines included
                offset = char_len[-1]
                
        x_max = self.S_data.MeteoceanConditions['x'].max()-offset
        x_min = self.S_data.MeteoceanConditions['x'].min()+offset
        y_max = self.S_data.MeteoceanConditions['y'].max()-offset
        y_min = self.S_data.MeteoceanConditions['y'].min()+offset
        
        if lx_max > x_max or lx_min < x_min or ly_max > y_max or ly_min < y_min:
            if offset > 0:
                err_str = ("The velocity field does not cover entirely the lease area, ",
                           "Including the turbine interdistance",
                              "The execution is terminated")
            else:
                err_str = ("The velocity field does not cover entirely the lease area. ",
                              "The execution is terminated")
        
            raise IOError(err_str)
        
    def irrotational_flow(self):
        U = self.S_data.MeteoceanConditions['U'].copy()
        V = self.S_data.MeteoceanConditions['V'].copy()
        x = self.S_data.MeteoceanConditions['x'].copy()
        y = self.S_data.MeteoceanConditions['y'].copy()
        u_sh = U.shape
        if not len(u_sh) == 3:
            raise IOError("The shape of the velocity field is incorrect",
                          "The matrix need to be nx, ny, n_seastate")
        if not u_sh[0] == len(x):
            if not u_sh[1] == len(x):
                raise IOError("The shape of the velocity field is incorrect",
                          "The matrix need to be nx, ny, n_seastate")
            U = np.transpose(U, (1, 0, 2))
            V = np.transpose(V, (1, 0, 2))
        
      
        for el in range(u_sh[2]):
            curl = self.__curl_2d(U[:,:,el], V[:,:,el], x, y)
            if np.sqrt(curl.max()**2+curl.min()**2) > 0.1:
                raise IOError("The input velocity field is highly rotational ",
                              "This violate the model assumption ",
                              "The execution is terminated")

    def __curl_2d(self, uu, vv, x, y):
        
        u = uu.copy() - np.nanmin(uu)
        if np.nanmax(u): u /= np.nanmax(u)
        u -= np.nanmean(u)
    
        v = vv.copy() - np.nanmin(vv)
        if np.nanmax(v): v /= np.nanmax(v)
        v -= np.nanmean(v)
    
        dx = np.mean(x[1:]-x[:-1])
        dy = np.mean(y[1:]-y[:-1])
        
        u = np.nan_to_num(u)
        v = np.nan_to_num(v)
                
        x_var = interpolate.interp2d(x, y, u.T, 'cubic')
        y_var = interpolate.interp2d(x, y, v.T, 'cubic')

        dvdx = (y_var(x + dx, y) - y_var(x - dx, y)) / (2 * dx)
        dudy = (x_var(x, y + dy) - x_var(x, y - dy)) / (2 * dy)
            
        return dvdx - dudy
        
    def getInstallationAreaConstraints(self):
        """
        getInstallationAreaConstraints: used to identify the nogo zones associated with the
        water depth installation constraints of the specific machine.

        Args:
            self (class)

        Note:
            the method update the following self.S_data attributes:
                .NogoAreas_bathymetry: create a list of polygons associated with unfeasible and feasible areas
                .Bathymetry: flatten the bathymetry for the wave case to the average value
        """
        Bathymetry = self.S_data.Bathymetry
        # Bathymetry = np.array([-10])
        # print(Bathymetry)
        if (len(Bathymetry) == 1
            and (not Bathymetry >= self.M_data.InstalDepth[0]
            or not Bathymetry <= self.M_data.InstalDepth[1])):
            errStr = ('Error[InstalDepth]:\nThe device installation '
                    'constraints do not fit the bathymetry of the area.'
                    "\nNo possible installation area has been found!")
            raise ValueError(errStr)  # there is no valid zone for the installation of the devices
            
        elif (len(Bathymetry) == 1
            and Bathymetry >= self.M_data.InstalDepth[0]
            and Bathymetry <= self.M_data.InstalDepth[1]):
            self.S_data.NogoAreas_bathymetry = None  # the whole lease area can be used
            if self.M_data.tidalFlag:
                X, Y = np.meshgrid(self.S_data.MeteoceanConditions['x'],
                                   self.S_data.MeteoceanConditions['y'])
                Z = X*0.+Bathymetry  # generate a flat bathymetry
                self.S_data.Bathymetry = np.vstack((X.ravel(), Y.ravel(), Z.ravel())).T
            pass
        else:
            (NoGo,
                unfeasible_points_mask) = get_unfeasible_regions(
                                            Bathymetry,
                                            self.M_data.InstalDepth,
                                            debug=self.debug,
                                            debug_plot=self.debug_plot)
            self.S_data.NogoAreas_bathymetry = NoGo                 

        if not self.M_data.tidalFlag:
            module_logger.warning('[Warning] The wave module cannot run with '
                                  'variable bathymetry\n'
                                  'The bathymetry is reduced to its average '
                                  'value.')
            module_logger.info('The averge bathymetry value is '
                               '{} m'.format(np.mean(Bathymetry[:,-1])))
            
            # calculate the average water depth withint he lease area and
            # outside the nogozones specified by the user
            active_area = Polygon(self.S_data.LeaseArea)
            mask_outofwater = Bathymetry[:,-1] <= 0  # true all points below swl
            mask_nan = np.logical_not(np.isnan(Bathymetry[:,-1]))  # true all valid points
            mask_lease = np.asarray([active_area.intersects(Point(el))
                                            for el in Bathymetry[:,:-1]])  # true all points inside
            mask_nogo = np.ones(mask_lease.shape, dtype='bool')
            
            if self.S_data.NogoAreas:
                if not isinstance(self.S_data.NogoAreas, list):
                    raise IOError('The nogo areas input by the user needs to '
                                  'be a list of numpy.ndarray')                        
                for el in self.S_data.NogoAreas:
                    nogo_poly = Polygon(el)
                    t = np.asarray([not nogo_poly.intersects(Point(grid_p))
                                            for grid_p in Bathymetry[:,:-1]])
                    mask_nogo *= t  # true all point outside
                        
            bathymetry_mask = (mask_lease *
                               mask_nan *
                               mask_outofwater *
                               mask_nogo)  # only true point satisfy the above conditions
            
            if Bathymetry[bathymetry_mask,-1].size == 0:
                raise ValueError("No valid points found within lease area. "
                                 "Check depths and exclusion zones")
            
            self.S_data._average_water_depth = np.mean(Bathymetry[
                                                        bathymetry_mask,-1])
                    
            return

    def getMainAngle(self):
        """
        getMainAngle: calculates the main angle used to orient the array layouts

        Args:
            self (class)

        Note:
            the method update the following self.S_data attributes:
                .mainAngle (float): orientation angle (rad) in the North-Easting global convention

        """
        if self.S_data.Main_Direction is None:
            # if the main direction is not given it is taken from the scatter diagram
            if self.M_data.tidalFlag:
                # for the tidal case the main direction is given by
                # the negative of the mean stream direction (U,V) of the most
                # representative sea condition
                ind_max_direction = np.argmax(self.S_data.MeteoceanConditions["p"])
                U = np.nanmean(self.S_data.MeteoceanConditions["U"]
                                                    [:, :, ind_max_direction])
                V = np.nanmean(self.S_data.MeteoceanConditions["V"]
                                                    [:, :, ind_max_direction])
                self.S_data.Main_Direction = -np.array([U,V])
            else:
                # for the wave case the main direction is given by
                # the direcion of the scatter diagram with highest probability of occurrence
                ind_max_direction = np.argmax(np.sum(self.S_data.MeteoceanConditions["p"],(0,1)))
                angle = self.S_data.MeteoceanConditions["B"][ind_max_direction]
                self.S_data.Main_Direction = np.array([np.cos(angle), np.sin(angle)],dtype = float)

        self.S_data.mainAngle = anglewrap(np.arctan2(self.S_data.Main_Direction[1],
                                                     self.S_data.Main_Direction[0],
                                                     ),
                                          conversion='r2r')

        if self.debug:
            strDEBUG = 'The main angle is \n'
            strDEBUG += str(self.S_data.mainAngle*180/np.pi)
            strDEBUG += '\n'+'*'*100
            strDEBUG += 'The main direction is \n'
            strDEBUG += str(self.S_data.Main_Direction)
            strDEBUG += '\n'+'*'*100
            module_logger.debug(strDEBUG)

    def change_angle_convention(self):
        """
        change_angle_convention: convert the angles given in the metocean convention into
            the East-North (UTM) coordinate system

        Args:
            self (class)

        Note:
            the method update the following self.S_data attributes:
                .MeteoceanConditions (dic): wave directions converted from South-West to East-North

            the method update the following self.M_data attributes:
                .SingleMachine (dic): wave directions converted from South-West to East-North
        """
        old_beta_array = self.S_data.MeteoceanConditions['B']
        new_beta_array = convertangle(old_beta_array)
        self.S_data.MeteoceanConditions['B'] = new_beta_array

    def compress_lease_area(self):
        
        char_len = self.M_data.Clen
        
         # MCT case, two turbines included
        if not np.size(char_len) == 1 and char_len[1] > 0: 
            MCT_buffer = char_len[-1]
        else:
            MCT_buffer = 0
            
        overall_buffer = max(MCT_buffer, abs(max(0,self.S_data.boundary_padding)))
            
        lease_pol = Polygon(self.S_data.LeaseArea)
        
        # always contraction
        lease_pol_buffer = lease_pol.buffer(-overall_buffer)
        
        
        
#        pol_max = pol.max(0)
#        pol_min = pol.min(0)
#        scale = (pol_max-pol_min)
#        offset = np.array([0.5, 0.5])
#        pol_norm = (pol-pol_min)/scale-offset
#        comp_norm = 1-(2*compression*1.05)/scale
#        pol_new = (pol_norm*comp_norm+offset)*scale+pol_min
        
        self.S_data.LeaseArea = np.c_[lease_pol_buffer.exterior.xy]
        self.S_data._LeaseArea = np.c_[lease_pol.exterior.xy]
    
    def rated_power_consistency_check(self):
        rated_power = self.M_data.RatedPowerDevice
        
        if self.M_data.tidalFlag:
            cp_norm = self.M_data.tidal_power_curve
            v = self.M_data.tidal_velocity_curve
            cIO = self.M_data.tidal_cutinout
            diam = self.M_data.Clen[0]
            area = np.pi * ((diam/2.0)**2.0)
            
            cp_bound = cp_norm[np.logical_and(v>=min(cIO), v<=max(cIO))]
            v_bound = v[np.logical_and(v>=min(cIO), v<=max(cIO))]
            
            cp = 0.5*1025.0*cp_bound*area*(v_bound**3.0)
            machine_power = max(cp)
        
        else:
            (Cfit,
             Kfit,
             CPTO,
             Kmooring,
             w_te,
             w_hm0,
             w_dirs,
             scatter_diagram, 
             power_matrix ) = reader.read_performancefit_solution(
                                                self.M_data.wave_data_folder)
            
             
            machine_power = power_matrix.max()
             
        if rated_power < machine_power:
            strPowerWarning = ("The rated power specified ({} W) "
                               "is smaller than the rated power calculated "
                               "from the machine data ({} W). The power will "
                               "be truncated at the rated "
                               "value.").format(rated_power,
                                                machine_power)
            module_logger.warning(strPowerWarning)
        
    def wec_power_matrix_check(self):
        if not self.M_data.tidalFlag:
            (Cfit,
             Kfit,
             CPTO,
             Kmooring,
             w_te,
             w_hm0,
             w_dirs,
             scatter_diagram, 
             power_matrix ) = reader.read_performancefit_solution(
                                                self.M_data.wave_data_folder)
            
            w_tp = convert_te2tp(
                            w_te,
                            self.S_data.MeteoceanConditions['specType'][0],
                            self.S_data.MeteoceanConditions['specType'][1])
             
            site = self.S_data.MeteoceanConditions
            s_tp = site['Tp']
            s_hm0 = site['Hs']
            s_dirs = site['B']
             
            if max(w_tp) < max(s_tp):
                strTpWarning = ("The range of Wave Periods specified in the machine power matrix",
                                "does not cover the given site",
                                "Due to the model linearity, this situation can bring unexpected/unrealistic results")
                module_logger.warning(strTpWarning)
                
            if max(w_hm0) < max(s_hm0):
                strHm0Warning = ("The range of Wave Heighs specified in the machine power matrix",
                                "does not cover the given site",
                                "Due to the model linearity, this situation can bring unexpected/unrealistic results")
                module_logger.warning(strHm0Warning)
            
            if max(w_dirs) < max(s_dirs):
                strDirsWarning = ("The range of Wave Directions specified in the machine power matrix",
                                "does not cover the given site",
                                "Due to the model linearity, this situation can bring unexpected/unrealistic results")
                module_logger.warning(strDirsWarning)
    
             
def convert_te2tp(te, specType, gamma):
    
    coeff = np.array([[  1.22139232e+00],
                      [ -7.26257028e-02],
                      [  1.74397331e-02],
                      [ -2.19288663e-03],
                      [  1.07357912e-04]])
                    
    # convert Te to Tp for the Metocean condition relative to the deployment site
    conversion_factor = 1.16450471
    
    if specType == 'Jonswap':
        if gamma > 7 or gamma < 1:
            module_logger.warning('The gamma value of the JONSWAP spectrum '
                                  'in the metocean data specification is out '
                                  'of the confident range [1-7].')
        
        conversion_factor = coeff[0] + \
                            coeff[1] * gamma + \
                            coeff[2] * gamma**2 + \
                            coeff[3] * gamma**3 + \
                            coeff[4] * gamma**4
                            
    tp = te * conversion_factor
    
    return tp
