# -*- coding: utf-8 -*-

#    Copyright (C) 2016 Francesco Ferri
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

"""
This module contains the class used to interface the tidal module developed
by Thomas Roc (ITPower).

.. module:: array
   :platform: Windows
   :synopsis: Tidal interface module for DTOcean WP2

.. moduleauthor:: Francesco Ferri <ff@civil.aau.dk>
.. moduleauthor:: Mathew Topper <mathew.topper@dataonlygreater.com>
"""

import logging

import numpy as np
from dtocean_hydro.output import ReducedOutput

from .main import wp2_tidal
from .utils.misc import bearing_to_radians, radians_to_bearing

# Start logging
module_logger = logging.getLogger(__name__)


class CallTidal:
    """
    CallTidal: Interface class generated for the integration of the tidal
    module into the overall WP2 package.

    Args:
         WP2Hydro (class): aggregation of the hydro class instantiated in the WP2 object
         WP2Input (class): aggregation of the input class instantiated in the WP2 object
         cfd_data (dict): cfd database as dictionary of pandas tables

    Optional args:
         debug (boolean): switch on and off the command line printing option
         debug_plot (boolean): switch on and off the plotting option

    Attributes:
      Pubblic:
        coord (numpy.ndarray): UTM coordinates (Easting, Northing) of the machines in the current layout
        hub_height (numpy.ndarray)[m]: represents the position of the hub from the machine (reference) CS.
        n_seastate (int)[-]: number of sea state representing the flow fields of the specific site
        power_prod_perD_perS (numpy.ndarray)[W]: average power production of each device per each sea
                                                    states WITH interaction
        power_prod_perD_perS_ni (numpy.ndarray)[W]: average power production of each device  per each sea
                                                        states WITHOUT interaction
        power_prod_perD (numpy.ndarray)[W]: average power production of each device over the sea states
        power_prod_array (numpy.ndarray)[W]: average power production of the array WITH interaction over
                                                the sea states
        power_prod_array_no_int (numpy.ndarray)[W]: average power production of the array WITHOUT interaction
                                                        over the sea states
        Nbodies (int)[-]: number of machine in the current array layout
        TI (list)[-]: list of dictionaries gathering the information of the turbulence
                        intensity of each device for each sea states.
                        The keys of each dictionary are:
                            "turbineID": where the ID is an integer with base count 0 up to Nb-1
        Resource_reduction (tuple)[Wh/Wh]: rate between in and out energies for the worst sea states
        AEP_perD (numpy.ndarray)[Wh]: average annual energy production of each device over the different sea states
        AEP_array (float)[Wh]: average annual energy production of the array over the different sea states
        q_array (float)[-]: average q-factor of the array over the different sea states
        q_perD (numpy.ndarray)[-]: average q-factor of each device over the different sea states
        Device_Model (unused)[-]: empty dictionary
        Device_Position (dict)[m]: dictionary gathering the information of the position of each machine in the array.
                                    For each key a tuple of UTM coordinates (Easting, Northing) is assigned.
                                    The dictionary keys are:
                                        "machineID": where the ID is an integer with base count 0 up to Nb-1
      Private:
        __turbines (dict)[m]: temporary data container - same format of the Device_Position
        __U (numpy.ndarray): x-contribution (Easting) to the velocity field at each grid node
        __V (numpy.ndarray): y-contribution (Northing) to the velocity field at each grid node
        __mct_flag (bool)[-]: True is the __mct_dist is greater than zero (MCT case)
        __mct_dist (float)[m]: distance of the hub from the center line
                                used in case the machine is composed by two parallel turbines
        __prob (numpy.ndarray)[-]: normalised probability of occurence of each sea states
        __u_mean (float)[m/s]: mean Easting velocity
        __v_mean (float)[m/s]: mean Northing velocity
        __TI (numpy.ndarray)[-]: turbulence intensity at each grid node for each sea state
        __SSH (numpy.ndarray)[m]: Sea Surface Height wrt the bathymetry datum, at each grid node for each sea state
        __base_feature (dict)[-]: contains the base feature of each turbine that need to be repeated for each body.
                                The dictionary is only used to reduce the line of code.
                                Keys:
                                    Ct (list)[-]:
                                        index 0-> (numpy.ndarray)[m/s]: Vector containing the stream velocity
                                        index 1-> (numpy.ndarray)[-]: Non-dimentional trust curve function of the stream velocity
                                    Cp (list)[W]:
                                        index 0-> (numpy.ndarray)[m/s]: Vector containing the stream velocity
                                        index 1-> (numpy.ndarray)[-]: Non-dimentional power curve function of the stream velocity
                                    Diam (numpy.ndarray)[m]: turbine diameter
                                    cutIO (numpy.ndarray)[m/s]: contain the cut_in and cut_out velocity of the turbine.
                                                        Outside the cut IN/OUT velocity range the machine will not produce
                                                        power. The generator is shut down, but the machine will still interact
                                                        with the others.
                                    floating (bool)[-]: defines whether the machine is floating (True) or not (False)
                                    HAS (float)[deg]: Yaw angle span of the machine.
                                    OA (float)[deg]: defines the orientation angle of the array.
                                                The angle is calculated from the north in a clockwise direction
                                                     0 - north
                                                     90 - east
                                                     180 - south
                                                     270 - west
                                    2way (bool): bidirectional working principle of the turbine
        __data (dict)[-]: contains the site data to be passed to the tidal solver, for the specific sea state
                        Keys:
                            U: see __U
                            V: see __V
                            TI: see __TI
                            SSH: see __SSH
                            bathy (numpy.ndarray) [m]: Describes the vertical profile of the sea bottom
                                                    at each (given) UTM coordinate.
                                                    Expressed as [X(Northing),Y(Easting),Z(Down)]
                            geophy (numpy.ndarray) [-]: Describes the sea bottom geophysic characteristic
                                                    at each (given) UTM coordinate.
                                                    Expressed as [X(Northing),Y(Easting),Geo]
                            X (numpy.ndarray)[m]: Vector containing the easting coordinate of the grid nodes
                            Y (numpy.ndarray)[m]: Vector containing the northing coordinate of the grid nodes
                            lease (numpy.ndarray) [m]: UTM coordinates of the lease area poligon
                                                        expressed as [X(Northing),Y(Easting)].
                            BR (float) [-]: describes the ratio between the lease area surface over the site
                                            area surface enclosed in a channel.
                                                    1. - closed channel
                                                    0. - open sea
                            beta (float): bed roughness
                            alpha (float): power law exponent
    
    Returns:
        CallTidal (class): return the current class object.
    """    
    def __init__(self,
                 WP2Hydro,
                 WP2Input,
                 cfd_data,
                 debug=False,
                 debug_plot=False):
        
        self.cfd_data = cfd_data
        self.Device_Model = None
        self.AEP_perD = None
        self.AEP_array = None
        self.q_array = None
        self.q_perD = None
        self.Device_Position = None
        self.__mct_flag = False
        self.__mct_dist = 0.
        self.hub_height = -WP2Input.M_data.lCS[-1]
        self.debug = debug
        self.debug_plot = debug_plot
        self.coord = None
        
        # machine features
        Cp = [WP2Input.M_data.tidal_velocity_curve,
              WP2Input.M_data.tidal_power_curve]
        Ct = [WP2Input.M_data.tidal_velocity_curve,
              WP2Input.M_data.tidal_thrust_curve]
        
        char_len = WP2Input.M_data.Clen
        Diam = char_len[0]
        
        # mapping the orientation angle in agreement with the tidal convention
        OA = radians_to_bearing(WP2Input.S_data.mainAngle)
        cutIO = WP2Input.M_data.tidal_cutinout
        floaT = WP2Input.M_data.Floatflag
        
        # the factor 2 is used to take the half span to the full span
        HAS = 2 * WP2Input.M_data.YawAngle * 180 / np.pi
        Bidirection = WP2Input.M_data.tidal_bidirectional
        
        # site features
        x = WP2Hydro.x
        y = WP2Hydro.y
        nx = len(x)
        ny = len(y)
        self.__U = WP2Hydro.U
        self.__V = WP2Hydro.V
        TI = WP2Hydro.TI
        SSH = WP2Hydro.wdepth
        self.__prob = WP2Hydro.p
        
        self.__u_mean = np.nanmean(self.__U)
        self.__v_mean = np.nanmean(self.__V)
        self.n_seastate = len(self.__prob)
        
        Bathy = points_to_grid(WP2Hydro.bathy, x, y).T
        rated_power = WP2Input.M_data.RatedPowerDevice
        
        # By default the MCT is not considered.
        # patch for the MCT case
        if not np.size(char_len) == 1:
            if char_len[1] > 0:  # MCT case, two turbines included
                # this flag is used to trigger the case of the MCT on/off.
                self.__mct_flag = True
                self.__mct_dist = char_len[1]
                rated_power /=  2.
        
        # conditional reconstruction of some input
        if len(TI) == 1:
            TI = TI * np.ones((ny, nx, self.n_seastate))
            
        if SSH.size == self.n_seastate:
            SSH = np.multiply(np.ones((ny, nx, self.n_seastate)), SSH)
        
        self.__TI = TI
        self.__SSH = SSH
        
        self.__base_feature = {'Ct': Ct,
                               'Cp': Cp,
                               'Diam': Diam,
                               'cutIO': cutIO,
                               'floating': floaT,
                               'HAS': HAS,
                               'OA': OA,
                               '2way': Bidirection,
                               'Rating': rated_power}
        
        self.__data = {'U': None,
                       'V': None,
                       'TI': None,
                       'SSH': None,
                       'bathy': Bathy,
                       'X': x,
                       'Y': y,
                       'lease': WP2Input.S_data._LeaseArea,
                       'BR': WP2Input.S_data.BR,
                       'beta': WP2Hydro.beta,
                       'alpha': WP2Hydro.alpha}
    
    def energy(self, coord):
        """
        method to assess the energy production of the array specified in the
        input argument.
        
        Args:
             coord (numpy.ndarray):
                 UTM coordinates (Easting, Northing) of the machines in the
                 current layout
             
        Returns:
            ReducedOutput Object: contains the minimum set of information
            required by the WP2Output class.
        """
        
        n = self.n_seastate
        self.__set_coordinates(coord)
        (nb, fea, turb) = self.__set_turbine_and_features()
        
        # initialise some support variables
        pow_perf_dev_state = np.zeros((nb, n))
        pow_perf_dev_state_ni = np.zeros((nb, n))
        pow_perf_array = np.zeros((n))
        pow_perf_array_no_int = np.zeros((n))
        resource_reduction_state = np.zeros((n))
        ti_dev_state = []
        
        U_dict = self.cfd_data['dfU'].to_dict()
        V_dict = self.cfd_data['dfV'].to_dict()
        TKE_dict = self.cfd_data['dfTKE'].to_dict()
        
        # iterate over the different sea states
        for seastate_id in range(self.n_seastate):
            
            module_logger.debug("Evaluating sea-state {}".format(seastate_id))
            
            self.__set_data(seastate_id)
            
            # evaluate the array performance for the current array layout and
            # sea state
            (pow_perf_dev_no_int,
             pow_perf_dev,
             pow_perf_array_no_int[seastate_id],
             pow_perf_array[seastate_id],
             resource_reduction_state[seastate_id],
             ti) = wp2_tidal(self.__data,
                             turb,
                             fea,
                             self.cfd_data,
                             U_dict,
                             V_dict,
                             TKE_dict,
                             debug=self.debug,
                             debug_plot=self.debug_plot)
            
            ti_dev_state.append(ti)
            
            pow_perf_dev_sorted = sorted(pow_perf_dev.items())
            pow_perf_dev_no_int_sorted = sorted(pow_perf_dev_no_int.items())
            
            for ii in range(nb):
                
                new_pow_perf_dev = pow_perf_dev_sorted[ii][1]
                new_pow_perf_dev_ni = pow_perf_dev_no_int_sorted[ii][1]
                
                pow_perf_dev_state[ii, seastate_id] = new_pow_perf_dev
                pow_perf_dev_state_ni[ii, seastate_id] = new_pow_perf_dev_ni
        
        # for the case of the MCT, the number of turbine has been doubled
        # Now the values are mapped back to the original number of devices
        if self.__mct_flag:
            (pow_perf_dev_state,
             pow_perf_dev_state_ni,
             ti_dev_state,
             turb, nb) = _map_mct_results(pow_perf_dev_state,
                                          pow_perf_dev_state_ni,
                                          ti_dev_state,
                                          turb)
        
        self.power_prod_perD_perS = pow_perf_dev_state * 1.0e6
        self.power_prod_perD_perS_ni = pow_perf_dev_state_ni * 1.0e6
        self.power_prod_array = pow_perf_array * 1.0e6
        self.power_prod_array_no_int = pow_perf_array_no_int * 1.0e6
        self.Nbodies = nb
        self.__turbines = turb
        self.TI = ti_dev_state
        self.Resource_reduction = np.nanmax(resource_reduction_state)
        
        return self.__output_tidal()
    
    def __output_tidal(self):
        """
        output_tidal: the method is used to map the tidal output to the WP2
        object format.
        
        Returns:
            the methods return the object itself with mapped attributes
        """
        
        machine = {}
        n_digits = len(str(self.Nbodies))
        
        for jj in range(self.Nbodies):
            
            stro = 'turbine{:0{width}d}'.format(jj, width=n_digits)
            pos = self.__turbines[stro]['position']
            strn = 'Device{}'.format(jj)
            machine[strn] = (pos[0], pos[1])
            
        year_hours = 365 * 24
        
        pow_dev = np.sum(self.power_prod_perD_perS * self.__prob, 1)
        aep_dev = year_hours * pow_dev
        aep_ar = np.sum(aep_dev)
        
        q_ar = np.sum(self.power_prod_array) / \
                                          np.sum(self.power_prod_array_no_int)
        q_dev = np.sum(self.power_prod_perD_perS, 1) / \
                                        np.sum(self.power_prod_perD_perS_ni, 1)
        
        # ReducedOutput arguments:
        #
        # aep_ar (float)[Wh]: annual energy production of the whole array
        # aep_dev (numpy.ndarray)[Wh]: annual energy production of each device
        #   within the array
        # q_dev (numpy.ndarray)[]: q-factor for each device, calculated as
        #   energy produced by the device within the array over the energy
        #   produced by the device without interaction
        # q_ar (float)[]: q-factor for the array, calculated as energy
        #   produced by the array over the energy produced by the device 
        #   without interaction times the number of devices.
        # pow_dev (numpy.ndarray)[W]: average power production of each
        #   device within the array
        # pow_dev_state (numpy.ndarray)[W]: average power production
        #   of each device within the array for each sea state
        # nb (float)[]: number of devices in the array
        # pos (numpy.ndarray)[m]: UTM coordinates of each device in
        #   the array. NOTE: the UTM coordinates do not consider different UTM
        #   zones. The maping in the real UTM coordinates is done at a higher
        #   level.
        # res_red (float)[]: ratio between absorbed and incoming
        #   energy.
        # ti (float)[TIDAL ONLY]: turbulence intensity within the array
        
        res = ReducedOutput(aep_ar,
                            aep_dev,
                            q_ar,
                            q_dev,
                            pow_dev,
                            self.power_prod_perD_perS,
                            self.Nbodies,
                            machine,
                            self.Resource_reduction,
                            self.TI)
        
        return res
    
    def __set_turbine_and_features(self):
        """
        Generates the required features and machine specification for each sea 
        state.
        
        Returns:
            (tuple):
                number of bodies (int)[-]:
                    number of bodies specified in the coord attribute
                machine features (dict)[-]:
                    repeat the __base_feature dictionary for each machine in 
                    the array
                machine position (dict)[-]:
                    contains the xy position ("position") of each machine 
                    ("turbineID") of the array.
        """
        
        turbines = {}
        features = {}
        
        if len(self.coord.shape) == 1:
            nb = 1
        else:
            nb = len(self.coord)
        
        n_digits = len(str(nb))
        
        for ii in range(nb):
            strn = 'turbine{:0{width}d}'.format(ii, width=n_digits)
            posxy = self.coord[ii, :]
            pxyz = np.array([posxy[0], posxy[1], self.hub_height])
            turbines[strn] =  {'position': pxyz}
            
            # Shallow copy is OK here as any mutables should be invariant
            # Worth checking if such features can be modified in the future.
            features[strn] = self.__base_feature.copy()
        
        return (nb, features, turbines)
    
    def __set_data(self, ss_id):
        """
        updates the __data dictionary for the specific sea state, specified in 
        the input argument (index)
        
        Args:
            ss_id (int)[-]: index of the current sea state.
        """
        nx = len(self.__data['X'])
        ny = len(self.__data['Y'])
        self.__data['U'] = row_major(nx, ny, self.__U[:,:,ss_id])
        self.__data['V'] = row_major(nx, ny, self.__V[:,:,ss_id])
        self.__data['TI'] = row_major(nx, ny, self.__TI[:,:,ss_id])
        self.__data['SSH'] = row_major(nx, ny, self.__SSH[:,:,ss_id])
    
    def __set_coordinates(self, coord):
        """
        Updates the coord argument for the specific array. If the 
        __mct_flag is True the number of machine is doubled and the 
        position of the machine is placed symmetrically wrt the original
        machine position.
        
        Args:
            coord (numpy.ndarray):
                UTM coordinates (Easting, Northing) of the  machines in the 
                current layout
            
        """
        
        if self.__mct_flag:
            self.coord = _get_perp_positions(self.__base_feature['OA'],
                                             self.__mct_dist,
                                             coord)
        else:
            self.coord = coord
        
        return


def _map_mct_results(pow_perf_dev_state,
                     pow_perf_dev_state_ni,
                     ti_dev_state,
                     turbines):
    """
    modifies the output generated by the tidal module, in order to account
    for the MCT case.
    
    Args:
        pow_perf_dev_state (numpy.ndarray)[W]:
            average power production of each device per each sea states 
            WITH interaction
        pow_perf_dev_state_ni (numpy.ndarray)[W]:
            average power production of each device  per each sea states
            WITHOUT interaction
        ti_dev_state (list)[-]:
            list of dictionaries gathering the information of the
            turbulence intensity of each device for each sea states.
        turbines (dict)[m]:
            dictionary gathering the information of the position of each
            machine in the array. For each key a tuple of UTM coordinates
            (Easting, Northing) is assigned.
    
    Returns:
        (tuple) (pow_perf_dev_state,
                 pow_perf_dev_state_ni,
                 ti_dev_state,
                 turbines,
                 nb)
            
            Same as input arguments but with machine coupled in pares (MCT)
    
    """
    (nb, n) = pow_perf_dev_state.shape
    nb /= 2
    ti_dev_state_n = []
    turbines_n = {}
    
    n_digits_old = len(str(nb * 2))
    n_digits_new = len(str(nb))
    
    for kk in range(nb):
        
        left_turb_name = 'turbine{:0{width}d}'.format(kk,
                                                      width=n_digits_old)
        right_turb_name = 'turbine{:0{width}d}'.format(kk + nb,
                                                       width=n_digits_old)
        turb_name = 'turbine{:0{width}d}'.format(kk, width=n_digits_new)
        
        pow_perf_dev_state[kk, ] = (pow_perf_dev_state[kk, ] +
                                              pow_perf_dev_state[kk+nb, ])
        pow_perf_dev_state_ni[kk, ] = (pow_perf_dev_state_ni[kk, ] +
                                         pow_perf_dev_state_ni[kk+nb, ])
        
        turb_pos = 0.5 * (turbines[left_turb_name]['position'] +
                          turbines[right_turb_name]['position'])
        
        turbines_n[turb_name] = {'position': turb_pos}
    
    for ll in range(n):
        
        local_dict = {}
        
        for kk in range(nb):
            
            left_turb_name = 'turbine{:0{width}d}'.format(
                                                    kk,
                                                    width=n_digits_old)
            right_turb_name = 'turbine{:0{width}d}'.format(
                                                    kk + nb,
                                                    width=n_digits_old)
            turb_name = 'turbine{:0{width}d}'.format(kk,
                                                     width=n_digits_new)
            
            turb_ti = 0.5 * (ti_dev_state[ll][left_turb_name] +
                             ti_dev_state[ll][right_turb_name])
            
            local_dict[turb_name] = turb_ti
        
        ti_dev_state_n.append(local_dict)
    
    ti_dev_state = ti_dev_state_n
    pow_perf_dev_state = np.delete(pow_perf_dev_state,
                                   np.s_[-nb:],
                                   0)
    pow_perf_dev_state_ni = np.delete(pow_perf_dev_state_ni,
                                      np.s_[-nb:],
                                      0)
    turbines = turbines_n
    
    return (pow_perf_dev_state,
            pow_perf_dev_state_ni,
            ti_dev_state,
            turbines,
            nb)


# Utilities
def row_major(nx, ny, q):
    """
    row_major: Converts the matrix q to row major if needed
    
    Args:
      - nx: x axis dimension, (int)
      - ny: y axis dimension, (int)
      - q: 2d matrix, (numpy.ndarray)
    Returns:
      - q: 2d matrix, dimension (ny, nx), (numpy.ndarray)
    """
    if q.shape == (nx, ny):
        q = q.T
    return q


def points_to_grid(points, xs, ys, fill_value=np.nan):
    
    #    The original reshaping did not take into account the posibility of
    #    non-rectangular domains.
    
    y_idx = get_indices(points[:, 1], ys)
    x_idx = get_indices(points[:, 0], xs)
    
    grid = np.empty(xs.shape + ys.shape)
    grid.fill(fill_value)
        
    grid[x_idx, y_idx] = points[:, 2]
              
    return grid


def get_indices(search, base):
    
    index = np.argsort(base)
    sorted_base = base[index]
    sorted_search = np.searchsorted(sorted_base, search)
    
    searchindex = np.take(index, sorted_search, mode="clip")
    mask = base[searchindex] != search
    
    result = np.ma.array(searchindex, mask=mask)
        
    return result


def _get_perp_positions(main_bearing, separation_len, centre_coord):
    
    perp_bearing = (main_bearing + 90) % 360
    perp_angle = bearing_to_radians(perp_bearing)
    perp_vector = np.array([np.cos(perp_angle), np.sin(perp_angle)])
    separation_vector = perp_vector * separation_len
    
    turb_coords = np.vstack((centre_coord + separation_vector,
                             centre_coord - separation_vector))
    
    return turb_coords
