#!/usr/bin/python2.7
# encoding: utf-8
"""
This module contains the main class that performs the prepocess and 
calls the tidal and wave submodule.

.. module:: input
   :platform: Windows
   :synopsis: WP2 main module to DTOcean
   
.. moduleauthor:: Francesco Ferri <ff@civil.aau.dk>
"""
from __future__ import division

# Start logging
import logging
module_logger = logging.getLogger(__name__)

# External Package
import numpy as np

# Local package
from dtocean_wave.WEC import wec
from dtocean_wave.MultiWEChydro import MultiBody
from dtocean_tidal.interface import CallTidal
from dtocean_tidal.submodel.ParametricWake import read_database
from dtocean_hydro.utils.set_wdirs_multibody import set_wdirs_multibody
import dtocean_wave.utils.read_bem_solution as read_wec_sol

# Relative imports
from .hydro import Hydro_pkg
from .array import Array_pkg
from .output import WP2output
from .utils import optimiser
from .configure import get_install_paths


class WP2:
    """
    WP2: main classes of the DTOcean WP2 package. 
    The class is composed of other classes specifically 
    Hydro, Array, Input, WEC, MultiBody, CallTidal, optimiser.
    The inetent is to delegate the work to isolated objects, and the WP2 class 
    objective is to wrap the different objects.

    Args:
        WP2input (WP2input class): WP2 input class. 

    Optional args:
        Cfit (numpy.ndarray): Only for the wave case. User defined value of the
                              damping. Used in the PM fitting method. The input
                              can be used to enter an additional damping in the
                              numerical model
        Kfit (numpy.ndarray): Only for the wave case. User defined value of the
                              stiffness. Used in the PM fitting method. The
                              input can be used to enter an additional
                              stiffness in the numerical model

        NOTE: the Cfit and Kfit are obsolete, they are now specified in the
        wec_solution.h5 file

        debug (boolean): if set to True, plots and additional command line
                         outputs are issued.
     
    Attributes:
            bin_dir (str): path name of the bin. The directory stores the
                executable files for the BEM solver.
            include_dir (str): path name of the include directory. The
                directory stores the wave and tidal calculation and data
            tidal_include (str): path name of the tida directory, where the
                wake database is located.
            wave_include (str): path name of the wave directory, where the BEM
                solver results are located.
            iInput (WP2input class): copy of the input argument.
            iArray (Array_pkg class): contains the array related features as
                lease area, active area, and generates the spatial arrangement
                of the devices.
            iHydro (Hydro_pkg class): hydro class for the different arrays
            iHydroWEC (Hydro_pkg class): hydro class specific for the isolated
                WEC (only waves)
            _debug (boolean): if set to True, plots and additional comand line
                outputs are issued.
            iWEC (WEC class): (only for arrays of WECs). The iWEC class 
                containts the hydrodynamic features of the single body
                evaluated from the BEM solver.
            iMB (MultiBody class): (only for arrays of WECs ). The IMB class
                contains all the hydrodynamic features of the array, evaluated
                using the direct matrix method.
    
    Returns:
        WP2output (WP2output class): the class return an Output object. If no
        array layout is possible the returns is -1.
  
    """
    def __init__(self, WP2input,
                       Cfit=None,
                       Kfit=None,
                       pickup=False,
                       debug=False):
        
        # The input object is passed for use in the optimisation loop method
        self.iInput = WP2input

        # The array class assess the active area during the instantianisation
        self.iArray = Array_pkg(WP2input.S_data.LeaseArea,
                                WP2input.S_data._LeaseArea,
                                WP2input.M_data.MinDist,
                                WP2input.S_data.mainAngle,
                                WP2input.S_data.NogoAreas_bathymetry)
        self._debug = debug

        if not WP2input.internalOptim:
            module_logger.info("The user provided an external map of the "
                               "hydrodynamic interaction. No internal model "
                               "is issued")
        else:
            bathy = WP2input.S_data.Bathymetry

            if WP2input.M_data.tidalFlag:
                devType = 'T'
                V = WP2input.S_data.MeteoceanConditions['V']
                U = WP2input.S_data.MeteoceanConditions['U']
                p = WP2input.S_data.MeteoceanConditions['p']
                TI = WP2input.S_data.MeteoceanConditions['TI']
                x = WP2input.S_data.MeteoceanConditions['x']
                y = WP2input.S_data.MeteoceanConditions['y']
                SSH = WP2input.S_data.MeteoceanConditions['SSH']
                self.iHydro = Hydro_pkg(devType, V, U, p, TI, x, y, SSH, bathy)

            else:
                devType = 'W'
                freqs = read_wec_sol.read_freq(
                                            WP2input.M_data.wave_data_folder)
                cfreqs = 2.*np.pi*freqs
                B = WP2input.S_data.MeteoceanConditions['B']

                mA = WP2input.S_data.mainAngle
                yA = WP2input.M_data.YawAngle
                specType = WP2input.S_data.MeteoceanConditions['specType']
                wdirs = set_wdirs_multibody(B.copy(),
                                            specType[-1],
                                            mA,
                                            yA,
                                            debug=debug)

                B = WP2input.S_data.MeteoceanConditions['B']
                Hs = WP2input.S_data.MeteoceanConditions['Hs']
                Tp = WP2input.S_data.MeteoceanConditions['Tp']
                ScatDiag = WP2input.S_data.MeteoceanConditions['p']

                SSH = WP2input.S_data.MeteoceanConditions['SSH']
                wdepth = np.abs(self.iInput.S_data._average_water_depth + SSH)

                # Instance of Hydro for WEC
                #Tp_WEC = WP2input.M_data.SingleDevice['x']
                #Hs_WEC = WP2input.M_data.SingleDevice['y']
                #B_WEC = WP2input.M_data.SingleDevice['z']
                #ScatDiag_WEC = WP2input.M_data.SingleDevice['Add']['p']
                #specType_WEC = WP2input.M_data.SingleDevice['Add']['specType']
                #self.iHydroWEC= Hydro_pkg(devType, B_WEC, Hs_WEC, Tp_WEC,
                #     ScatDiag_WEC, wdepth, freqs, wdirs_single, specType_WEC)
                # Instance of Hydro for array (the only difference is wdirs)
                self.iHydroMB = Hydro_pkg(devType,
                                          B,
                                          Hs,
                                          Tp,
                                          ScatDiag,
                                          wdepth,
                                          cfreqs,
                                          wdirs,
                                          specType)

            # for the wave case the single body model need created from the
            # hydrodynamic data and then compared with the user input power
            # matrix
            if WP2input.M_data.tidalFlag:
                # Load the tidal database
                if WP2input.M_data.tidal_data_folder is None:
                    # Get path to data dirs through configuration
                    path_dict = get_install_paths()
                    self.cfd_data = read_database(path_dict["tidal_include"])
                else:
                    self.cfd_data = read_database(
                                        WP2input.M_data.tidal_data_folder)
            else:
                self.iWEC = wec(WP2input.M_data.wave_data_folder,
                                wdepth,
                                self.iInput.M_data.RatedPowerDevice,
                                debug=False)
                self.iWEC.load_single_machine_model(
                        WP2input.S_data.MeteoceanConditions['specType'])

                # for the wave case if any of the array distance is smalled
                # than the circumscribing cylinder the implemented theory is
                # violated and the results are unrealistic. The Dmin inputed by
                # the user is therefore  updated
                if min(self.iArray.Dmin)/2. < self.iWEC.radius:
                    Dmin_new = ()
                    for D in self.iArray.Dmin:
                        Dmin_new += (max(D/2,self.iWEC.radius*1.001)*2,)
                        module_logger.warning("The minimum distance specified "
                                              "by the user does not fit the "
                                              "requirements of the interacton "
                                              "theory. The minimum distance "
                                              "is modified.")
                        module_logger.info(
                                    "Old input: {}".format(self.iArray.Dmin))
                        module_logger.info(
                                    "Interaction theory min requirement: "
                                    "{}".format(self.iWEC.radius*2))
                        module_logger.info(
                                    "New min distance: {}".format(Dmin_new))
                    self.iArray.Dmin = Dmin_new

    def optimiseExternalTab(self):
        """
        optimiseExternalTab: reads the array layouts given in the dictionary
        and return the optimal configuration under the given constraints.
        
        Returns: 
            Ouput obj
        """
        userIn = self.iInput.M_data.UserOutputTable.copy()
        q_min = self.iInput.M_data.OptThreshold
        AEP_arr = np.zeros((len(userIn.keys())))
        q_arr = np.zeros((len(userIn.keys())))

        for conf in userIn.keys():
            ID = int(conf[4:])
            q_arr[ID] = userIn[conf]["q_factor_Array"]
            layout = userIn[conf]["Array_layout"]
            Nb = userIn[conf]["Nbodies"]
            coord = np.zeros((Nb,2))
            for key in layout.keys():
                IDb = int(key[6:])
                coord[IDb,0] = layout[key][0]
                coord[IDb,1] = layout[key][1]

            self.iArray.coord = coord
            self.iArray.checkMinDist()

            if self._debug:
                inside = self.iArray.checkout(
                                    nogo_list=self.iInput.S_data.NogoAreas)
                self.iArray.show(inside)
            AEP = userIn[conf]["Annual_Energy_Production_perD"]

            # it has been decided to do not remove the unfeasible nodes
            # this step is left to the user, it only used in the debug mode
            # to visualise the coherence between the user grid and the
            # actual feasible area
            # the only operation performed is the check on the minimum distance
            if not self.iArray.minDist_constraint:
                AEP_arr[ID] = AEP.sum()
            # if not np.any(inside):
            #    AEP_arr[ID] = -1
            # else:
            #    AEP_arr[ID] = AEP[inside].sum()

        mask_q = q_arr >= q_min
        if np.any(mask_q):
            indM = np.argmax(AEP_arr[mask_q])
            IDs = np.arange(len(userIn.keys()))[mask_q]
    
            OptConf = userIn["conf{}".format(IDs[indM])]

            return WP2output(OptConf["Annual_Energy_Production_Array"],
                        OptConf["power_prod_perD_perS"],
                        OptConf["Annual_Energy_Production_perD"],
                        OptConf["power_prod_perD"],
                        OptConf["Array_layout"],
                        OptConf["Nbodies"],
                        OptConf["Resource_Reduction"],
                        OptConf["Hydrodynamic_Parameters"],
                        OptConf["q_factor_Per_Device"],
                        OptConf["q_factor_Array"],
                        OptConf["TI"])
        else:
            raise ValueError("No feasible configurations found in the given "
                             "result table")

    def optimisationLoop(self):
        """
        optimisationLoop: the method calls iteratively either the tidal or wave
            sub-models, in order to identify the optimal array configuration
            that satisfy the given constraints. The optimisation criteria is
            the maximisation of the annual energy production of the whole
            array.
        
        Note:
            if the attribute self.internalOptim is False, the WP2 is using a
            precompiled list of outputs dictionaries and no internal
            calculation is performed. This option can be actually used to run
            a batch run and speed up the overall DTOcean optimisation.
            
            Depending on the self.iInput.M_data.UserArray['Option'] different
            behavior of the method are achieved:
                Opt 1 --> Normal use, optimisation based on the prebuilted
                          parameterised array layouts
                Opt 2 --> Fixed array layout, no optimisation is done, only the
                          evaluation of the specific array layout performance
                          is issued.
                Opt 3 --> User defined array layout, optimisation of the given
                          layout, by stretching/compressing of the device
                          position.
            
        """
        
        warning_str = (
               'The given BEM solution for the isolated device is calculated '
               'at a water depth of {}m while the average bathymetry of the '
               'lease area minus the input nogo-zones is {}m and the average '
               'installation depth of the different machines with the given '
               'array layout is {}m. The BEM software should be run with a '
               'water depth close to the average water depth of the lease '
               'area to minimise the errors introduced with the flat '
               'bathymetry approximation. With the actual hydrodynamic model '
               'the accuracy of the results is reduced. Consider '
               're-evaluating the hydrodynamic model of the isolated device '
               'for the actual average bathymetry [{}m]. This value is '
               'calculated by averaging the bathymetry gridpoints that satify '
               'the following three conditions:\n'
               '\t --> point inside the lease area\n'
               '\t --> point outside the user-given nogo-zones\n'
               '\t --> point below the considered water level datum.')

        if not self.iInput.internalOptim:
            return self.optimiseExternalTab()

        stopRun = False
        Opt = self.iInput.M_data.UserArray['Option']
        Value = self.iInput.M_data.UserArray['Value']

        # initialise either the tidal or wave object
        if self.iInput.M_data.tidalFlag:
            hyd_obj = CallTidal(self.iHydro,
                                self.iInput,
                                self.cfd_data,
                                debug=self._debug,
                                debug_plot=self._debug)
        else:
            hyd_obj = MultiBody(self.iHydroMB,
                                self.iWEC,
                                cylamplitude=True)

        if Opt == 2:
            self.iArray.coord = Value
            self.iArray.checkMinDist()
        else:
            opt_obj = optimiser.SearchOptimum(hyd_obj,
                                              self.iArray,
                                              Value, Opt,
                                              self.iInput.M_data.MaxNumDevices,
                                              self.iInput.M_data.OptThreshold,
                                              self.iInput.S_data.NogoAreas,
                                              debug=False)

            stat = opt_obj.eval_optimal_layout(opt_method_id=1)
            if stat < 0:
                stopRun = True
            # regenerate the optimal array layout

        if not stopRun:
            
            module_logger.info('Finishing the WP2 task: Evaluation of the '
                               'final array layout interaction')
            inside = self.iArray.checkout(
                                        nogo_list=self.iInput.S_data.NogoAreas)
            if self._debug: self.iArray.show(inside)
            if not np.any(inside):
                return -1
            
            if Opt == 2 and not inside.all():
                exc_strings = ["({}, {})".format(xy[0], xy[1])
                                        for xy in self.iArray.coord[~inside]]
                exc_string = ", ".join(exc_strings)
                infoStr = ("Devices at positions {} have been excluded. "
                           "Check lease area boundary and no-go "
                           "areas").format(exc_string)
                module_logger.info(infoStr)

            if not self.iInput.M_data.tidalFlag:
                
                bem_depth = self.iWEC.water_depth
                site_depth = self.iWEC.depth
                layout = self.iArray.coord[inside]
                SSH = np.nanmean(self.iInput.S_data.MeteoceanConditions['SSH'])
                    
                dev_depth = get_device_depths(self.iInput.S_data.Bathymetry,
                                              layout)
                device_average_depth = -np.mean(dev_depth) + SSH
                                              
                if (np.any(np.logical_or(bem_depth <= site_depth * 0.99,
                                         bem_depth >= site_depth * 1.01)) or 
                    np.any(np.logical_or(bem_depth <=
                                                 device_average_depth * 0.99,
                                         bem_depth >=
                                                 device_average_depth * 1.01))
                    ):
                    
                    module_logger.warning(warning_str.format(
                                                        bem_depth,
                                                        site_depth,
                                                        device_average_depth,
                                                        site_depth))

                
            # solve the hydrodynamic interaction of the array
            hyd_res = hyd_obj.energy(self.iArray.coord[inside])
            norm_dir = (self.iInput.S_data.Main_Direction / 
                            np.linalg.norm(self.iInput.S_data.Main_Direction))
            
            power_matrix_dims = {
                       "te": self.iInput.S_data.MeteoceanConditions['Te'],
                       "hm0": self.iInput.S_data.MeteoceanConditions['Hs'],
                       "dirs":  self.iInput.S_data.MeteoceanConditions['B']}

            result = WP2output(hyd_res.AEP_array,
                               hyd_res.power_prod_perD_perS,
                               hyd_res.AEP_perD,
                               hyd_res.power_prod_perD,
                               hyd_res.Device_Position,
                               hyd_res.Nbodies,
                               hyd_res.Resource_reduction,
                               hyd_res.Device_Model,
                               hyd_res.q_perD,
                               hyd_res.q_array,
                               norm_dir,
                               hyd_res.TI,
                               hyd_res.power_matrix_machine,
                               power_matrix_dims)

            result.remap_res(self.iInput.S_data.electrical_connection_point)
            result.logRes()
            
        else:
            
            result = -1
            
        return result
    

def get_device_depths(bathymetry, layout):
    
    if len(layout.shape) != 2:
        raise IndexError("Layout array does not have two dimensions")
        
    if layout.shape[1] != 2:
        raise IndexError("Layout array second dimension does not have two "
                         "elements")
    
    xr = bathymetry[:,0]
    yr = bathymetry[:,1]
    zr = bathymetry[:,2]
    
    # Get absolute tolerance for x
    uniquex = np.unique(xr.round(decimals=4))                
    dx = uniquex[1] - uniquex[0]
    xtol = 0.1 * dx
        
    dev_depth = []
    
    for ndev in range(layout.shape[0]):
        
        # Find the nearest matching x value
        idx = (np.abs(xr - layout[ndev][0])).argmin()
        x_near = xr[idx]
        
        # Find indexes with matching xnear and filter yr & zr
        xmap = np.isclose(xr, x_near, rtol=0, atol=xtol)
        yrx = yr[xmap]
        zrx = zr[xmap]
        
        # Find nearest matching y and the depth
        idx = (np.abs(yrx - layout[ndev][1])).argmin()
        ndev_depth = zrx[idx]
        
        dev_depth.append(ndev_depth)
        
    return dev_depth
