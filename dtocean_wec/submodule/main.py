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
This module contains the main classes used to obtain the solution of the hydrodynamic problem for a given wave energy converter

.. module:: output
   :platform: Windows
   :synopsis: wec_external_module module to DTOcean

.. moduleauthor:: Francesco Ferri <ff@civil.aau.dk>
.. moduleauthor:: Mathew Topper <mathew.topper@dataonlygreater.com>
 """

import os

import numpy as np

import dtocean_wave.utils.hdf5_interface as h5i

from .nemoh_reader import NemohReader
from .wamit_reader import WamitReader
from .db_reader import DbReader
from .nemoh_run import NemohExecute
from .power_matrix_fitting import PowerMatrixFit
from .utils.write_results import write_to_raw_results


class BemSolution(object):
    """
    bem_solution: the class is used to sort the different input choices and generate a result folder
                with the structure used in the hydrodynamic module of the DTOcean tool.

    Args:
        input_type (int): identify the input type selected by the user.
            Four different input types are possible:
                1 (default) - precompiled_case selection. This option reload a solution specified via the ID from a WEC database shipped with the DTOcean tool.
                2 - data folder containing the inputs to run the embedded bem software Nemoh
                3 - data folder containing a pre-solved Nemoh problem
                4 - data folder containing a pre-solved WAMIT problem
         project_folder (string)[]: a project folder name. If not present the folder will be create.
         precompiled_case (int,optional)[]: number of the precompiled case to be selected from the WEC database.
         data_folder (string,optional)[]: path name of the folder holding the required input files.
         debug (boolean, string)[]: debug flag
    Attributes:
           same as input arguments

    Returns: None
    """
    def __init__(self,
                 dataobj,
                 db_folder,
                 bin_folder,
                 debug=False):
        
        self.debug = debug
        self.data_folder = dataobj.data_folder
        self.project_folder = dataobj.project_folder
        self.precompiled_case = dataobj.precompiled_case
        self.__get_array_matrix = dataobj.get_array_mat
        self.dataobj = dataobj
        self.db_folder = db_folder
        self.bin_folder = bin_folder

        self.input_type = dataobj.input_type
        self.reader = None
        


    def call_module(self):
        """
       
        """
        if self.input_type == 1:  # read data from wec database
            self.reader = DbReader(self.precompiled_case,
                                   self.db_folder,
                                   debug=self.debug)

        elif self.input_type == 2:  # solve the bem problem
            bem_obj = NemohExecute(self.project_folder,
                                   self.dataobj.general_inputs,
                                   self.dataobj.body_inputs,
                                   get_array_mat=self.__get_array_matrix,
                                   debug=self.debug)
                                   
            # generate the nemoh results
            bem_obj.gen_path()
            bem_obj.gen_mesh_files()
            bem_obj.gen_multibody_structure()
            bem_obj.gen_hdyn_files()
            bem_obj.run_nemoh(self.bin_folder)
            bem_obj.run_hst()
            
            self._bem_obj = bem_obj
            self.reader = NemohReader(  self.data_folder,
                                        self.dataobj.general_inputs,
                                        self.dataobj.body_inputs,
                                        get_array_mat=self.__get_array_matrix,
                                        debug=self.debug)

        elif self.input_type == 3:
            self.reader = NemohReader(  self.data_folder,
                                        self.dataobj.general_inputs,
                                        self.dataobj.body_inputs,
                                        get_array_mat=self.__get_array_matrix,
                                        debug=self.debug)

        else:
            self.reader = WamitReader(  self.data_folder,
                                        self.dataobj.general_inputs,
                                        get_array_mat=self.__get_array_matrix, 
                                        debug=self.debug)
                                    
        self.reader.load_data()
#        try:
#            self.reader.load_data()
#        except Exception as e:
#            raise ValueError(e)

    def pm_fitting(self,
                   machine_spec=None,
                   site_spec=None,
                   ext_k=None,
                   ext_d=None):
        """identifies additional damping and stiffness via an error minimisation between
        the certified power matrix and the automatic generated one.

        certified_machine (dict): dictionary containing the specification of the single machine.
            Dictionary keys:
                'power_matrix' (numpy.ndarray)[W]:
                'pto_damping' (numpy.ndarray)[N/(m/s) or Nm/(rad/s)]:
                'mooring_stiffness' (numpy.ndarray)[N/m or Nm/rad]:
                'scatter_diagrma' (numy.ndarray): Scatter diagram of the certified power matrix, defined through the probability of occurence of each sea state
                'specType' (tup): description of the spectral shape recorded at the site (spectrum name, gamma factor, spreading parameter)
                'te' (numpy.ndarray)[s]: Vector containing the wave energy periods
                'hm0' (numpy.ndarray)[m]: Vector containing the significant wave heights
                'dirs' (numpy.ndarray)[rad]: Vector containing the wave directions

        ext_k (numpy.ndarray): external/additional stiffness matrix to be applied in the numerical model of the WEC.
        ext_d (numpy.ndarray): external/additional damping matrix to be applied in the numerical model of the WEC.
        :return:
        """
        per_fit = PowerMatrixFit(self.reader.periods,
                                    self.reader.directions,
                                    self.reader.m_m,
                                    self.reader.m_add,
                                    self.reader.c_rad,
                                    ext_d,
                                    self.reader.k_hst,
                                    ext_k,
                                    self.reader.force_tr_mat,
                                    self.reader.pto_dof-1,
                                    self.reader.moor_dof-1,
                                    self.reader.order,
                                    self.db_folder,
                                    debug=self.debug)

        if self.input_type==1:
            per_fit.load_fitting_data(self.precompiled_case)
        else:
            if machine_spec == None or site_spec == None:
                raise IOError("The performance fitting cannot be run without the machine and site specification",
                              "The execution is aborted.")

            per_fit.perf_fitting(machine_spec, site_spec)

        self.per_fit = per_fit
        
        return

    def generate_outputs(self):
        # writing to a ASCII files is not really usefull.
        # TODO: talk with Mat and Vincenzo for the matter
        # for now just save it to a pickle
        #if not self.reader.n_dof == self.reader.loose_dof:
        # move the row data from Nemoh run into a raw_result folder
        # move the hydrodynamic data to the raw_results folder.
        # This applies only for the case 2 with nemoh run
        write_to_raw_results(self.project_folder)
        modes = self.reader.modes
       
        temp_dic = {'m_m': self.reader.m_m,
                    'm_add': self.reader.m_add,
                    'c_rad': self.reader.c_rad,
                    'f_ex': self.reader.f_ex,
                    'periods': self.reader.periods,
                    'directions': self.reader.directions,
                    'k_hst': self.reader.k_hst,
                    'diffraction_tr_mat': self.reader.diffraction_tr_mat,
                    'force_tr_mat': self.reader.force_tr_mat,
                    'amplitude_coefficient_radiation': self.reader.amplitude_coefficient_radiation,
                    'water_depth': self.reader.water_depth,
                    'cyl_radius': self.reader.cyl_radius,
                    'modes': modes,
                    'c_ext': self.per_fit.c_ext,
                    'k_ext': self.per_fit.k_ext,
                    'c_fit': self.per_fit.c_fit,
                    'k_fit': self.per_fit.k_fit,
                    'c_pto': self.per_fit.c_pto,
                    'k_mooring': self.per_fit.k_mooring,
                    'te': self.per_fit.te,
                    'hm0': self.per_fit.hm0,
                    'wave_dir': self.per_fit.wave_dir,
                    'scatter_diagram': self.per_fit.scatter_diagram,
                    'max_order': self.reader.order,
                    'truncation_order': self.reader.truncation_order,
                    'mooring_dof': self.reader.moor_dof,
                    'pto_dof': self.reader.pto_dof}

        h5i.save_dict_to_hdf5(temp_dic, os.path.join(self.project_folder, 'wec_solution.h5'))
        # h5pi.save_dict_to_pandas_and_hdf5(temp_dic, os.path.join(self.project_folder, 'test2.h5'))
        # pkli.save_dict_to_pickle(temp_dic, os.path.join(self.project_folder, 'test3.pkl'))


if __name__ == "__main__":
    # usage of the external wave module
    import sys
    sys.path.append(r"C:\Users\francesco\Desktop\test_gui\utils")
    
    from data_interface import DataStructure
    data = h5i.load_dict_from_hdf5(r"C:\Users\francesco\Desktop\test_gui\test_prj\test_prj_data_collection.hdf5")
    data['inputs_hydrodynamic']['general_input']['data_folder'] = r"C:\Users\francesco\Desktop\test_gui\mtatest\hydrodynamic"
    dataobj = DataStructure(data)
    if dataobj.check_inputs(): 
        wec_obj = BemSolution(dataobj,                
              debug=False)
        wec_obj.call_module()
    
    additional_stiffness = None
    additional_damping = None

    site_spec = {'spec_shape':'Jonswap'}
    site_spec['spec_gamma'] = 1.0
    site_spec['spec_spreading'] = -1
    site_spec['te'] = np.linspace(3,10,1)
    site_spec['hm0'] = np.linspace(0.5, 3.5, 1)
    site_spec['wave_angles'] = np.linspace(0,0,1)
    site_spec['probability_of_occurence'] = np.ones((1,1,1))

    machine_spec = {'c_pto': np.zeros((1,1,1,1,1))}
    machine_spec['yaw'] = 0.
    machine_spec['k_mooring'] = np.zeros((1,1,1,1,1))
    machine_spec['power_matrix'] = np.zeros((1,1,1))
    wec_obj.pm_fitting(machine_spec, site_spec,
                       ext_k=additional_stiffness,
                       ext_d=additional_damping)
#    wec_obj.generate_outputs()
#
#    # TODO: check the shape of the cpto and c_ext in input
#    # the visualiser will be hooked at the BEM object, which unifies the 4 possible cases
#    plotter = Visualiser(wec_obj)
#    plotter.show_rao(0,0,0,0,0)







