# -*- coding: utf-8 -*-

#    Copyright (C) 2016 Pau Mercadez Ruiz
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
This module contains the class used to generate the numerical model of the WEC
using a given WAMIT solution.

.. module:: hyd_WAMIT
   :platform: Windows
   :synopsis: Numerical model of WEC builder

.. moduleauthor:: Pau Mercadez Ruiz <pmr@civil.aau.dk>
"""

import os

import numpy as np

from hydrostatics import Hydrostatics_Nemohcal
import utils.input_control as checker
import utils.file_utilities as f_util
from utils.transfers import transfers
from utils.multibody_analysis import MultiBodyAnalysis
from utils.general import block_diag


class NemohReader():
    def __init__(self,
                 data_folder,
                 general_inputs,
                 body_inputs,
                 get_array_mat=True,
                 debug=False):
        """
        self.phi_s = None
        self.phi_r = None
        self.f_ex = None
        self.m_add = None
        self.c_rad = None
        self.k_hst = None
        self.m_m = None
        self.cyl_r = None
        self.cyl_z = None
        self.cyl_t = None
        self.directions = None
        self.periods = None
        self.n_dof = None
        self.modes = None
        self.pto_dof = None
        self.moor_dof = None

        """
        self.rho = 1025
        self.g = 9.82
        self._debug = debug
        self.__get_array_mat = get_array_mat
        self.g_d = general_inputs
        self.b_d = body_inputs
        self.water_depth = self.g_d['water_depth']
        self.cyl_nzeta = self.g_d['cyl_nzeta']
        self.cyl_ntheta = self.g_d['cyl_ntheta']
        self.tot_dof = self.g_d['ndof']
        self.pto_dof = self.g_d['pto_dof']
        self.moor_dof = self.g_d['mooring_dof']
        self.point_application = self.b_d['local_cs']
        self.frequency_def = self.g_d['frequency_def']
        self.angle_def = self.g_d['angle_def']
        
        self.shared_dof_binary = self.b_d['shared_dof']
        
        bodies = self.b_d['body'].values()  # list of bodies
        self.bodies = bodies
        self.n_bodies = len(bodies)

        self._path_prj_hdy = data_folder
        self._path_prj_dyn_res = os.path.join(self._path_prj_hdy, 'results')
        
        
        self.periods = 2*np.pi/(np.linspace(self.frequency_def[1],self.frequency_def[2],self.frequency_def[0]))
        if self.angle_def[1]== 0. and self.angle_def[2] == 360.0 :
            self.directions = np.linspace(self.angle_def[1],self.angle_def[2],self.angle_def[0], endpoint=False)
        else:
            self.directions = np.linspace(self.angle_def[1],self.angle_def[2],self.angle_def[0])
            
        self.cyl_t = None
        self.cyl_z = None
        
        self.loose_dof = len(self.bodies)*len(np.where(self.shared_dof_binary)[0])
        
  
        self.m_m = None                               
        self.vel_map = None
        self.m_add = None 
        self.c_rad = None 
        self.f_ex = None
        self.k_hst = None
        self.diffraction_tr_mat = None
        self.force_tr_mat = None
        self.amplitude_coefficient_radiation = None
        self.order = None
        self.truncation_order = None
        

     
        self.cyl_radius = None
       
        self.n_dof = None
        self.modes = None

    def check_inputs(self):
        """check_inputs perform a control of the data structure for the case where
        a nemoh result is provided or where the bem problem need to be solved

        input_res: true--> results are given
        input_res_ false--> results need to be calculated
        """

        status = checker.check_nemoh_results(self._path_prj_hdy)

        if not status[0]:
            raise ValueError(status[1])

        return 0

    def load_data(self):
        
        mb_obj = MultiBodyAnalysis(self.bodies, self.shared_dof_binary, self.point_application)
        mb_obj.set_multibody_data()
        mb_obj.get_multibody_topology()
        mb_obj.set_joints_position_orientation()
        mb_obj.eval_velocity_transformation_matrix()
        bodies_dofs = mb_obj.get_dofs_floating_bodies()
        self.vel_map = mb_obj.b_mat
 
        with open(os.path.join(self._path_prj_hdy,'nemoh.cal'), 'r') as fhand:
            self.water_depth = float(fhand.readlines()[3].split('!')[0])
            
        self.cyl_radius = self.get_cylinder_dimension()    
        self.modes = bodies_dofs[0]
        # at this point either the solution is provided or calculated
        # the project folder is populated with the standard structure
        (phi_s, phi_r) = self.read_results()
        

        # calculate diffraction and force transfer matrix
        if self.__get_array_mat: 
            (self.diffraction_tr_mat,
            self.force_tr_mat,
            self.amplitude_coefficient_radiation,
            self.order,
            self.truncation_order) = transfers(self.water_depth,
                                      self.directions*np.pi/180.,
                                      self.periods,
                                      (self.cyl_radius, self.cyl_t, self.cyl_z),
                                      phi_s,
                                      phi_r,
                                      self.f_ex)
        else:
            (self.diffraction_tr_mat,
            self.force_tr_mat,
            self.amplitude_coefficient_radiation,
            self.order,
            self.truncation_order) = np.array([0]), np.array([0]), np.array([0]), np.array([0]), np.array([0])
        
        # the multibody is added to the reader for debugging and verification porpouses
        # TODO: remove this in the final version
        self._multibody_object = mb_obj
        del phi_s, phi_r
    

        
    def map_generalise_dof(self, k_hst, c_rad, m_add, f_ex, phi_r):
        B = self.vel_map
        self.m_add = B.T.dot(m_add.transpose(2,1,0)).transpose(2,0,1).dot(B)
        self.c_rad = B.T.dot(c_rad.transpose(2,1,0)).transpose(2,0,1).dot(B)
        self.f_ex = B.T.dot(f_ex.transpose(1,2,0)).transpose(2,1,0)
        self.k_hst = B.T.dot(k_hst).dot(B)   
        
        # generate a block mass matrix
        M = np.zeros((0))
        for el in range(self.n_bodies):
            mass = np.zeros((6,6))
            mass[[0,1,2],[0,1,2]] = self.bodies[el]['mass']
            mass[3:,3:] = self.bodies[el]['inertia']
            if el == 0:
                M = block_diag( mass[np.ix_(np.where(self.shared_dof_binary)[0],
                                            np.where(self.shared_dof_binary)[0])] )
            else:
                M = block_diag(M, mass[np.ix_(np.where(self.shared_dof_binary)[0],
                                              np.where(self.shared_dof_binary)[0])] )
                    
        self.m_m = B.T.dot(M).dot(B)
        
        if not self.__get_array_mat:
            phi_r = np.zeros(c_rad.shape[:-1]+(1,1))
            
        
        phi_r_m = B.T.dot(phi_r.transpose(0,2,1,3)).transpose(1,0,2,3)
         
        return phi_r_m


    def get_cylinder_dimension(self):
        with open(os.path.join(self._path_prj_hdy,'nemoh.cal'),'r') as fid:
            lines = fid.readlines()

        # TODO: this is a weak point, increase the realiability of this search
        for i_lin in range(len(lines)):
            if not lines[-i_lin-1]=='\n':
                break

        cylinder = f_util.split_string(lines[-i_lin-2])
        cyl_nzeta = int(cylinder[2])
        cyl_ntheta = int(cylinder[1])
        cyl_radius = cylinder[0]
        
        return cyl_radius


    def read_results(self):
        
        k_hst = self.__read_hst()

        (freq_ca, ca) = self.__read_ca()
        (freq_cm, cm) = self.__read_cm()

        (freq_fex, angle_fex, fex) = self.__read_fex()

        freq_ca, freq_cm, freq_fex, angle_fex = [vect.reshape(-1) for vect in [freq_ca, freq_cm, freq_fex, angle_fex]]

        # perform some final check on the data format
        if not np.allclose(freq_ca,freq_cm) or not np.allclose(freq_cm,freq_fex,rtol=1e-3):
            raise AttributeError("The frequency retrieved from the nemoh results does not match each other.",
                                 "Check the consistency of the result folder")

        if not np.allclose(freq_ca,2.*np.pi/self.periods, rtol=1e-3):
            raise ValueError("The frequency vector retrieved from the nemoh results does not match the one",
                             " given in the input.")

        if not np.allclose(angle_fex,self.directions, rtol=1e-3):
            raise ValueError("The angle vector retrieved from the nemoh results does not match the one",
                             " given in the input.")

        # read the heavy files
        (phi_s, phi_rad, self.cyl_t, self.cyl_z) = self.__read_cylsurface()

        # map the loose dofs into the generalised dofs space
        phi_r = self.map_generalise_dof(k_hst, ca, cm, fex, phi_rad)
        
        del phi_rad
        return (phi_s, phi_r)

    def __read_hst(self):
        directory = os.path.join(self._path_prj_hdy)
        k_hst = Hydrostatics_Nemohcal(directory)
        
        return k_hst

    def __read_ca(self):
        f_ca = open(os.path.join(self._path_prj_dyn_res,'CA.dat'), 'r')
        n_freq = int(float(f_ca.readline().split(':')[-1]))
        ca = np.zeros((n_freq, self.loose_dof, self.loose_dof))
        freq = np.zeros((n_freq,1))  # Hz

        # Nemoh split the line at the sixth element.
        # Therefore a system with 9 dof will have 1full line plus 3elements in the second row.
        # These last are still relative to the same i-index.
        # The number or row per degree of freedom are totDOF/6
        n_lines = int(np.ceil(self.loose_dof/6.))

        for i_freq in range(n_freq):
            freq[i_freq] = float(f_ca.readline().split()[0])
            for dof in range(self.loose_dof):
                ls_tmp = []
                for nls in range(n_lines):
                    ls_tmp += f_ca.readline().split()

                ca[i_freq,dof,:] = np.asarray(ls_tmp, 'f')

        f_ca.close()

        return freq, ca

    def __read_cm(self):
        f_cm = open(os.path.join(self._path_prj_dyn_res,'CM.dat'), 'r')
        n_freq = int(float(f_cm.readline().split(':')[-1]))
        cm = np.zeros((n_freq, self.loose_dof, self.loose_dof))
        freq = np.zeros((n_freq,1))  # Hz

        # Nemoh split the line at the sixth element.
        # Therefore a system with 9 dof will have 1full line plus 3elements in the second row.
        # These last are still relative to the same i-index.
        # The number or row per degree of freedom are totDOF/6
        n_lines = int(np.ceil(self.loose_dof/6.))

        for i_freq in range(n_freq):
            freq[i_freq] = float(f_cm.readline().split()[0])
            for dof in range(self.loose_dof):
                ls_tmp = []
                for nls in range(n_lines):
                    ls_tmp += f_cm.readline().split()

                cm[i_freq,dof,:] = np.asarray(ls_tmp, 'f')

        f_cm.close()

        return freq, cm

    def __read_fex(self):
        f_Ex = open(os.path.join(self._path_prj_dyn_res,'ExcitationForce.tec'), 'r')
        n_directions = int(self.angle_def[0])
        n_freq = int(self.frequency_def[0])

        fex = 1j*np.zeros((n_freq,n_directions,self.loose_dof))
        angles = np.zeros((n_directions,1))
        frequencies = np.zeros((n_freq,1))
        n_lines = int(np.ceil(((self.loose_dof*2)+1)/80.))
        for skip_header in range(self.loose_dof+1):
            f_Ex.readline()

        for i_angle in range(n_directions):
            line = f_Ex.readline()
            angles[i_angle] = float(line.split("=")[2].split()[0])
            for i_freq in range(n_freq):
                ls_tmp = []
                for nls in range(n_lines):
                    ls_tmp += f_Ex.readline().split()

                frequencies[i_freq] = float(ls_tmp[0])
                MagPha = np.asarray(ls_tmp[1:], 'f')
                fex[i_freq,i_angle,:] = MagPha[0:-1:2]*np.exp(1j*MagPha[1::2])

        f_Ex.close()

        return (frequencies, angles, np.conj(fex))

    def __read_cylsurface(self):
        if not self.__get_array_mat:        
            return (0, 0, 0, 0)
        n_freq = len(self.periods)
        n_directions = len(self.directions)

        if not self.cyl_radius == 0.:

            n_problems = self.loose_dof*n_freq+n_directions*n_freq

            Pressure = 1j*np.zeros((self.cyl_ntheta, self.cyl_nzeta, n_problems))
            Location = np.zeros((self.cyl_ntheta, self.cyl_nzeta, 3))

            for problems in range(n_problems):
                #print problems
                f_prob = open(os.path.join(self._path_prj_dyn_res,'cylsurface.{:5d}.dat'.format(problems+1)), 'r')
                f_prob.readline()
                f_prob.readline()

                for n_zeta in range(self.cyl_nzeta):
                    for n_theta in range(self.cyl_ntheta):
                        temp = np.asarray(f_prob.readline().split(), 'f')
                        Location[n_theta, n_zeta, :] = temp[0:3]
                        Pressure[n_theta, n_zeta, problems] = temp[3]*np.exp(1j*temp[4])
                f_prob.close()

            Scattering = 1j*np.zeros((n_freq, n_directions, self.cyl_nzeta, self.cyl_ntheta))
            Radiation =  1j*np.zeros((n_freq, self.loose_dof, self.cyl_nzeta, self.cyl_ntheta))

            problem_index = -1
            for i_cfreq in range(n_freq):
                for i_angle in range(n_directions):
                    problem_index +=1
                    Scattering[i_cfreq, i_angle,:,:] = np.transpose(Pressure[:,:, problem_index])
                for dof in range(self.loose_dof):
                    problem_index +=1
                    Radiation[i_cfreq,dof,:,:] = np.transpose(Pressure[:,:, problem_index])
        else:
            # Dummy data when no field points are selected
            raise ValueError("The get_array_matrix flag is true but the cylinder radius is zero. ", 
                             "The execution will be terminated")
            #Scattering = 1j*np.zeros((n_freq,n_directions,1,1))
            #Radiation =  1j*np.zeros((n_freq,self.loose_dof,1,1))
            #Location = np.zeros((1,1,3))

        # Cast nan values from Scattering and Radiation to 0
        Scattering[np.isnan(Scattering)] = 0.
        Radiation[np.isnan(Radiation)] = 0.

        cfreq = 2.*np.pi/self.periods
        g = self.g
        for ind, freq in enumerate(cfreq):
            Scattering[ind] *= g/1j/freq
            Radiation[ind] *= g/1j/freq*(-1j*freq) # *(-1j*freq) Radiation will be due to unit amplitude motion

        Scattering = np.conj(Scattering)
        Radiation = np.conj(Radiation)

        # cylindrical discretization
        cyl_t = np.linspace(0, 2*np.pi, self.cyl_ntheta, endpoint = False)
        cyl_z = Location[0,:,2]

        return (Scattering, Radiation, cyl_t, cyl_z)


if __name__ == "__main__":
    import sys
    sys.path.append(r"C:\Users\francesco\Desktop\test_gui\utils")
    from data_interface import DataStructure
    import dtocean_wave.utils.hdf5_interface as h5i
    data = h5i.load_dict_from_hdf5(r"C:\Users\francesco\Desktop\test_gui\test_prj\test_prj_data_collection.hdf5")
    data_path = r"C:\Users\francesco\Desktop\test_gui\test_prj\hydrodynamic"
    dataobj = DataStructure(data)
    dataobj.body_inputs['body']['body0']['mesh'] = os.path.join(u'C:\\Users\\francesco\\Desktop\\test_gui', dataobj.body_inputs['body']['body0']['mesh'])
    dataobj.body_inputs['body']['body1']['mesh'] = os.path.join(u'C:\\Users\\francesco\\Desktop\\test_gui', dataobj.body_inputs['body']['body1']['mesh'])

    reader = NemohReader(dataobj.project_folder,
                         data_path,
                          dataobj.general_inputs,
                          dataobj.body_inputs,
                          get_array_mat=False,
                          debug=False)
    #reader.check_inputs()
    reader.load_data()
    
