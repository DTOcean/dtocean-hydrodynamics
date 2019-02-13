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
This module contains the main classes used to obtain the solution of the hydrodynamic problem for a given wave energy converter

.. module:: output
   :platform: Windows
   :synopsis: wec_external_module module to DTOcean

.. moduleauthor:: Francesco Ferri <ff@civil.aau.dk>
.. moduleauthor:: Mathew Topper <mathew.topper@dataonlygreater.com>
"""
import os

import numpy as np
from scipy import optimize

import dtocean_wave.utils.hdf5_interface as h5i
from dtocean_wave.utils.StrDyn import EnergyProduction
from dtocean_wave.utils.WatWaves import len2

from utils.conversion_utils import *


class PowerMatrixFit():
    def __init__(self, periods,
                       directions,
                       m_m,
                       m_add,
                       c_rad,
                       c_ext,
                       k_hst,
                       k_ext,
                       f_ex,
                       force_matrix,
                       pto_dof,
                       moor_dof,
                       order,
                       db_folder,
                       debug=False):
        self.m_m  = m_m
        self.m_add = m_add
        self.c_rad = c_rad
        if c_ext is None:
            c_ext = 0.
        self.c_ext = c_ext
        if k_ext is None:
            k_ext = 0.
        self.k_ext = k_ext
        self.k_hst = k_hst
        self.force_matrix = force_matrix
        self.n_dof = m_add.shape[-1]
        self.periods = periods
        self.directions = directions
        self.pto_dof = pto_dof
        self.moor_dof = moor_dof
        self.order = order
        self.db_folder = db_folder
        self.debug = debug
        self._fex = f_ex

        self.c_pto = None
        self.k_mooring = None
        self.c_fit = None
        self.k_fit = None
        self.tp = None
        self.te = None
        self.hm0 = None
        self.wave_dir = None
        self.scatter_diagram = None
        self.wec_power_matrix = None
        self.wec_original_power_matrix = None

    def load_fitting_data(self, caseID):
        f_n = os.path.join(self.db_folder,
                           'case{}'.format(caseID),
                           'case{}_pmfit.h5'.format(caseID))
        dic = h5i.load_dict_from_hdf5(f_n)
        self.c_ext = dic['c_ext']
        self.k_ext = dic['k_ext']
        self.c_fit = dic['c_fit']
        self.k_fit = dic['k_fit']
        self.c_pto = dic['c_pto']
        self.k_mooring = dic['k_mooring']
        self.te = dic['te']
        self.tp = convert_te2tp(self.te.copy())
        self.hm0 = dic['hm0']
        self.wave_dir = convert_angle(dic['wave_dir'].copy())
        self.scatter_diagram = dic['scatter_diagram']

    def __check_matrix_dimension(self, mat, sd_shape, application_dof):
        ndof = self.n_dof
        if isinstance(mat, float):  # check for a simple case
            mat_diag = np.zeros((ndof))
            mat_diag[application_dof] = 1.
            sea_state_mat = np.diag(mat_diag)
            out_matrix = np.resize(sea_state_mat,sd_shape+(ndof,ndof))
            out_matrix *= mat
        elif type(mat) == np.ndarray:
            if mat.shape == (ndof,ndof):
                out_matrix = np.resize(mat,sd_shape+(ndof,ndof))
            elif mat.shape  == sd_shape + (ndof, ndof):
                out_matrix = mat
            else:
                raise ValueError("Data shape not understood")
        else:
            raise ValueError("Data type not understood")

        return out_matrix

    def perf_fitting(self, machine_spec, site_spec ):
        """
        PerfFitting: the method is used to minimise the error between the internally generated numerical model and
                        certified power matrix of the isolated WEC.
                        The error is minimised by adding damping and stiffness to the system which values is obtained
                        via a bounded optimisation problem.

        Args:
            single_device (dic): dictionary gathering the information of the single machine

        Optional args:
            Cfit (numpy.ndarray): damping fitting specified as input to the method
            Kfit (numpy.ndarray): stiffness fitting specified as input to the methos
            pickup (boolean): pickup the solution from a pickle file if available

        """
        sp_distr = site_spec['probability_of_occurence']
        sp_type = site_spec['spec_shape']
        sp_gamma = site_spec['spec_gamma']
        sp_spreading = site_spec['spec_spreading']
        te = site_spec['te']
        tp = convert_te2tp(te.copy(), sp_type, sp_gamma)
        hm0 = site_spec['hm0']
        wave_angles = convert_angle(site_spec['wave_angles'].copy())

        y_angle_range = machine_spec['yaw']
        m_angle = wave_angles[sp_distr.sum((0,1)).argmax()]

        wdirs_yaw = set_wdirs_with_yaw(wave_angles.copy(),
                                       sp_spreading,
                                       m_angle,
                                       y_angle_range)

        self.tp = tp
        self.te = te
        self.hm0 = hm0
        self.wave_dir = wave_angles
        self.scatter_diagram = sp_distr
        self.scatter_diagram_spec = (sp_type, sp_gamma, sp_spreading)
        sd_shape = sp_distr.shape

        ndof = self.n_dof
        periods = self.periods  # wave periods analysed in bem

        self.c_ext = self.__check_matrix_dimension(self.c_ext, sd_shape, range(ndof))
        self.k_ext = self.__check_matrix_dimension(self.k_ext, sd_shape, range(ndof))
        CPTO_sm = self.__check_matrix_dimension(machine_spec['c_pto'], sd_shape, self.pto_dof)  # pto damping matrix for the single machine
        # overrule the user input by setting the index priority to the one read in the gnr file
        mask = np.ones_like(CPTO_sm, dtype=bool)
        mask[:,:,:,self.pto_dof, self.pto_dof] = False
        CPTO_sm[mask] = 0.

        Kmooring_sm = self.__check_matrix_dimension(machine_spec['k_mooring'], sd_shape, self.moor_dof )  # mooring stiffness damping matrix for the single machine
        # overrule the user input by setting the index priority to the one read in the gnr file
        mask = np.ones_like(Kmooring_sm, dtype=bool)
        mask[:,:,:,self.moor_dof, self.moor_dof] = False
        Kmooring_sm[mask] = 0.

        self.c_pto = CPTO_sm
        self.k_mooring = Kmooring_sm
        Kfit_sm = np.zeros(Kmooring_sm.shape)
        Cfit_sm = np.zeros(CPTO_sm.shape)
        wec_power_matrix = np.zeros(sd_shape)
        
        if np.all(machine_spec['power_matrix']==0):
            self.c_fit = Cfit_sm
            self.k_fit = Kfit_sm
            if self.force_matrix.shape == (1L,):
                self.wec_power_matrix = self.__get_wec_power_matrix(wdirs_yaw,
                                                                wec_power_matrix,
                                                                no_yaw=True) 
            else:
                self.wec_power_matrix = self.__get_wec_power_matrix(wdirs_yaw,
                                                                wec_power_matrix) 
                                                        
            self.wec_original_power_matrix = self.wec_power_matrix.copy()
            
            return -1
        
        Nm = np.array(range(-self.order.max(), self.order.max()+1), dtype=float)
        Gp = np.transpose(self.force_matrix, axes = (0,2,1))

        ind = 0
        max_iter = len(wave_angles)*len(te)*len(hm0)
        my_bounds = [(-1,None),]*ndof*2
        x0 = [np.zeros((ndof*2),"f")]  # x2 because two parameter are assesed
        for h_ind in range(len(hm0)):
            for t_ind in range(len(te)):
                for d_ind in range(len(wave_angles)):
                    if self.debug: ind += 1; print('sea state {} over {}'.format(ind, max_iter))
                    for orient, thetas in zip(wdirs_yaw[::2], wdirs_yaw[1::2]):
                        if not wave_angles[d_ind] in thetas: continue
                        ths = np.array([th for th in thetas] , dtype=float)
                        g_rot = np.transpose(Gp*np.exp(-1j*Nm*orient), axes = (0,2,1))
                        # Excitation Force isolated WEC
                        (mths, mode) = np.meshgrid(ths, Nm, indexing = 'ij')
                        AP = np.exp(-1j*mode*(np.pi/2+mths))
                        f_ex = np.zeros((len(periods),len(thetas),ndof), dtype = complex)
                        for f in range(len2(self.periods)):
                            f_ex[f,:,:]= np.dot(AP,g_rot[f])
                        ref_power = machine_spec['power_matrix'][t_ind, h_ind, d_ind]
                        lin_fit_b = optimize.fmin_l_bfgs_b(self.__power_fit, x0,
                                                           args=((t_ind, h_ind, d_ind), thetas, ref_power, f_ex),
                                                           approx_grad=1,
                                                           bounds=my_bounds,
                                                           pgtol=1e-12,
                                                           iprint=int(self.debug))
                        Cfit_sm[t_ind, h_ind, d_ind] = np.diag(lin_fit_b[0][:ndof])
                        Kfit_sm[t_ind, h_ind, d_ind] = np.diag(lin_fit_b[0][ndof:])
                        wec_power_matrix[t_ind, h_ind, d_ind] = self.__power_fit(lin_fit_b[0], 
                                                                                (t_ind, h_ind, d_ind), 
                                                                                thetas, 0, f_ex)
        self.c_fit = Cfit_sm
        self.k_fit = Kfit_sm
        self.wec_power_matrix = np.sqrt(wec_power_matrix)
        wec_original_power_matrix = np.zeros(wec_power_matrix.shape)
        self.wec_original_power_matrix = self.__get_wec_power_matrix(wdirs_yaw,
                                                                     wec_original_power_matrix,
                                                                     skip_fit=True)

    def __power_fit(self, x, seastate_id, thetas, target_pow, f_ex):
        """
        power_fit: calls the power_wec to assess the power production of the WEC, given a frequency domain model and the sea states

        Args:
            x (list): flat list of Cfit and Kfit coefficients
            dir_id (int): index of the current wave direction
            Sp_dens (numpy.ndarray): spectral density of the current sea states
            PMpow (float): user given power associated with the current sea state
            Cpto (numpy.ndarray): pto damping matrix associated with the current sea state
            Kmooring (numpy.ndarray): mooring stiffness matrix associated with the current sea state

        Returns:
            (float): square error between user given and calculated power associated with the current sea state
        """
        Nx = len(x)/2
        
        ScatDiag_mod = [np.ones((1,1,1)), self.scatter_diagram_spec]
        Cpto_mod = self.c_pto[seastate_id[0], seastate_id[1], seastate_id[2]].reshape((1, 1, 1, Nx, Nx))
        Kmoor_mod = self.k_mooring[seastate_id[0], seastate_id[1], seastate_id[2]].reshape((1, 1, 1, Nx, Nx))
        
        Kext_mod = self.k_ext[seastate_id[0], seastate_id[1], seastate_id[2]].reshape((1, 1, 1, Nx, Nx))
        Cext_mod = self.c_ext[seastate_id[0], seastate_id[1], seastate_id[2]].reshape((1, 1, 1, Nx, Nx))
        
        # normalised by Khst
        k_fit = x[Nx:]
        k_fit = np.diag(k_fit).reshape(Cpto_mod.shape) + Kext_mod
        
        # normalised by Cpto
        c_fit = x[:Nx]
        c_fit = np.diag(c_fit).reshape(Cpto_mod.shape) + Cext_mod
        
        _, abs_power = EnergyProduction(1,
                                        [self.wave_dir[seastate_id[2]]],
                                        [self.hm0[seastate_id[1]]],
                                        [self.tp[seastate_id[0]]],
                                         thetas,
                                         self.periods,
                                         ScatDiag_mod,
                                         self.m_m,
                                         self.m_add,
                                         Cpto_mod,
                                         self.c_rad,
                                         Kmoor_mod,
                                         self.k_hst,
                                         f_ex,
                                         k_fit,
                                         c_fit)
       
        return ((target_pow - abs_power[0,0,0,0]))**2

    def __get_wec_power_matrix(self, wdirs_yaw, wec_power_matrix, skip_fit=False, no_yaw=False):
        if no_yaw:
            x_scale = 0.0
            for h_ind in range(len(self.hm0)):
                for t_ind in range(len(self.te)):
                    for d_ind in range(len(self.wave_dir)):
                        ind = np.argmin(np.abs(self.wave_dir[d_ind]-self.directions))
                        f_ex = self._fex[:,ind,:].reshape((len(self.periods), 1, self.n_dof))
                        x = np.r_[np.diag(self.c_fit[t_ind, h_ind, d_ind]),np.diag(self.k_fit[t_ind, h_ind, d_ind])]
                        wec_power_matrix[t_ind, h_ind, d_ind] = self.__power_fit(x*x_scale,(t_ind, h_ind, d_ind), [self.wave_dir[d_ind]], 0., f_ex)
            
            return np.sqrt(wec_power_matrix)
        else:
            Nm = np.array(range(-self.order.max(), self.order.max()+1), dtype=float)
            Gp = np.transpose(self.force_matrix, axes = (0,2,1))
            
            x_scale = 1.0
            if skip_fit: x_scale = 0.0;
            for h_ind in range(len(self.hm0)):
                for t_ind in range(len(self.te)):
                    for d_ind in range(len(self.wave_dir)):
                        for orient, thetas in zip(wdirs_yaw[::2], wdirs_yaw[1::2]):
                            if not self.wave_dir[d_ind] in thetas: continue
                            ths = np.array([th for th in thetas] , dtype=float)
                            g_rot = np.transpose(Gp*np.exp(-1j*Nm*orient), axes = (0,2,1))
                            # Excitation Force isolated WEC
                            (mths, mode) = np.meshgrid(ths, Nm, indexing = 'ij')
                            AP = np.exp(-1j*mode*(np.pi/2+mths))
                            f_ex = np.zeros((len(self.periods),len(thetas),self.n_dof), dtype = complex)
                            for f in range(len2(self.periods)):
                                f_ex[f,:,:]= np.dot(AP,g_rot[f])
                            
                            x = np.r_[np.diag(self.c_fit[t_ind, h_ind, d_ind]),np.diag(self.k_fit[t_ind, h_ind, d_ind])]
                            wec_power_matrix[t_ind, h_ind, d_ind] = self.__power_fit(x*x_scale, 
                                                                                    (t_ind, h_ind, d_ind), 
                                                                                    thetas, 0., f_ex)
    
            return np.sqrt(wec_power_matrix)
        
if __name__ == "__main__":


    import pickle
    pkl_file = open('tests/test_data.pkl', 'r')

    site_spec = {'spec_shape':'Jonswap'}
    site_spec['spec_gamma'] = 1.0
    site_spec['spec_spreading'] = -1
    site_spec['te'] = np.linspace(3,10,5)
    site_spec['hm0'] = np.linspace(0.5, 3.5, 3)
    site_spec['wave_angles'] = np.linspace(0,270,1)
    site_spec['probability_of_occurence'] = np.ones((5,3,1))

    machine_spec = {'c_pto': np.ones((5,3,1,7,7))*1.0e6}
    machine_spec['k_mooring'] = np.zeros((5,3,1,7,7))
    machine_spec['power_matrix'] = np.zeros((5,3,1))
    machine_spec['yaw'] = 0.0

    (per, dirs, m_m, m_add, c_rad, k_hst, force_tr_mat,order,pto_dof,mooring_dof,fex) = pickle.load(pkl_file)

    per_fit = PowerMatrixFit(per, dirs, m_m, m_add, c_rad, np.ones((5,3,1,7,7)), k_hst, np.ones((5,3,1,7,7)), force_tr_mat, pto_dof-1, np.array([0]), order,"pippo")
    per_fit.perf_fitting(machine_spec, site_spec )
