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

import dtocean_wave.utils.hdf5_interface as h5i


class DbReader():
    def __init__(self, case, db_folder, debug=False):
        self.wec_ID = case
        self.db_folder = db_folder
        
    def load_data(self):
        f_n = os.path.join(self.db_folder, 'case{}'.format(self.wec_ID), 'case{}_hyd.h5'.format(self.wec_ID)) 
        dic = h5i.load_dict_from_hdf5(f_n)
        self.m_m = dic['m_m']
        self.m_add = dic['m_add']
        self.c_rad = dic['c_rad']
        self.f_ex = dic['f_ex']
        self.periods = dic['periods']
        self.directions = dic['directions']
        self.k_hst = dic['k_hst']
        self.diffraction_tr_mat = dic['diffraction_tr_mat']
        self.force_tr_mat = dic['force_tr_mat']
        self.amplitude_coefficient_radiation = dic['amplitude_coefficient_radiation']
        self.water_depth = dic['water_depth']
        self.cyl_radius = dic['cyl_radius']
        self.modes = dic['modes']
        self.order = dic['max_order']
        self.truncation_order = dic['truncation_order']
        self.pto_dof = dic['pto_dof']
        self.moor_dof = dic['mooring_dof']               
        
if __name__ == "__main__":
    reader = DbReader(1, db_folder=r"..\src\wec_db")
    reader.load_data()