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
This module contains the methods used to load the hydrodynamic module of the isolated WEC
using the Nemoh software.

.. module:: read_bem_solution
   :platform: Windows
   :synopsis: WEC Numerical model loader

.. moduleauthor:: Francesco Ferri <ff@civil.aau.dk>
.. moduleauthor:: Mathew Topper <mathew.topper@dataonlygreater.com>
"""

import os
import hdf5_interface as h5i


def read_hydrodynamic_solution(data_folder):
    dic = h5i.load_dict_from_hdf5(__check_filename(data_folder))
    m_m = dic['m_m']
    m_add = dic['m_add']
    c_rad = dic['c_rad']
    f_ex = dic['f_ex']
    periods = dic['periods']
    directions = dic['directions']
    k_hst = dic['k_hst']
    diffraction_tr_mat = dic['diffraction_tr_mat']
    force_tr_mat = dic['force_tr_mat']
    amplitude_coefficient_radiation = dic['amplitude_coefficient_radiation']
    water_depth = dic['water_depth']
    cyl_radius = dic['cyl_radius']
    modes = dic['modes']
    order = dic['max_order']
    truncation_order = dic['truncation_order']

    return (m_m, m_add, c_rad, k_hst, f_ex, periods, directions,
            diffraction_tr_mat, force_tr_mat, amplitude_coefficient_radiation,
            water_depth, cyl_radius, modes, order, truncation_order)


def read_freq(data_folder):
    f_n = __check_filename(data_folder)
    data = h5i.load_key_from_hdf5(f_n, 'periods')
    return 1./data


def read_directions(data_folder):
    f_n = __check_filename(data_folder)
    return h5i.load_key_from_hdf5(f_n, 'directions')


def read_performancefit_solution(data_folder):

    dic = h5i.load_dict_from_hdf5(__check_filename(data_folder))
    c_ext = dic['c_ext']
    k_ext = dic['k_ext']
    c_fit = dic['c_fit']
    k_fit = dic['k_fit']
    c_pto = dic['c_pto']
    k_mooring = dic['k_mooring']
    te = dic['te']
    hm0 = dic['hm0']
    wave_dir = dic['wave_dir']
    scatter_diagram = dic['scatter_diagram']
    power_matrix = dic['power_matrix']

    return (c_fit+c_ext,
            k_fit+k_ext,
            c_pto,
            k_mooring,
            te,
            hm0,
            wave_dir,
            scatter_diagram,
            power_matrix)


def __check_filename(data_folder, fn="wec_solution.h5"):
    f_n = os.path.join(data_folder, 'wec_solution.h5')
    if not os.path.isfile(f_n):
        raise ValueError("Missing the wec_solution.h5 file in the specified data folder [{}] ".format(data_folder))

    return f_n

