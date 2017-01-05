#!/usr/bin/python2.7
# encoding: utf-8
"""
This module contains the methods used to load the hydrodynamic module of the isolated WEC
using the Nemoh software.

.. module:: read_bem_solution
   :platform: Windows
   :synopsis: WEC Numerical model loader

.. moduleauthor:: Francesco Ferri <ff@civil.aau.dk>
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
    tp = dic['tp']
    hm0 = dic['hm0']
    wave_dir = dic['wave_dir']
    scatter_diagram = dic['scatter_diagram']
    power_matrix = dic['power_matrix']

    return (c_fit+c_ext, k_fit+k_ext, c_pto, k_mooring, tp, hm0, wave_dir, scatter_diagram, power_matrix)


def __check_filename(data_folder, fn="wec_solution.h5"):
    f_n = os.path.join(data_folder, 'wec_solution.h5')
    if not os.path.isfile(f_n):
        raise ValueError("Missing the wec_solution.h5 file in the specified data folder [{}] ".format(data_folder))

    return f_n

