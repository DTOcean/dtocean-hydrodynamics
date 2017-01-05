# -*- coding: utf-8 -*-
"""
Created on Wed Apr 20 18:19:28 2016

@author: frafe
"""

from pandas import read_hdf
import os
# this query selects the columns A and B
# where the values of A is greather than 0.5
output_folder = r'C:\Users\francesco\Desktop\test_prj1'
f_n = os.path.join(output_folder,'wec_object.h5')
m_add = read_hdf(f_n,'Added Mass')
m_m = read_hdf(f_n,'Mass Matrix')
periods = read_hdf(f_n,'Periods')
directions = read_hdf(f_n,'Directions')
k_hst = read_hdf(f_n,'Hydrostatic')
c_rad = read_hdf(f_n,'Radiation Damping')
general = read_hdf(f_n,'General')
diffraction_matrix = read_hdf(f_n,'Diffraction_matrix')
force_matrix = read_hdf(f_n,'Force_matrix')
f_ex = read_hdf(f_n,'Excitation')
c_ext = read_hdf(f_n,'External damping')
k_ext = read_hdf(f_n,'External stifness')
c_fit = read_hdf(f_n,'Fit damping')
k_fit = read_hdf(f_n,'Fit stifness')
c_pto = read_hdf(f_n,'PTO damping')
k_moor = read_hdf(f_n,'Mooring stifness')