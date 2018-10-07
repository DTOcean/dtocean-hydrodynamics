# -*- coding: utf-8 -*-

#    Copyright (C) 2016 Francesco Ferri
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
Created on Mon Apr 25 15:16:44 2016

.. moduleauthor:: Francesco Ferri <ff@civil.aau.dk>
"""

import pandas as pd
import numpy as np


def save_dict_to_pandas_and_hdf5(dic, filename):
    """
    ....
    """
    h5file = pd.HDFStore(filename)
    modes = pd.DataFrame(dic['modes'], columns= ['type','ax_x','ax_y','ax_z','p_x','p_y','p_z'])
    modes.index.name = 'dof_ID'

    General = pd.Series([dic['water_depth'], dic['cyl_radius'], 
                         modes, dic['order'], 
                         dic['truncation_order']], 
                         ['water_depth','cylinder_radius','dofs','order','truncation_order'])

    Periods = pd.Series(dic['periods'])
    Dirs = pd.Series(dic['directions'])

    M_m = pd.DataFrame(dic['m_m'])
    M_m.columns.name = ['dof_j']
    M_m.index.name = ['dof_i']

    K_hst = pd.DataFrame(dic['k_hst'])
    K_hst.columns.name = ['dof_j']
    K_hst.index.name = ['dof_i']

    M_add = send_to_pandas(dic['m_add'], 'dof_i', 'dof_j')
    C_rad = send_to_pandas(dic['c_rad'], 'dof_i', 'dof_j')
    F_ex = send_to_pandas(dic['f_ex'], 'dof_i', 'direction')
    Diffraction_matrix = send_to_pandas(dic['diffraction_tr_mat'], 'mode_i', 'mode_j')
    Force_matrix = send_to_pandas(dic['force_tr_mat'], 'mode_i', 'dof_i')
    Amplitude_coeff = send_to_pandas(dic['amplitude_coefficient_radiation'], 'dof_i', 'mode_i')
           
    C_ext = send_to_pandas(dic['c_ext'], 'te', 'hm0')
    K_ext = send_to_pandas(dic['k_ext'], 'te', 'hm0')
    C_fit = send_to_pandas(dic['c_fit'], 'te', 'hm0')
    K_fit = send_to_pandas(dic['k_fit'], 'te', 'hm0')
    C_pto = send_to_pandas(dic['c_pto'], 'te', 'hm0')
    K_moor = send_to_pandas(dic['k_mooring'], 'te', 'hm0')
    
    Scatter_diagram = send_to_pandas(dic['scatter_diagram'], 'te', 'hm0')
    Te = send_to_pandas(dic['te'])
    Hm0 = send_to_pandas(dic['hm0'])
    Wave_dir = send_to_pandas(dic['wave_dir'])
    

    h5file.put('Mass Matrix', M_m)
    h5file.put('Added Mass', M_add)
    h5file.put('Periods', Periods)
    h5file.put('Directions', Dirs)
    h5file.put('Hydrostatic', K_hst)
    h5file.put('Radiation Damping', C_rad)
    h5file.put('Excitation', F_ex)
    h5file.put('General', General)
    h5file.put('Diffraction_matrix', Diffraction_matrix)
    h5file.put('Force_matrix', Force_matrix)
    h5file.put('Amplitude_coeff', Amplitude_coeff)
    

    h5file.put('External damping', C_ext)
    h5file.put('External stifness', K_ext)
    h5file.put('Fit damping', C_fit)
    h5file.put('Fit stifness', K_fit)
    h5file.put('PTO damping', C_pto)
    h5file.put('Mooring stifness', K_moor)
    h5file.put('scatter_diagram', Scatter_diagram)
    h5file.put('te', Te)
    h5file.put('hm0', Hm0)
    h5file.put('wave_dir', Wave_dir)
    
    h5file.close()
    

def send_to_pandas(array, ind_name='ax1', col_name='ax2'):
    if len(array.shape) == 1:
        Array = pd.Series(array)
    elif len(array.shape) == 2:
        Array = pd.DataFrame(array)
    elif len(array.shape) == 3:
        Array = pd.Panel(array)
    elif len(array.shape) == 4:
        array_ls = []
        array_mt = []
        ax1 = array.shape[0]
        ax2 = array.shape[1]
        for ind in range(ax1):
            for ind2 in range(ax2):
                temp = pd.DataFrame(array[ind, ind2, :, :])
                temp.columns.name = ind_name
                temp.index.name = col_name
                array_ls.append(temp)
            array_mt.append(array_ls)
            array_ls = []
        Array = pd.DataFrame(array_mt)
    else:
        raise IOError("Invalid array shape")
    
    return Array

def load_dict_from_pandas_and_hdf5(filename):
    
    outputs = {}
    outputs['amplitude_coefficient_radiation'] = np.asarray(pd.read_hdf(filename,'Amplitude_coeff'))
    outputs['k_hst'] = np.asarray(pd.read_hdf(filename,'Hydrostatic'))
    outputs['m_add'] = np.asarray(pd.read_hdf(filename,'Added Mass'))
    outputs['m_m'] = np.asarray(pd.read_hdf(filename,'Mass Matrix'))
    outputs['periods'] = np.asarray(pd.read_hdf(filename,'Periods'))
    outputs['directions'] = np.asarray(pd.read_hdf(filename,'Directions'))
    outputs['k_hst'] = np.asarray(pd.read_hdf(filename,'Hydrostatic'))
    outputs['c_rad'] = np.asarray(pd.read_hdf(filename,'Radiation Damping'))
    general = pd.read_hdf(filename,'General')
    outputs['water_depth'] = general.water_depth
    outputs['cyl_radius'] = general.cylinder_radius
    outputs['modes'] = general.dofs
    outputs['truncation_order'] = general.truncation_order
    outputs['order'] = general.order
    outputs['diffraction_tr_matrix'] = np.asarray(pd.read_hdf(filename,'Diffraction_matrix'))
    outputs['force_tr_matrix'] = np.asarray(pd.read_hdf(filename,'Force_matrix'))
    outputs['f_ex'] = np.asarray(pd.read_hdf(filename,'Excitation'))
    outputs['c_ext'] = np.asarray(pd.read_hdf(filename,'External damping'))
    outputs['k_ext'] = np.asarray(pd.read_hdf(filename,'External stifness'))
    outputs['c_fit'] = np.asarray(pd.read_hdf(filename,'Fit damping'))
    outputs['k_fit'] = np.asarray(pd.read_hdf(filename,'Fit stifness'))
    outputs['c_pto'] = np.asarray(pd.read_hdf(filename,'PTO damping'))
    outputs['k_mooring'] = np.asarray(pd.read_hdf(filename,'Mooring stifness'))
    outputs['te'] = np.asarray(pd.read_hdf(filename, 'te'))
    outputs['hm0'] = np.asarray(pd.read_hdf(filename, 'hm0'))
    outputs['wave_dir'] = np.asarray(pd.read_hdf(filename, 'wave_dir'))
    outputs['scatter_diagram'] = np.asarray(pd.read_hdf(filename, 'scatter_diagram'))
     
    return outputs
    

if __name__ == "__main__":
    outputs = load_dict_from_pandas_and_hdf5(r'C:\Users\francesco\Desktop\test_prj1\test2.h5')
    