# -*- coding: utf-8 -*-

#    Copyright (C) 2016 Francesco Ferri
#    Copyright (C) 2019 Mathew Topper
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
This module contains the methods used to write and read dictionaries to hdf5
files

.. module:: hdf5_interface
   :platform: Windows
   :synopsis: hdf5 reader and writer

.. moduleauthor:: Francesco Ferri <ff@civil.aau.dk>
.. moduleauthor:: Mathew Topper <mathew.topper@dataonlygreater.com>
"""

import os

import h5py
import numpy as np


def create_empty_project(folder, prj_name):

    general_inputs = {'angle_def':[],
                      'frequency_def':[],
                      'mooring_dof':[],
                      'pto_dof':[],
                      'ndof':[],
                      'water_depth':0,
                      'get_array_mat':0,
                      'cyl_nzeta':[],
                      'cyl_ntheta':[],
                      'data_folder': '',
                      'input_type': 1}
    
    body = -1
            
    bodies = body
    body_inputs = {'shared_dof':[0,0,0,0,0,0],
                    'local_cs': [],
                    'body': bodies}
    filename = os.path.join(folder,
                            ''.join([prj_name, '_data_collection.hdf5']))

    dic = {'inputs_hydrodynamic': {'general_input':general_inputs,
                                   'body_inputs': body_inputs},
           'prj_filename': filename,
           'prj_folder': folder,
           'prj_name': prj_name}
    
    save_dict_to_hdf5(dic, filename)
     
    return load_dict_from_hdf5(filename)


def save_dict_to_hdf5(dic, filename):
    
    with h5py.File(filename, 'w') as h5file:
        recursively_save_dict_contents_to_group(h5file, '/', dic)


def recursively_save_dict_contents_to_group(h5file, path, dic):

    for key, item in dic.items():
        if isinstance(item, (np.ndarray,
                             np.int64,
                             np.float64,
                             float,
                             int,
                             list,
                             str,
                             bytes)):
            h5file[path + key] = item
        elif isinstance(item, dict):
            recursively_save_dict_contents_to_group(h5file,
                                                    path + key + '/', item)
        else:
            raise ValueError('Cannot save {} type'.format(item))


def load_key_from_hdf5(filename, key_str):
    
    with h5py.File(filename, 'r') as h5file:
        data = np.asarray(h5file.get(key_str).value)

    return data


def load_dict_from_hdf5(filename):

    with h5py.File(filename, 'r') as h5file:
        return recursively_load_dict_contents_from_group(h5file, '/')


def recursively_load_dict_contents_from_group(h5file, path):

    ans = {}
    
    for key, item in h5file[path].items():
        if isinstance(item, h5py._hl.dataset.Dataset):
            ans[key] = item.value
        elif isinstance(item, h5py._hl.group.Group):
            ans[key] = \
                    recursively_load_dict_contents_from_group(h5file,
                                                              path + key + '/')
    
    return ans


if __name__ == '__main__':

    data = {'x': 'astring',
            'y': np.arange(10),
            'd': {'z': np.ones((2,3)),
                  'b': b'bytestring'}}
    print(data)
    filename = 'test.h5'
    save_dict_to_hdf5(data, filename)
    dd = load_dict_from_hdf5(filename)
    print(dd)
