# -*- coding: utf-8 -*-
"""
Created on Fri May 27 16:41:27 2016

@author: francesco
"""
import utils.hdf5_interface as h5i 
import os

path = r'C:\Users\francesco\Desktop'
filename = 'test_rm3_sterling_data_collection.hdf5'

dic = h5i.load_dict_from_hdf5(os.path.join(path, filename))


# h5i.save_dict_to_hdf5(dic, os.path.join(path, filename))
