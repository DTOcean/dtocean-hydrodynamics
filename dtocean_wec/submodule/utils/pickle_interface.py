# -*- coding: utf-8 -*-
"""
Created on Mon Apr 25 15:45:08 2016

@author: francesco
"""
import pickle

def save_dict_to_pickle(dic, filename):
    with open(filename,'wb') as f_hyd:
        pickle.dump(dic, f_hyd)
        
        
def load_dict_from_pickle(filename):
    return pickle.load(open(filename,'rb'))
