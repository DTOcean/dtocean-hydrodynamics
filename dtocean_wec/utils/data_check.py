# -*- coding: utf-8 -*-
"""
Created on Fri May 27 19:07:01 2016

@author: francesco
"""
import numpy as np
import os

def check_inputs_hydrodynamic(hyd_d, prj_path):
    g_i = hyd_d['general_input']
    ndof = g_i['ndof']
    pto_dof = g_i['pto_dof']
    mooring_dof = g_i['mooring_dof']
    freq_def = g_i['frequency_def']
    ang_def = g_i['angle_def']
    
    i_t = g_i['input_type']
    err = []
    if not isinstance(ndof, (list, np.ndarray)):
        err.append('dof: type error')
    elif not len(ndof) == 1:
        err.append('dof: length error')
    if not isinstance(pto_dof, (int, float, list, np.ndarray)):
        err.append('pto_dof: type error')
    if not isinstance(mooring_dof, (int, float, list, np.ndarray)):
        err.append('mooring_dof: type error')
    if i_t == 2:
        if not isinstance(freq_def, (list, np.ndarray)):
            err.append('freq_def: type error')
        elif not len(freq_def) == 3:
            err.append('freq_def: length error')
        elif not freq_def[0] > 0:
            err.append("the number of frequency need to be a positive integer")
        if not isinstance(ang_def, (list, np.ndarray)):
           err.append('ang_def: type error')
        elif not len(ang_def) == 1:
            err.append('ang_def: length error')
        elif not ang_def[0] > 0:
            err.append("the number of angles need to be a positive integer")
    elif i_t ==3:
        if not isinstance(freq_def, (list, np.ndarray)):
            err.append('freq_def: type error')
        if not isinstance(ang_def, (list, np.ndarray)):
            err.append('ang_def: type error')
       
    if not np.all(pto_dof <= ndof):
        err.append("pto dof greater then the total number of dof")
    if not np.all(mooring_dof <= ndof):
        err.append("mooring dof greater then the total number of dof")
    
    if g_i['water_depth'] == 0 and g_i['get_array_mat'] == 1:
        err.append("the water depth cannot be zero when the array interaction matrixes need to be solved.")
        
    
    b_i = hyd_d['body_inputs']
    if not isinstance(b_i['local_cs'], (list, np.ndarray)):
        err.append("local_cs: type error")
    elif not len(b_i['local_cs']) == 3:
        err.append("local_cs: length error")
        
    if not isinstance(b_i['shared_dof'], (list, np.ndarray)):
        err.append("shared_dof: type error")
    elif not len(b_i['shared_dof']) == 6:
        err.append("shared_dof: lenght error")
    elif not any(np.asarray(b_i['shared_dof']) > 0):
        err.append("No degrees of freedom has been selected")
    elif np.any([True for el in b_i['shared_dof'] if el>1]) or np.any([True for el in b_i['shared_dof'] if el<0]):
        err.append("shared_dof: value error")
    
    if 'body' in b_i.keys():
        bodies = b_i['body'].values()
        if bodies:
            for bd in bodies:
                err += check_mass(bd['mass'])
                err +=check_inertia(bd['inertia'])
                err += check_ID(bd['ID'])
                err += check_mesh(bd['mesh'], os.path.join(prj_path, "raw_data"))
                err += check_cog(bd['cog'])
                err += check_child(bd['child_dof_position'])
                err += check_parent(bd['parent_body'])
                err += check_dof(bd['dof_with_parent'])
                err += check_cs(bd['axis_angles'])
        else:
            err.append("body tab: no body specified")
    else:
            err.append("body tab: no body specified")
    return err
        
def check_mass(el):
    if not isinstance(el, (list, np.ndarray)):
        return ['mass: type error']
    elif not len(el) == 1:
        return ['mass: length error']
    
    return []
    
def check_inertia(el):
    if not isinstance(el, (list, np.ndarray)):
        return ['inertia: type error']
    elif not len(el) == 3:
        return 'inertia: length error'
    elif len(el[0]) != 3 or len(el[1]) != 3 or len(el[2]) != 3:
        return ['inertia: length error']
     
    return []
    
def check_ID(el):
    if not isinstance(el, (list, np.ndarray)):
        return ['ID: type error']
    elif not len(el) == 1:
        return ['ID: length error']
    return []
    
def check_mesh(el, path):
    if not isinstance(el, (str, unicode)):
		print('mesh error')
		print(el)
		print(type(el))
		return ['mesh: type error']
    elif not os.path.isfile(os.path.join(path,el)):
		return ['mesh: value error']
    return []
    
def check_cog(el):
    if not isinstance(el, (list, np.ndarray)):
        return ['cog: type error']
    elif not len(el) == 3:
        return ['cog: length error']  
    return []
    
def check_child(el):
    if not isinstance(el, (list, np.ndarray, int)):
        return ['child: type error']
    if el == -1:
        return []
    else:
        nch = len(el)
        for eli in range(nch):
            if not isinstance(el[eli], (list, np.ndarray)):
                return ['child: type error']
            elif not el[eli]:
                pass
            elif not len(el[eli]) == 4:
                return ['child: length error']   
    return []
    
def check_parent(el):
    if not isinstance(el, (list, np.ndarray)):
        return ['parent: type error']
    elif not len(el) == 1:
        return ['parent: length error']
    return []
    
def check_dof(el):
    if not isinstance(el, (list, np.ndarray, int)):
        return ['dof: type error']
    if el == -1:
        return []
    else:
        ndof = len(el)
        for eli in range(ndof):
            if not isinstance(el[eli], list):
                return ['dof: type error']
            elif not el[eli]:
                pass
            elif not len(el[eli]) == 4:
                return ['dof: length error']
            elif not el[eli][0] == 1 and not el[eli][0] == 2:
                return ['dof: format error']
    return [] 
    
def check_cs(el):
    if not isinstance(el, (list, np.ndarray)):
        return ['cs: type error']
    elif not len(el) == 3:
        return ['cs: length error']
    
    return []

    