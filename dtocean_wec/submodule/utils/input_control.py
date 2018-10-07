# -*- coding: utf-8 -*-

#    Copyright (C) 2016 Francesco Ferri, Pau Mercadez Ruiz
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

import os

import file_utilities as f_util


def check_hydrodynamic_files(folder, get_array_mat=True):
    # using a sequential apporach this check comes after the check_folder.
    # therefore there is no need to re-check the directory existence
    base_file_list = ('nemoh.cal',)
    for el in base_file_list:
        if not os.path.isfile(os.path.join(folder, el)):
            return (False, 'Missing the {} file in the results folder'.format(el))

    hdy_folder = os.path.join(folder)
    base_file_list = ('ProblemDescription.txt', 'CA.dat', 'CM.dat', 'ExcitationForce.tec')
    for el in base_file_list:
        if not os.path.isfile(os.path.join(hdy_folder, 'results', el)):
            return (False, 'Missing the {} file in the results folder'.format(el))
    
    with open(os.path.join(hdy_folder,  'results', 'ProblemDescription.txt')) as fid:
        n_problems = int(float(fid.readline().split()[-1]))
    
    if get_array_mat:
        for problem in range(n_problems):
            cyl_pot = 'cylsurface.{:5d}.dat'.format(problem+1)
            if not os.path.isfile(os.path.join(hdy_folder,  'results', cyl_pot)):
                return (False, 'Missing the {} file in the results folder'.format(el))
            
    return (True, '')
        

def check_hydrostatic_files(folder, n_bodies):
    # using a sequential apporach this check comes after the check_folder.
    # therefore there is no need to re-check the directory existence
    hst_folder = os.path.join(folder, 'hydrostatic')
    for body in range(n_bodies):
        if not os.path.isfile(os.path.join(hst_folder,
                                           'body{}'.format(body), 'mesh', 'KH.dat')):
            return (False, 'Missing KH.dat file for body {}'.format(body))
            
    return (True, '')

def check_folder_tree(folder):
    # level 0
#    level0 = ('hydrodynamic', )
#    sub_l0 = next(os.walk(folder))[1]
#    for el in level0:
#        if not el in sub_l0:
#            return (False, 'Missing {} folder'.format(el))
    
    # level1
#    level1_st = tuple(['body{}'.format(el) for el in range(n_bodies)])
#    sub_l1_st = next(os.walk(os.path.join(folder, 'hydrostatic')))[1]
#    for el in sub_l1_st:
#        if not el in level1_st:
#            return (False, 'Missing {} folder'.format(el))
      
    level1_dy = ('results','mesh')
    sub_l1_dy = next(os.walk(os.path.join(folder)))[1]
    for el in sub_l1_dy:
        if not el in level1_dy:
            return (False, 'Missing {} folder'.format(el))
    
    # level 2
#    for body in range(n_bodies):
#        sub_l2_st = next(os.walk(os.path.join(folder, 'hydrostatic', 'body{}'.format(body))))[1]
#        if not 'mesh' in sub_l2_st:
#            return (False, 'Missing mesh folder for body {}'.format(body))
            
    return (True, '')

def check_nemoh_results(data_folder, get_array_mat=True):
    #status_files = check_nemoh_inputs(data_folder)
    #if not status_files[0]:
    #    raise ValueError(status_files[1])
#    # hydrostatic: search for KH.dat in each body folder
#    f_mds = f_util.load_file(data_folder, "WEC_mds.wp2")
#    ind = 1
#    n_dof_sh = f_util.split_string(f_mds[ind], num_type=int)[0]
#    ind += n_dof_sh+2
#    n_bodies = f_util.split_string(f_mds[ind], num_type=int)[0]
    
    # check base folder structure
    # hydrostatic, hydrostatic/body#, hydrostatic/body#/mesh
    # hydrodynamic, hydrodynamic/results, hydrodynamic/mesh
    if not os.path.isdir(data_folder):
        return (False, 'The specified data folder is invalid')
    status_tree = check_folder_tree(data_folder)
    if not status_tree[0]:
        return status_tree
        
    # check the required files
    # status_st_files = check_hydrostatic_files(data_folder, n_bodies)
    #if not status_st_files[0]:
    #    return status_st_files
    status_dy_files = check_hydrodynamic_files(data_folder, get_array_mat=get_array_mat)
    if not status_dy_files[0]:
        return status_dy_files

    return (True, '')

def check_nemoh_inputs(data_folder):
    """check_nemoh_inputs searches for the three input files 'WEC_gnr.wp2', 'WEC_mmx.wp2', 'WEC_mds.wp2'
    and relative mesh file(s), required by the software to generate the solution of the bem problem. 

    """
    input_ls = ['WEC_gnr.wp2', 'WEC_mmx.wp2', 'WEC_mds.wp2']
    for el in input_ls:
        if not os.path.isfile(os.path.join(data_folder, el)):
            return (False, 'Missing {}'.format(el))
   
    # f_gnr = f_util.load_file(data_folder, "WEC_gnr.wp2")
    # f_mmx = f_util.load_file(data_folder, "WEC_mmx.wp2")
   
    # check mesh files 
    f_mds = f_util.load_file(data_folder, "WEC_mds.wp2")
    ind = 3
    n_bodies = f_util.split_string(f_mds[ind], num_type=int)[0]
    ind += 3
    mesh_files = []
    for body in range(n_bodies):
        mesh_files += [f_mds[ind].split()[0]]
        ind += 4
        ndof_p = f_util.split_string(f_mds[ind], num_type=int)[0]
        ind += ndof_p+1
        nchild = f_util.split_string(f_mds[ind], num_type=int)[0]
        ind += nchild+3
    
    for el in mesh_files:
        if not os.path.isfile(os.path.join(data_folder, el)):
            return (False, 'Missing {}'.format(el))
    
    return (True, '')

def check_wamit_results(data_folder, get_array_mat=True):
    gnr_input_ls = []
    wamit_req_ls = [".cfg", ".pot", ".mmx", ".1", ".2", ".hst", ".fpt"]
    for el in gnr_input_ls:
        if not os.path.isfile(os.path.join(data_folder, el)):
            return (False, 'Missing {}'.format(el))
    try:
        map(check_filetype_in_folder, wamit_req_ls, [data_folder,]*len(wamit_req_ls))
    except Exception as e:
        return (False, "WAMIT requirements not met. {}".format(e))
    
    check_file_type = [el for el in os.listdir(data_folder) if el.endswith('.6p')]
    if get_array_mat and not check_file_type:
        return (False, "The specified filename (.6) does not exist in the folder {}".format( data_folder))
    elif len(check_file_type) > 1:  # this cndition cannot exist, just keep it here in case the comparison is bugged
        return (False, "The specified filename (.6) exists in multiple copies in the folder {}. ".format(data_folder),
                         "Remove the unused files")
        
    return (True, )

def check_filetype_in_folder(ftype, data_folder):
    check_file_type = [el for el in os.listdir(data_folder) if el.endswith(ftype)]
    if not check_file_type:
        raise NameError("The specified filename ({}) does not exist in the folder {}".format(ftype, data_folder))
    elif len(check_file_type) > 1:  # this cndition cannot exist, just keep it here in case the comparison is bugged
        raise ValueError("The specified filename ({}) exists in multiple copies in the folder {}. ".format(ftype, data_folder),
                         "Remove the unused files")


    
if __name__ == "__main__":
    print(check_nemoh_inputs(r"C:\Users\francesco\Desktop\Nemoh_project"))
    print(check_nemoh_results(r"C:\Users\francesco\Desktop\Nemoh_project"))
    