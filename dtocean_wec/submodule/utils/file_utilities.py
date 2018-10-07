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
import sys
import shutil
from distutils.util import strtobool


def query_yes_no(question):
    # query the user on the terminal for a y or n answer
    sys.stdout.write('%s [y/n]\n' % question)
    while True:
        try:
            return strtobool(raw_input().lower())
        except ValueError:
            sys.stdout.write('Please respond with \'y\' or \'n\'.\n')

def copy_result_to_project(source, destination, file_fl=False):
    src_st = os.path.join(source,'hydrostatic')
    src_dy = os.path.join(source,'hydrodynamic')
    dst_st = os.path.join(destination,'hydrostatic')
    dst_dy = os.path.join(destination,'hydrodynamic')
    
    try:
        shutil.copytree(src_st, dst_st)
        shutil.copytree(src_dy, dst_dy)
    except Exception as e:
            return (False, e)
    
    return (True, '')

def clean_prj_folder(path, exept=None):
    
    # remove all files and subfolders from the project folder
    folder_ls = os.listdir(path)
    if not exept is None:
        folder_ls = [el for el in folder_ls if el != exept]
    
    for el in folder_ls:
        el_path = os.path.join(path, el)
        try:
            if os.path.isfile(el_path):
                os.unlink(el_path)
            elif os.path.isdir(el_path):
                shutil.rmtree(el_path)
        except Exception as e:
            return (False, e)
            
    return (True, '') 

def split_string(string, num_type=float):
    list_el = string.split()
    list_out = []
    for el in list_el:
        try:
            list_out.append(num_type(el))
        except:
            break
    return list_out

def load_file(path, file_name):
    """
    the method is case sensitive, it distinguishes from capital and lower size letters
    :param path:
    :param file_name:
    :return:
    """
    check_file_name = [el for el in os.listdir(path) if el==file_name]
    if not check_file_name:
        raise NameError("The specified filename ({}) does not exist in the folder {}".format(file_name, path))
    elif len(check_file_name) > 1:  # this cndition cannot exist, just keep it here in case the comparison is bugged
        raise ValueError("The specified filename ({}) exists in multiple copies in the folder {}. ".format(file_name, path),
                         "Remove the unused files")

    with open(os.path.join(path, file_name)) as fid:
        lines = fid.readlines()

    return lines