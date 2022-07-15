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
Created on Wed Jun 15 09:16:39 2016

.. moduleauthor:: Francesco Ferri <ff@civil.aau.dk>
"""
import os

from PyQt4.QtCore import *
from PyQt4.QtGui import *

import dtocean_wave.utils.hdf5_interface as h5i

from .form_utils import *
from .read_db_form import Ui_Form as Ui_T1


class ReadDb(QWidget, Ui_T1):
    trigger_results = pyqtSignal([dict])
    trigger_save = pyqtSignal([dict])
    trigger_save_prj = pyqtSignal()
    
    def __init__(self, parent=None):
        QWidget.__init__(self, parent)
        self.setupUi(self)
        
        self.btn_load_data_t1.clicked.connect(self.load_data)
        self._wec_db_folder = os.path.join(parent.wec_share_path, "wec_db")
        self.trigger_results.connect(parent.set_wec_db_results)
        self.trigger_save.connect(parent.task_save_hydrodynamic)
        self.trigger_save_prj.connect(parent.save_project)
        self._init_flag = True
        
    def set_data(self, data):
        if not 'hyd' in data.keys():
            self._init_flag = True
        else:
            self._init_flag = False
            
        self._data = data
        self.populate_prj()
        self.btn_load_data_t1.setEnabled(True)
        
        if 'hyd' in data.keys() and not self._init_flag :
            self.trigger_results.emit(self._data)
        
    def load_data(self):
        current_case = str(self.cb_wec_t1.currentText()).replace(' ', '_')
        file_name = os.path.join(self._wec_db_folder,
                                 current_case,
                                 '{}_results.hdf5'.format(current_case))
        description_file = os.path.join(self._wec_db_folder,
                                        current_case,
                                        '{}_description.html'.format(current_case)).replace("\\", "/")
        data_db = h5i.load_dict_from_hdf5(file_name)
        
        data_db['prj_filename'] = self._data['prj_filename']
        data_db['prj_folder'] = self._data['prj_folder']
        data_db['prj_name'] = self._data['prj_name']
        data_db['inputs_hydrodynamic']['general_input']['data_folder'] = os.path.join(self._wec_db_folder, current_case)
        print(os.path.join(self._wec_db_folder, current_case))
        stat = check_wec_db_data(data_db)
        self.save_project()
        if stat[0]: 
            self._data = data_db
            self.ww_wec.load(QUrl("file:///"+description_file))
            print("Data loaded. It is possible to plot the data.")
            self.trigger_save.emit(self._data['inputs_hydrodynamic'])
            self.trigger_results.emit(self._data)
        else:
            print(stat[1])
        
    def save_project(self):
        current_case = str(self.cb_wec_t1.currentText()).replace(' ', '_')
        data_f = os.path.join(self._wec_db_folder,
                              current_case)
       
        input_type = 1
                
        self._data['inputs_hydrodynamic']['general_input']['input_type'] = input_type
        self._data['inputs_hydrodynamic']['general_input']['data_folder'] = data_f
        
        self.trigger_save.emit(self._data['inputs_hydrodynamic'])
        
    
    def populate_prj(self):
        folder = self._wec_db_folder
        sub_ls = []
        if os.path.isdir(folder):
            for sub in os.walk(folder).next()[1]:
                sub_a = os.path.join(folder, sub)
                if '{}_results.hdf5'.format(sub) in os.listdir(sub_a) and '{}_description.html'.format(sub) in os.listdir(sub_a):
                    sub_ls.append(sub)
            if not sub_ls:
                print("Ops, something went wrong. The wec db folder seems not valid anymore. Check the src folder path")
        
        if not self._init_flag:
            data_folder = self._data['inputs_hydrodynamic']['general_input']['data_folder']
            if os.path.isdir(data_folder):
                case = os.path.basename(os.path.normpath(data_folder))
                if '{}_results.hdf5'.format(case) in os.listdir(data_folder) and '{}_description.html'.format(case) in os.listdir(data_folder):
                    if case in sub_ls:
                        sub_ls.remove(case)
                        sub_ls.insert(0, case)
                
                
        self.cb_wec_t1.clear()
        self.cb_wec_t1.addItems(sub_ls)
        
        if not self._init_flag:
            self.load_data()