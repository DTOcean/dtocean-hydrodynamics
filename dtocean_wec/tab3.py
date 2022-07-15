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
Created on Wed Jun 15 09:15:29 2016

.. moduleauthor:: Francesco Ferri <ff@civil.aau.dk>
"""
import os
import shutil

from PyQt4.QtCore import *
from PyQt4.QtGui import *

import utils.data_check as dck
from .form_utils import *
from .read_nemoh_form import Ui_Form as Ui_T3
from .submodule.utils import input_control as in_ck
from .utils.file_utilities import split_string, split_string_multilist


class ReadNemoh(QWidget, Ui_T3):
    trigger_results = pyqtSignal([dict])
    trigger_save = pyqtSignal([dict])
    trigger_reset_forms = pyqtSignal()
    trigger_mesh_view = pyqtSignal([dict])

    
    def __init__(self, parent=None):
        QWidget.__init__(self, parent)
        self.setupUi(self)
        self.btn_browse_t3.clicked.connect(self.browse_data)        
        self.btn_mesh_f_t3.clicked.connect(self.browse_mesh)
        self.btn_submit_t3.clicked.connect(self.submit_input_data)
        self.btn_load_data_t3.clicked.connect(self.load_data)
        self.btn_load_data_t3.setEnabled(False)
        
        QToolTip.setFont(QFont("SansSerif", 11))
        self.label_25.setToolTip("The number of independent dofs is the problem dimension with only generalised dofs.\nThe generalised dofs are obtained from the free body dofs after the\napplication of the mechanical constraints \n(only for multibody wec, e.g. RM3, Pelamis)")
        self.label_28.setToolTip("dofs associated with the energy extraction (damping only).\nFor the description of the index rule refer to the User Manual, section Degrees of Freedom")
        self.label_26.setToolTip("dofs associated with the mooring load (stiffness only).\nFor the description of the index rule refer to the User Manual, section Degrees of Freedom")
        
        for ish, chl in enumerate(self.groupBox_5.children()):
            if ish == 0: continue
            chl.stateChanged.connect(self.shared_dof_handles)
        
        self.btn_add_body_t3.clicked.connect(self.add_data_tab_body)
        self.btn_remove_body_t3.clicked.connect(self.remove_data_tab_body)
        self.trigger_results.connect(parent.set_hydrodynamic_results)
        self.trigger_save.connect(parent.task_save_hydrodynamic)
        self.trigger_reset_forms.connect(parent.task_reset_forms)
        
        self.db_folder = os.path.join(parent.wec_share_path, "wec_db")
        self.bin_folder = parent.bin_path
        
        self.trigger_mesh_view.connect(parent.task_show_mesh)
        
        self._raw = "raw_data"
        self._prj = parent._data["prj_folder"]
        if not os.path.isdir(os.path.join(self._prj, self._raw)):
            os.mkdir(os.path.join(self._prj, self._raw))
        
    def set_data(self, data):
        self._data = data
        self.populate_prj()
        self.btn_load_data_t3.setEnabled(False)
        
        if 'hyd' in data.keys():
            self.trigger_results.emit(self._data['hyd'])
        
    def browse_data(self):
        folder = QFileDialog.getExistingDirectory(self, "Select the folder that contains the Nemoh results and inputs")
        self.le_data_t3.setText(folder)
        
    def clear_tab_inputs(self):
        self.mesh_f_t2.setText('')
        self.id_body.setText('')
        self.cog_body.setText('')
        self.cs_body.setText('')
        self.parent_body.setText('')
        self.mass_body.setText('')
        self.inertia_body.setText('')
        self.dof_body.setText('')
        self.chil_body.setText('')
    
    def check_body_data(self):
        mesh_body = self.mesh_f_t3.text()
        if not os.path.isfile(mesh_body):
            mesh_body = os.path.join(self._prj, self._raw, os.path.basename(mesh_body))
        
        id_body = split_string(self.id_body_t3.text(), int)
        cog_body = split_string(self.cog_body_t3.text(), sep=',')
        cs_body = split_string(self.cs_body_t3.text(), sep=',')
        dof_body = split_string_multilist(self.dof_body_t3.text(), sep=',', sep_multi=';')
        chil_body = split_string_multilist(self.chil_body_t3.text(), sep=',', sep_multi=';')
        parent_body = split_string(self.parent_body_t3.text(), int)
        mass_body = split_string(self.mass_body_t3.text(), float)
        inertia_body = split_string_multilist(self.inertia_body_t3.text(), float, sep=',', sep_multi=';')

        
        data = [id_body, mass_body, inertia_body, mesh_body, cog_body, cs_body, dof_body, chil_body, parent_body]
        
        data_st = []
        data_st.append(raise_error([el for el in range(1) if os.path.isfile(mesh_body)], self.mesh_f_t3, 1, self.add_err_t3))
        data_st.append(raise_error(id_body, self.id_body_t3, 1, self.add_err_t3))
        data_st.append(raise_error(cog_body, self.cog_body_t3, 3, self.add_err_t3))
        data_st.append(raise_error(cs_body, self.cs_body_t3, 3, self.add_err_t3))
        data_st.append(raise_error(dof_body, self.dof_body_t3, 0, self.add_err_t3, gt=True))
        data_st.append(raise_error(chil_body, self.chil_body_t3, 0, self.add_err_t3, gt=True))
        data_st.append(raise_error(parent_body, self.parent_body_t3, 1, self.add_err_t3))
        data_st.append(raise_error(mass_body, self.mass_body_t3, 1, self.add_err_t3))
        data_st.append(raise_error([item for sublist in inertia_body for item in sublist], self.inertia_body_t3, 9, self.add_err_t3))
        
        return data, data_st
         
    def add_data_tab_body(self):
        # first check if the data is actually available in the edit line
        data, data_st = self.check_body_data()
        if -1 in data_st:
            pass
        else:
            #copy mesh file to a local raw_data folder in order to use relative path
            src = data[3]
            if not os.path.split(src)[0]=="raw_data":
                dst = os.path.join(self._prj, self._raw)
                if not dst==os.path.split(src)[0]:
                    shutil.copy(src, dst)
                
            fn = str(os.path.basename(src))    
            
            row_position = self.tab_body_t3.rowCount()
            self.tab_body_t3.insertRow(row_position)
            self.tab_body_t3.setItem(row_position, 0, QTableWidgetItem("".join([str(x) for x in data[0]])))
            self.tab_body_t3.setItem(row_position, 1, QTableWidgetItem("".join([str(x) for x in data[1]])))
            self.tab_body_t3.setItem(row_position, 2, QTableWidgetItem(";".join([str(x) for x in data[2]]).replace('[','').replace(']','')))
            self.tab_body_t3.setItem(row_position, 3, QTableWidgetItem(fn))  # data[3]
            self.tab_body_t3.setItem(row_position, 4, QTableWidgetItem(",".join([str(x) for x in data[4]])))
            self.tab_body_t3.setItem(row_position, 5, QTableWidgetItem(",".join([str(x) for x in data[5]])))
            self.tab_body_t3.setItem(row_position, 6, QTableWidgetItem(";".join([str(x) for x in data[6]]).replace('[','').replace(']','')))
            self.tab_body_t3.setItem(row_position, 7, QTableWidgetItem(";".join([str(x) for x in data[7]]).replace('[','').replace(']','')))
            self.tab_body_t3.setItem(row_position, 8, QTableWidgetItem("".join([str(x) for x in data[8]])))
    
    def remove_data_tab_body(self, clear_all=False):
        if clear_all:
            self.tab_body_t3.setRowCount(0)
        else:
            rows = sorted(set(index.row() for index in
                          self.tab_body_t3.selectedIndexes()))
            for row in rows:
                self.tab_body_t3.removeRow(row)
        
    def load_data(self):
        self.btn_load_data_t3.setEnabled(False)
        self.btn_submit_t3.setEnabled(False)
        
        stat = send_data_to_bem_interface(self._data,
                                          self.db_folder,
                                          self.bin_folder)
        
        self.btn_submit_t3.setEnabled(True)

        if stat[0]:
            self._data['hyd'] = stat[1]
            self.trigger_results.emit(self._data['hyd'])

            
    def check_input_data(self):
        folder_stat = in_ck.check_nemoh_results(self.le_data_t3.text(), get_array_mat=self.cb_gen_array_mat_t3.isChecked())
        if folder_stat[0]:
            status = dck.check_inputs_hydrodynamic(self._data['inputs_hydrodynamic'], self._prj)
            if status:
                print('\n'.join([x for x in status]))
                return (False, )
            else:
                return (True, )
        else:
            print('\n'.join([x for x in folder_stat[1]]))
            return (False, )
            
    def browse_mesh(self):
        folder = QFileDialog.getOpenFileName(self, "Select the mesh file of the body")
        self.mesh_f_t3.setText(folder)
    
    def show_mesh(self):
        filename = str(self.mesh_f_t3.text())
        if not os.path.isfile(filename):
            fn = os.path.basename(filename)
            filename = os.path.join(self._prj, self._raw, fn)
            if not os.path.isfile(filename):
                self.browse_folder()
			

        path, file_name = os.path.split(filename)
        self.trigger_mesh_view.emit({'path':path,'f_n':file_name})
    
    def submit_input_data(self):
        self.save_project()
        status = self.check_input_data()
        if status[0]:
            self.btn_load_data_t3.setEnabled(True)
            self._data['inputs_hydrodynamic']['general_input']['ndof'] = split_string(self.ndof_t3.text(), int)
            self._data['inputs_hydrodynamic']['general_input']['pto_dof'] = split_string(self.pto_dof_t3.text(), int)
            self._data['inputs_hydrodynamic']['general_input']['mooring_dof'] = split_string(self.moor_dof_t3.text(), int)
            self._data['inputs_hydrodynamic']['general_input']['input_type'] = 3
            self._data['inputs_hydrodynamic']['general_input']['data_folder'] = self.le_data_t3.text()
            self._data['inputs_hydrodynamic']['general_input']['get_array_mat'] = int(self.cb_gen_array_mat_t3.isChecked())
            
        else:
            self.btn_load_data_t3.setEnabled(False)
    
    def pull_data_from_form(self):
        bodies = {}
        for el in range(self.tab_body_t3.rowCount()):
            
            id_body = split_string(self.tab_body_t3.item(el,0).text(), int)
            mass_body = split_string(self.tab_body_t3.item(el,1).text())
            inertia_body = split_string_multilist(self.tab_body_t3.item(el,2).text(), float, sep=',', sep_multi=';')
            mesh_body = self.tab_body_t3.item(el,3).text()
            
            cog_body = split_string(self.tab_body_t3.item(el,4).text(), float, sep=',')
            cs_body = split_string(self.tab_body_t3.item(el,5).text(), float, sep=',')
            dof_body = split_string_multilist(self.tab_body_t3.item(el,6).text(), float, sep=',', sep_multi=';')
            if not dof_body[0]:
                dof_body = -1
            chil_body = split_string_multilist(self.tab_body_t3.item(el,7).text(), float, sep=',', sep_multi=';')
            if not chil_body[0]:
                chil_body = -1
            parent_body = split_string(self.tab_body_t3.item(el,8).text(), int)
            
            str_bi = 'body{}'.format(id_body[0])
            bodies[str_bi] =  {'mass': mass_body, 'inertia': inertia_body, 'ID': id_body, 'mesh': mesh_body,
                                                    'cog': cog_body, 'child_dof_position': chil_body,
                                                    'dof_with_parent': dof_body, 'number_child':1,
                                                    'axis_angles':cs_body, 'parent_body': parent_body}
            

        ndof = split_string(self.ndof_t3.text(), int)
        pto_dof = split_string(self.pto_dof_t3.text(), int, sep=',')
        moor_dof = split_string(self.moor_dof_t3.text(), int)
        local_cs = split_string(self.local_cs_t3.text(), float, sep=',')
        
        shared_dof = list(self._data['inputs_hydrodynamic']['body_inputs']['shared_dof'])
        get_array_mat = self.cb_gen_array_mat_t3.isChecked()
       
        input_type = 3
        data_f = self.le_data_t3.text()

        general_inputs = {'ndof': ndof, 'pto_dof': pto_dof, 'mooring_dof': moor_dof,
                          'frequency_def': [], 'angle_def': [],  'data_folder': data_f,
                          'water_depth': [], 'get_array_mat': int(get_array_mat), 
                          'cyl_nzeta': [], 'cyl_ntheta': [],'input_type':input_type}
        
        
        body_inputs = {'local_cs': local_cs, 'shared_dof': shared_dof, 'body': bodies}
        inputs_hydrodynamic = {'general_input': general_inputs, 'body_inputs': body_inputs}
        
        return inputs_hydrodynamic
        
    def save_project(self):
        inputs_hydrodynamic = self.pull_data_from_form()
        
        dic_cmp = compare_dictionaries(self._data['inputs_hydrodynamic'], inputs_hydrodynamic,'d1','d2')
        
        if dic_cmp != "":
            self._data['inputs_hydrodynamic'] = inputs_hydrodynamic
            if 'hyd' in self._data.keys():
                print("the hydrodynamic inputs have been modified, the hydrodynamic results and the power fit results will be deleted!")
                del self._data['hyd']
            self.trigger_save.emit(inputs_hydrodynamic)
            self.trigger_reset_forms.emit()
            
        else:
            self.trigger_save.emit(inputs_hydrodynamic)
            if 'hyd' in self._data.keys():
                self.trigger_results.emit(self._data['hyd'])

        
    def shared_dof_handles(self, checked):
        checkBox_ind = int(self.sender().objectName()[2:-3])
        if checked:
            self._data['inputs_hydrodynamic']['body_inputs']['shared_dof'][checkBox_ind] = 1
        else:
            self._data['inputs_hydrodynamic']['body_inputs']['shared_dof'][checkBox_ind] = 0
        
    def populate_prj(self):
        
        in_hy = self._data['inputs_hydrodynamic']
        in_g = in_hy['general_input']
        self.le_data_t3.setText(in_g['data_folder'])
        self.ndof_t3.setText(",".join([str(el) for el in in_g['ndof']]))
        self.pto_dof_t3.setText(",".join([str(el) for el in in_g['pto_dof']]))
        self.moor_dof_t3.setText(",".join([str(el) for el in in_g['mooring_dof']]))
        
        in_b = in_hy['body_inputs']
        self.local_cs_t3.setText(",".join([str(el) for el in in_b['local_cs']]))
        shared_dof = in_b['shared_dof']
        sh_ch = self.groupBox_5.children()
        for iel, el in enumerate(shared_dof):
            if el==1:
                sh_ch[iel+1].setChecked(True)
            else:
                sh_ch[iel+1].setChecked(False)
        
        bodies = []
        if 'body' in in_b.keys():
            bodies = in_b['body']
        if bodies == -1:
            pass
        else:
            for bi in bodies:
                di = bodies[bi]
                
                mesh_f = di['mesh']
                fn = os.path.basename(mesh_f)                    
                
                self.mesh_f_t3.setText(fn)
                self.id_body_t3.setText(",".join([str(el) for el in di['ID']]))
                self.cog_body_t3.setText(",".join([str(el) for el in di['cog']]))
                self.cs_body_t3.setText(",".join([str(el) for el in di['axis_angles']]))
                self.parent_body_t3.setText(",".join([str(el) for el in di['parent_body']]))
                self.mass_body_t3.setText(",".join([str(el) for el in di['mass']]))
                
                dof = di['dof_with_parent'].tolist()
                string = ''
                if not dof == -1:
                    for iel, el in enumerate(dof):
                        if not len(el)==0:
                            if not iel == 0: string += ';' 
                            string += ",".join(str(el) for el in el)          
                self.dof_body_t3.setText(string)
                    
                ch = di['child_dof_position'].tolist()
                string = ''
                if not ch == -1:
                    for iel, el in enumerate(ch):
                        if not len(el)==0:
                            if not iel == 0: string += ';' 
                            string += ",".join(str(el) for el in el)
                self.chil_body_t3.setText(string)
                inertia_ls = di['inertia'].tolist()
                self.inertia_body_t3.setText(";".join([str(el) for el in inertia_ls]).replace('[', '').replace(']',''))
                                            
                self.add_data_tab_body()

def compare_dictionaries(dict_1, dict_2, dict_1_name, dict_2_name, path=""):
    """Compare two dictionaries recursively to find non mathcing elements

    Args:
        dict_1: dictionary 1
        dict_2: dictionary 2

    Returns:
    for key in keys:
    val1, val2 = ax1[key], ax2[key]
    are_different = val1 != val2
    if isinstance(val1, np.ndarray):
        are_different = are_different.any()

    if are_different:
        print(key,ax1[key])

    """
    err = ''
    key_err = ''
    value_err = ''
    old_path = path
    for k in dict_1.keys():
        path = old_path + "[%s]" % k
        if not dict_2.has_key(k):
            key_err += "Key %s%s not in %s\n" % (dict_2_name, path, dict_2_name)
        else:
            if isinstance(dict_1[k], dict) and isinstance(dict_2[k], dict):
                err += compare_dictionaries(dict_1[k],dict_2[k],'d1','d2', path)
            elif isinstance(dict_1[k], np.ndarray):
                if not np.all(dict_1[k] == dict_2[k]):
                    value_err += "Value of %s%s (%s) not same as %s%s (%s)\n"\
                        % (dict_1_name, path, dict_1[k], dict_2_name, path, dict_2[k])
            else:
                if dict_1[k] != dict_2[k]:
                    value_err += "Value of %s%s (%s) not same as %s%s (%s)\n"\
                        % (dict_1_name, path, dict_1[k], dict_2_name, path, dict_2[k])

    for k in dict_2.keys():
        path = old_path + "[%s]" % k
        if not dict_1.has_key(k):
            key_err += "Key %s%s not in %s\n" % (dict_2_name, path, dict_1_name)

    return key_err + value_err + err