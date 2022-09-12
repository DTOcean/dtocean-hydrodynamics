# -*- coding: utf-8 -*-

#    Copyright (C) 2016 Francesco Ferri
#    Copyright (C) 2017-2018 Mathew Topper
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
Created on Wed Jun 15 09:15:30 2016

.. moduleauthor:: Francesco Ferri <ff@civil.aau.dk>
.. moduleauthor:: Mathew Topper <mathew.topper@dataonlygreater.com>
"""

import os
import shutil

from PyQt4.QtCore import *
from PyQt4.QtGui import *

import utils.data_check as dck
from .form_utils import *
from .run_nemoh_form import Ui_Form as Ui_T2
from .utils.file_utilities import split_string, split_string_multilist


class RunNemoh(QWidget, Ui_T2):
    trigger_results = pyqtSignal([dict])
    trigger_save = pyqtSignal([dict])
    trigger_reset_forms = pyqtSignal()
    trigger_mesh_view = pyqtSignal([dict])
    
    def __init__(self, parent=None):
        QWidget.__init__(self, parent)
        self.setupUi(self)
        
        self.btn_mesh_f_t2.clicked.connect(self.browse_folder)
        self.btn_view_mesh.clicked.connect(self.show_mesh)
        self.btn_submit_t2.clicked.connect(self.submit_input_data)
        
        self.btn_add_body.clicked.connect(self.add_data_tab_body)
        self.btn_remove_body.clicked.connect(self.remove_data_tab_body)
        self.btn_calculate_t2.clicked.connect(self.link_to_bem_interface)
        
        QToolTip.setFont(QFont("SansSerif", 11))
        self.label_2.setToolTip("The number of independent dofs is the problem dimension with only generalised dofs.\nThe generalised dofs are obtained from the free body dofs after the\napplication of the mechanical constraints \n(only for multibody wec, e.g. RM3, Pelamis)")
        self.label_3.setToolTip("dofs associated with the energy extraction (damping only).\nFor the description of the index rule refer to the User Manual, section Degrees of Freedom")
        self.label_4.setToolTip("dofs associated with the mooring load (stiffness only).\nFor the description of the index rule refer to the User Manual, section Degrees of Freedom")

        for ish, chl in enumerate(self.groupBox_2.children()):
            if ish == 0: continue
            chl.stateChanged.connect(self.shared_dof_handles)
    
        self.btn_calculate_t2.setEnabled(False)
        self.btn_stop_calculation_t2.setEnabled(False)
              
        
        self.cb_gen_array_mat.stateChanged.connect(self.gen_array_handles)
        self.gr_cyl_spec.setEnabled(False)        
        
        # self._data = parent._data
        self.trigger_results.connect(parent.set_hydrodynamic_results)
        self.trigger_save.connect(parent.task_save_hydrodynamic)
        self.trigger_reset_forms.connect(parent.task_reset_forms)
        self.trigger_mesh_view.connect(parent.task_show_mesh)
        self.db_folder = os.path.join(parent.wec_share_path, "wec_db")
        self.bin_folder = parent.bin_path
        
        self._raw = "raw_data"
        self._prj = parent._data["prj_folder"]
        if not os.path.isdir(os.path.join(self._prj, self._raw)):
            os.mkdir(os.path.join(self._prj, self._raw))
         
		
    def show_mesh(self):
        filename = str(self.mesh_f_t2.text())
        if not os.path.isfile(filename):
            fn = os.path.basename(filename)
            filename = os.path.join(self._prj, self._raw, fn)
            if not os.path.isfile(filename):
                self.browse_folder()
			

        path, file_name = os.path.split(filename)
        self.trigger_mesh_view.emit({'path':path,'f_n':file_name})
	
    def set_data(self, data):
        self._data = data
        self.populate_prj()
        self.btn_calculate_t2.setEnabled(False)
        self.btn_stop_calculation_t2.setEnabled(False)
        
        if 'hyd' in data.keys():
            self.trigger_results.emit(self._data['hyd'])

    def gen_array_handles(self, checked):
        if checked:
            self.gr_cyl_spec.setEnabled(True)
        else:
            self.gr_cyl_spec.setEnabled(False) 
            
    def link_to_bem_interface(self):
        read_nemoh_flag = False
        if 'hydrodynamic' in os.listdir(self._data['prj_folder']):
            msgBox = QMessageBox()
            msgBox.setText('The project folder already contains a BEM result folder. Do you want to overwrite it, load another project or start a new project?')
            msgBox.addButton(QPushButton('Overwrite'), QMessageBox.YesRole)
            msgBox.addButton(QPushButton('Load'), QMessageBox.NoRole)
            msgBox.addButton(QPushButton('New'), QMessageBox.RejectRole)
            ret = msgBox.exec_();
            
            if ret == 0:
                try:
                    prj_fn = self._data['prj_filename']
                    clean_prj_folder(self._data['prj_folder'], exept=[prj_fn, "raw_data"])
                except:
                    print(Exception)
                    
            elif ret == 1:
                # Force the read_nemoh option in thge BemSolution class by changing the run flag
                # in the _data dictionary.
                QMessageBox.information(self, "Information", "It is important that the data contained in the hydrodynamic folder is in agreement with the gui input.")
                read_nemoh_flag = True
            else:
                print("Not implemented yet, use the create new project under the File menu!")
        
        self.btn_calculate_t2.setEnabled(False)
        self.btn_stop_calculation_t2.setEnabled(False)
        stat = send_data_to_bem_interface(self._data,
                                          self.db_folder,
                                          self.bin_folder,
                                          force_read_flag=read_nemoh_flag)
        if stat[0]:
            self._data['hyd'] = stat[1]
            self.btn_submit_t2.setEnabled(True)
            self.trigger_results.emit(stat[1])
    
    def load_nemoh_solution():
        """
        
        """
    
    def browse_folder(self):
        folder = QFileDialog.getOpenFileName(
                self,
                "Select the mesh file of the body",
                "",
                "mesh files (*.dat *.gdf *.GDF)")
        self.mesh_f_t2.setText(folder)
        
    def submit_input_data(self):
        self.save_project()
        status = self.check_input_data()
        if status:
            self.btn_calculate_t2.setEnabled(True)
            self.btn_stop_calculation_t2.setEnabled(True)
        else:
            self.btn_calculate_t2.setEnabled(False)
            self.btn_stop_calculation_t2.setEnabled(False)
            
    def check_input_data(self):
        status = dck.check_inputs_hydrodynamic(self._data['inputs_hydrodynamic'], self._prj)
        if status:
            print(' \n '.join([x for x in status]))
            return False
        else:
            return True
            
            
    def pull_data_from_form(self):
        bodies = {}
        for el in range(self.tab_body.rowCount()):
            
            id_body = split_string(self.tab_body.item(el,0).text(), int)
            mass_body = split_string(self.tab_body.item(el,1).text())
            inertia_body = split_string_multilist(self.tab_body.item(el,2).text(), float, sep=',', sep_multi=';')
            mesh_body = str(self.tab_body.item(el,3).text())
            
            cog_body = split_string(self.tab_body.item(el,4).text(), float, sep=',')
            cs_body = split_string(self.tab_body.item(el,5).text(), float, sep=',')
            dof_body = split_string_multilist(self.tab_body.item(el,6).text(), float, sep=',', sep_multi=';')
            if not dof_body[0]:
                dof_body = -1
            chil_body = split_string_multilist(self.tab_body.item(el,7).text(), float, sep=',', sep_multi=';')
            if not chil_body[0]:
                chil_body = -1
            parent_body = split_string(self.tab_body.item(el,8).text(), int)
            
            str_bi = 'body{}'.format(id_body[0])
            bodies[str_bi] =  {'mass': mass_body, 'inertia': inertia_body, 'ID': id_body, 'mesh': mesh_body,
                                                    'cog': cog_body, 'child_dof_position': chil_body,
                                                    'dof_with_parent': dof_body, 'number_child':1,
                                                    'axis_angles':cs_body, 'parent_body': parent_body}
            

        ndof = split_string(self.ndof.text(), int)
        pto_dof = split_string(self.pto_dof.text(), int, sep=',')
        moor_dof = split_string(self.moor_dof.text(), int)
        freq_def = split_string(self.fre_def.text(), float, sep=',')
        angle_def = split_string(self.angles_def.text(), float)
        local_cs = split_string(self.local_cs.text(), float, sep=',')
        shared_dof = list(self._data['inputs_hydrodynamic']['body_inputs']['shared_dof'])
        water_depth = self.sb_water_depth.value()
        get_array_mat = self.cb_gen_array_mat.isChecked()
        nzeta, ntheta = (0,0)
        if get_array_mat:
            ntheta = self.sb_ntheta.value()
            nzeta = self.sb_nzeta.value()
        input_type = 2
        
        general_inputs = {'ndof': ndof, 'pto_dof': pto_dof, 'mooring_dof': moor_dof,
                          'frequency_def': freq_def, 'angle_def': angle_def, 
                          'water_depth': water_depth, 'get_array_mat': int(get_array_mat), 
                          'cyl_nzeta': nzeta, 'cyl_ntheta': ntheta,'input_type':input_type}
        
        
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
        checkBox_ind = int(self.sender().objectName()[-1])
        if checked:
            self._data['inputs_hydrodynamic']['body_inputs']['shared_dof'][checkBox_ind] = 1
        else:
            self._data['inputs_hydrodynamic']['body_inputs']['shared_dof'][checkBox_ind] = 0
        

    def check_body_data(self):
        mesh_body = str(self.mesh_f_t2.text())
        if not os.path.isfile(mesh_body):
            mesh_body = os.path.join(self._prj, self._raw, os.path.basename(mesh_body))
        
        id_body = split_string(self.id_body.text(), int)
        cog_body = split_string(self.cog_body.text(), sep=',')
        cs_body = split_string(self.cs_body.text(), sep=',')
        dof_body = split_string_multilist(self.dof_body.text(), sep=',', sep_multi=';')
        chil_body = split_string_multilist(self.chil_body.text(), sep=',', sep_multi=';')
        parent_body = split_string(self.parent_body.text(), int)
        mass_body = split_string(self.mass_body.text(), float)
        inertia_body = split_string_multilist(self.inertia_body.text(), float, sep=',', sep_multi=';')

        
        data = [id_body, mass_body, inertia_body, mesh_body, cog_body, cs_body, dof_body, chil_body, parent_body]
        
        data_st = []
        data_st.append(self.__raise_error([el for el in range(1) if os.path.isfile(os.path.join(self._prj, mesh_body))], self.mesh_f_t2, 1, self.add_err_t2))
        data_st.append(self.__raise_error(id_body, self.id_body, 1, self.add_err_t2))
        data_st.append(self.__raise_error(cog_body, self.cog_body, 3, self.add_err_t2))
        data_st.append(self.__raise_error(cs_body, self.cs_body, 3, self.add_err_t2))
        data_st.append(self.__raise_error(dof_body, self.dof_body, 0, self.add_err_t2, gt=True))
        data_st.append(self.__raise_error(chil_body, self.chil_body, 0, self.add_err_t2, gt=True))
        data_st.append(self.__raise_error(parent_body, self.parent_body, 1, self.add_err_t2))
        data_st.append(self.__raise_error(mass_body, self.mass_body, 1, self.add_err_t2))
        data_st.append(self.__raise_error([item for sublist in inertia_body for item in sublist], self.inertia_body, 9, self.add_err_t2))
        
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
            
            row_position = self.tab_body.rowCount()
            self.tab_body.insertRow(row_position)
            self.tab_body.setItem(row_position, 0, QTableWidgetItem("".join([str(x) for x in data[0]])))
            self.tab_body.setItem(row_position, 1, QTableWidgetItem("".join([str(x) for x in data[1]])))
            self.tab_body.setItem(row_position, 2, QTableWidgetItem(";".join([str(x) for x in data[2]]).replace('[','').replace(']','')))
            self.tab_body.setItem(row_position, 3, QTableWidgetItem(fn))  # data[3]
            self.tab_body.setItem(row_position, 4, QTableWidgetItem(",".join([str(x) for x in data[4]])))
            self.tab_body.setItem(row_position, 5, QTableWidgetItem(",".join([str(x) for x in data[5]])))
            self.tab_body.setItem(row_position, 6, QTableWidgetItem(";".join([str(x) for x in data[6]]).replace('[','').replace(']','')))
            self.tab_body.setItem(row_position, 7, QTableWidgetItem(";".join([str(x) for x in data[7]]).replace('[','').replace(']','')))
            self.tab_body.setItem(row_position, 8, QTableWidgetItem("".join([str(x) for x in data[8]])))
            
    
    def remove_data_tab_body(self, clear_all=False):
        if clear_all:
            self.tab_body.setRowCount(0)
        else:
            rows = sorted(set(index.row() for index in
                          self.tab_body.selectedIndexes()))
            for row in rows:
                self.tab_body.removeRow(row)
            
    def __raise_error(self, el, obj, el_sh, label_err, gt=True):
        label_err.setText("")
        label_err.setStyleSheet("QLabel {color:black}")
        obj.setStyleSheet("QLineEdit {border: 0px}")
            
        if not el:
            label_err.setText("invalid inputs")
            label_err.setStyleSheet("QLabel {color:red}")
            obj.setStyleSheet("QLineEdit {border: 1px solid red}")
            return -1
            
        elif not len(el) == el_sh:
            if gt:
                if not len(el) >= el_sh:
                    label_err.setText("invalid inputs")
                    label_err.setStyleSheet("QLabel {color:red}")
                    obj.setStyleSheet("QLineEdit {border: 1px solid red}")
                    
                    return -1
        return 0
    
    def populate_prj(self):
        in_hy = self._data['inputs_hydrodynamic']
        in_g = in_hy['general_input']
        self.ndof.setText(",".join([str(el) for el in in_g['ndof']]))
        self.pto_dof.setText(",".join([str(el) for el in in_g['pto_dof']]))
        self.moor_dof.setText(",".join([str(el) for el in in_g['mooring_dof']]))
        self.fre_def.setText(",".join([str(el) for el in in_g['frequency_def']]))
        self.angles_def.setText(",".join([str(el) for iel, el in enumerate(in_g['angle_def']) if iel==0]))
        self.sb_water_depth.setValue(in_g['water_depth'])
        
        self.cb_gen_array_mat.setChecked(in_g['get_array_mat']==1)
        if self.cb_gen_array_mat.isChecked():
            self.sb_ntheta.setValue(in_g['cyl_ntheta'])
            self.sb_nzeta.setValue(in_g['cyl_nzeta'])
        in_b = in_hy['body_inputs']
        self.local_cs.setText(",".join([str(el) for el in in_b['local_cs']]))
        shared_dof = in_b['shared_dof']
        sh_ch = self.groupBox_2.children()
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
                self.mesh_f_t2.setText(fn)
                self.id_body.setText(",".join([str(el) for el in di['ID']]))
                self.cog_body.setText(",".join([str(el) for el in di['cog']]))
                self.cs_body.setText(",".join([str(el) for el in di['axis_angles']]))
                self.parent_body.setText(",".join([str(el) for el in di['parent_body']]))
                self.mass_body.setText(",".join([str(el) for el in di['mass']]))
                
                dof = di['dof_with_parent'].tolist()
                string = ''
                if not dof == -1:
                    for iel, el in enumerate(dof):
                        if not len(el)==0:
                            if not iel == 0: string += ';' 
                            string += ",".join(str(el) for el in el)          
                self.dof_body.setText(string)
                    
                ch = di['child_dof_position'].tolist()
                string = ''
                if not ch == -1:
                    for iel, el in enumerate(ch):
                        if not len(el)==0:
                            if not iel == 0: string += ';' 
                            string += ",".join(str(el) for el in el)
                self.chil_body.setText(string)
                inertia_ls = di['inertia'].tolist()
                self.inertia_body.setText(";".join([str(el) for el in inertia_ls]).replace('[', '').replace(']',''))
                                            
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


