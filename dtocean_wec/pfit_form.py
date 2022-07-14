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
Created on Wed Jun 15 09:19:49 2016

.. moduleauthor:: Francesco Ferri <ff@civil.aau.dk>
.. moduleauthor:: Mathew Topper <mathew.topper@dataonlygreater.com>
"""

import os
import csv

import numpy as np
from matplotlib.figure import Figure
from matplotlib.backends.backend_qt4agg import (
    FigureCanvasQTAgg as FigureCanvas,
    NavigationToolbar2QT as NavigationToolbar)
from PyQt4.QtCore import *
from PyQt4.QtGui import *

import dtocean_wave.utils.hdf5_interface as h5i

from .performance_fit import Ui_Form as Ui_Fit
from .plot_performance_fit import VisualisePerformanceFit
from .submodule.power_matrix_fitting import PowerMatrixFit


class ActionPatch():
    def __init__(self, action):
        self.action = action
    
    def text(self):
        return self.action

class PowerPerformance(QWidget, Ui_Fit):
    trigger_results = pyqtSignal([dict])
    trigger_save = pyqtSignal([dict])
    
    def __init__(self, parent=None):
        QWidget.__init__(self, parent)
        self.setupUi(self)
        
        self._data = None
        self.btn_browse_pfit.clicked.connect(self.browse_data)
        self.btn_load_pfit.clicked.connect(self.load_data)
        self.btn_fitting.clicked.connect(self.fitting)
        self.btn_skip_fit.clicked.connect(self.skip_fitting)
        self.btn_plot_pfit.clicked.connect(self.plot_data)
        
        self.wwAngle.valueChanged.connect(self.__setYaw)
        self.uiSpreading.valueChanged.connect(self.__setSpreading)
        self.uiGamma.valueChanged.connect(self.__setGamma)
        
        self.btn_fitting.setEnabled(False)
        self.btn_skip_fit.setEnabled(False)
        self.btn_plot_pfit.setEnabled(False)
        self.btn_load_pfit.setEnabled(False)
        
        self.plotter = VisualisePerformanceFit()
        
        self.start_figure()
        
        self.trigger_results.connect(parent.task_end_pfit)  
        self.trigger_save.connect(parent.set_pfit_data)
        
        self.db_folder = os.path.join(parent.wec_share_path, "wec_db")
        
        self.__yaw = 0.0
        self.__spectrum = "Jonswap"
        self.__gamma = 1.0
        self.__spreading = -1
    
    def __setSpreading(self):
        self.__spreading = self.uiSpreading.value()
    
    def __setGamma(self):
        self.__gamma = self.uiGamma.value()        
    
    def __setYaw(self):
        self.__yaw = self.wwAngle.value()
    
    def closeEvent(self, event): 
        print "Closing" 
        event.accept()
       
    def start_figure(self):
        fig = Figure()
        self.addmpl(fig)
        
    def rmmpl(self,):
        self.mpl_vl.removeWidget(self.canvas)
        self.canvas.close()
        self.mpl_vl.removeWidget(self.toolbar)
        self.toolbar.close()
        
    def addmpl(self, fig):
        self.canvas = FigureCanvas(fig)
        self.mpl_vl.addWidget(self.canvas)
        self.canvas.draw()
        self.toolbar = NavigationToolbar(self.canvas, self.mpl_window, coordinates=True)
        self.mpl_vl.addWidget(self.toolbar)
        
    def set_dimensions(self, angle):
        angle = angle.tolist()
        angle = [str(x) for x in angle]
       
        self.cb_angle_index.clear()
        self.cb_angle_index.addItems(angle)
        self.plotter.set_data(self._data['inputs_pm_fit'])
        self.rb_sd.setChecked(True)
        
    def plot_data(self):
        # check which tab is active first
        active_rb = [ix for ix, x in enumerate(self.groupBox.children()) if ix>0 and x.isChecked()]
        
        if active_rb[0] == 1:
            self.plot_scatter_diagram()
        elif active_rb[0] == 2:
            self.plot_power_matrix()
        elif active_rb[0] == 3:
            self.plot_pto()
        elif active_rb[0] == 4:
            self.plot_mooring()
        elif active_rb[0] == 5:
            self.plot_ext_d()
        elif active_rb[0] == 6:
            self.plot_ext_k()
        elif active_rb[0] == 7:
            pass
        
    def plot_pto(self):
        fig = self.plotter.show_pto()
        self.rmmpl()
        self.addmpl(fig)
        
    def plot_mooring(self):
        fig = self.plotter.show_mooring()
        self.rmmpl()
        self.addmpl(fig)
        
    def plot_ext_d(self):
        fig = self.plotter.show_ext_d()
        self.rmmpl()
        self.addmpl(fig)
        
    def plot_ext_k(self):
        fig = self.plotter.show_ext_k()
        self.rmmpl()
        self.addmpl(fig)
        
    def plot_scatter_diagram(self):
        a_i  = self.cb_angle_index.currentIndex()
        fig = self.plotter.show_scatter_diagram(a_i)
        self.rmmpl()
        self.addmpl(fig)
        
    def plot_power_matrix(self):
        a_i  = self.cb_angle_index.currentIndex()
        fig = self.plotter.show_power_matrix(a_i)
        self.rmmpl()
        self.addmpl(fig)
        
    def save_data(self):
        prj_folder = self._data['prj_folder']
        prj_fn = self._data['prj_filename']
        project = os.path.join(prj_folder, prj_fn)
        h5i.save_dict_to_hdf5(self._data, project)
                
    def fitting(self):
        self.btn_fitting.setEnabled(False)
        self.btn_skip_fit.setEnabled(False)       
        per_fit = PowerMatrixFit(self._data['hyd']['periods'],
                                 self._data['hyd']['directions'],
                                self._data['hyd']['m_m'],
                                self._data['hyd']['m_add'],
                                self._data['hyd']['c_rad'],
                                self._data['inputs_pm_fit']['machine']['ext_d'],
                                self._data['hyd']['k_hst'],
                                self._data['inputs_pm_fit']['machine']['ext_k'],
                                self._data['hyd']['f_ex'],
                                self._data['hyd']['force_tr_mat'],
                                [el-1 for el in self._data['hyd']['pto_dof']],
                                [el-1 for el in self._data['hyd']['mooring_dof']],
                                self._data['hyd']['order'],
                                self.db_folder)
        
        per_fit.perf_fitting(self._data['inputs_pm_fit']['machine'],
                             self._data['inputs_pm_fit']['site'])
                             
        pfit = {'c_ext': per_fit.c_ext,
                'k_ext': per_fit.k_ext,
                'c_fit': per_fit.c_fit,
                'k_fit': per_fit.k_fit,
                'c_pto': per_fit.c_pto,
                'k_mooring': per_fit.k_mooring,
                'te': per_fit.te,
                'hm0': per_fit.hm0,
                'wave_dir': per_fit.wave_dir,
                'scatter_diagram': per_fit.scatter_diagram,
                'power_matrix': per_fit.wec_power_matrix,
                'original_power_matrix': per_fit.wec_original_power_matrix,
                'user_power_matrix': self._data['inputs_pm_fit']['machine']['power_matrix']}
                
        self._data['p_fit'] = pfit
        self.trigger_results.emit(pfit)
    
    def skip_fitting(self):
        old_pm = np.copy(self._data['inputs_pm_fit']['machine']['power_matrix'])
        self._data['inputs_pm_fit']['machine']['power_matrix'] *= 0
        
        self.fitting()
        self._data['inputs_pm_fit']['machine']['power_matrix'] = old_pm
        self._data['p_fit']['user_power_matrix'] = old_pm
        self.trigger_results.emit(self._data['p_fit'])
    
    def load_data(self):
        
        ndof = self._data['inputs_hydrodynamic']['general_input']['ndof'][0]
        dat1 = read_file(os.path.join(self.data_folder, 'performance_fit_scatter_diagram.csv'))
        dat2 = read_file(os.path.join(self.data_folder, 'performance_fit_power_matrix.csv'))
        dat3 = read_file(os.path.join(self.data_folder, 'performance_fit_pto.csv'))
        dat4 = read_file(os.path.join(self.data_folder, 'performance_fit_mooring.csv'))
        dat5 = read_file(os.path.join(self.data_folder, 'performance_fit_ext_damping.csv'))
        dat6 = read_file(os.path.join(self.data_folder, 'performance_fit_ext_stiffness.csv'))
        t, h, a, sd = scatter_diagram_data(dat1)
        pm = power_matrix_data(dat2, a, (len(t), len(h), len(a)))
        pto = pto_data(dat3, ndof)
        mooring = mooring_data(dat4, ndof)
        ext_d = ext_damping_data(dat5, ndof)
        ext_k = ext_stiffness_data(dat6, ndof)
                
        machine_spec = {'c_pto': pto,
                        'k_mooring': mooring,
                        'yaw': self.__yaw / 180.0 * np.pi,
                        'power_matrix': pm,
                        'ext_d': ext_d,
                        'ext_k': ext_k}
        
        site_spec = {'spec_shape': 'Jonswap',
                     'spec_gamma': self.__gamma,
                     'spec_spreading': self.__spreading,
                     'te': t,
                     'hm0': h,
                     'wave_angles': a / 180.0 * np.pi,
                     'probability_of_occurence': sd}
        
        stat = check_data_sanity(ndof, machine_spec, site_spec)
        
        if stat:
            print(stat)
            self.btn_fitting.setEnabled(False)
            self.btn_skip_fit.setEnabled(False)
        else:
            self._data['inputs_pm_fit'] = {'machine': machine_spec, 'site': site_spec, 'data_folder': self.data_folder}
            self.btn_skip_fit.setEnabled(True)
            self.btn_plot_pfit.setEnabled(True)
            self.set_dimensions(a)
            if not self._data['inputs_hydrodynamic']['general_input']['get_array_mat']:
                print("WARNING: The hydrodynamic data has been calculated without the foce transfer matrix option")
                print("Therefore the performance fitting cannot be performed.")
                print("To proceed further click on the skip performance fitting button.")
                self.btn_fitting.setEnabled(False)
            else:
                self.btn_fitting.setEnabled(True)
        
            self.trigger_save.emit(self._data['inputs_pm_fit'])

        
    def browse_data(self):
        #if os.path.isdir(self.le_pfit_data.text()):
        
        actual_folder = self.le_pfit_data.text()
            
        folder = QFileDialog.getExistingDirectory(self,
                    "Select the folder that contains the performance fit data",
                    actual_folder)
        self.le_pfit_data.setText(folder)
        stat = check_pfit_folder(folder)
        self.le_pfit_data.setStyleSheet("QLineEdit {border: 0px}")
        if stat:
            print(stat)
            self.le_pfit_data.setStyleSheet("QLineEdit {border: 1px solid red}")
            self.btn_load_pfit.setEnabled(False)
        else:
            self.data_folder = str(folder)
            self.btn_load_pfit.setEnabled(True)
            
        
    def set_data(self, data):
        self._data = data
        
        if 'inputs_pm_fit' in self._data.keys():
            if 'data_folder' in self._data['inputs_pm_fit'].keys():
                self.data_folder = str(self._data['inputs_pm_fit']['data_folder'])
                if not os.path.isdir(self.data_folder):
                    msg = "The data path stored in the project is not valid anymore. The system will attempt to load the data saved in the project. Do you want to load the data (Yes) or specify another data folder (No)"
                    choice = QMessageBox.warning(self, "Warning", msg,QMessageBox.Yes | QMessageBox.No)
                                                 
                    if choice == QMessageBox.Yes:
                        self.le_pfit_data.setText(self.data_folder)
                    else:
                        self.browse_data()
                        

                stat = check_pfit_folder(self.data_folder)
                self.le_pfit_data.setText(self.data_folder)
                if not stat:
                    if 'p_fit' in self._data.keys():
                        self.btn_plot_pfit.setEnabled(True)
                        self.set_dimensions(self._data['inputs_pm_fit']['site']['wave_angles']*180/np.pi)
                        self.trigger_results.emit(self._data['p_fit'])
                    else:
                        self.btn_load_pfit.setEnabled(True)
                        self.load_data()
                        return 0
            
        return -1
                
            


def check_pfit_folder(folder):
    if not os.path.isdir(folder):
        return "Invalid path"
    requirements = ["performance_fit_scatter_diagram.csv", 
                    "performance_fit_power_matrix.csv",
                    "performance_fit_pto.csv",
                    "performance_fit_mooring.csv",
                    "performance_fit_ext_damping.csv",
                    "performance_fit_ext_stiffness.csv"]
    miss_items = []              
    for req in requirements:
        if not req in os.listdir(folder):
            miss_items.append(req)
        
    return miss_items

def read_file(file_name):
    with open(file_name, 'r') as f:
        data = [row for row in csv.reader(f.read().splitlines())]
    return data

def scatter_diagram_data(data):
    te = np.asarray([float(x) for x in data[4] if not x is ''])
    hm0 = np.asarray([float(x) for x in data[6] if not x is ''])
    angles = np.asarray([float(x) for x in data[8] if not x is ''])
    
    if angles.shape[0] < 1:
        return -1
    scatter_diagram = np.zeros((len(te), len(hm0), len(angles)))
    ind = 10
    for angle in angles:
        angle_r = [float(x) for x in data[ind] if x != '' if x != 'Angle'][0]
        a_ind = np.where(angle_r==angles)[0][0]
        for lines in range(len(hm0)):
            scatter_diagram[:,lines,a_ind] = np.asarray([float(x) for x in data[ind+lines+1] if not x is ''])
    
        ind += lines+2
    
    return te, hm0, angles, scatter_diagram
    
def power_matrix_data(data, angles, shape):
    power_matrix = np.zeros(shape)
    ind = 0
    for angle in angles:
        angle_r = [float(x) for x in data[ind] if x != '' if x != 'Angle'][0]
        a_ind = np.where(angle_r==angles)[0][0]
        for lines in range(shape[1]):
            power_matrix[:,lines,a_ind] = np.asarray([float(x) for x in data[ind+lines+1] if not x is ''])
    
        ind += lines+2
    
    return power_matrix

def pto_data(data, ndof):
    pto = []
    for dof in range(ndof):
        pto.append([float(x) for x in data[dof] if not x is ''])
    
    return np.asarray(pto)

def mooring_data(data, ndof):
    mooring = []
    for dof in range(ndof):
        mooring.append([float(x) for x in data[dof] if not x is ''])
    
    return np.asarray(mooring)

def ext_damping_data(data, ndof):
    ext_d = []
    for dof in range(ndof):
        ext_d.append([float(x) for x in data[dof] if not x is ''])
    
    return np.asarray(ext_d)

def ext_stiffness_data(data, ndof):
    ext_k = []
    for dof in range(ndof):
        ext_k.append([float(x) for x in data[dof] if not x is ''])
    
    return np.asarray(ext_k)
    
    
def check_data_sanity(ndof, machine_spec, site_spec):
    error_s = []
    p1 = machine_spec['c_pto']
    p2 = machine_spec['k_mooring']
    p3 = machine_spec['power_matrix']
    p4 = machine_spec['ext_d']
    p5 = machine_spec['ext_k'] 
    
    p6 = site_spec['te']
    p7 = site_spec['hm0']
    p8 = site_spec['wave_angles']
    p9 = site_spec['probability_of_occurence']
    
    if not p1.shape == (ndof, ndof):
        error_s.append("The pto damping shape does not match the number of dof")
        error_s.append("pto damping shape: {} - dofs {}".format(p1.shape, ndof))
    if not p2.shape == (ndof, ndof):
        error_s.append("The mooring stiffness shape does not match the number of dof")
        error_s.append("pto damping shape: {} - dofs {}".format(p2.shape, ndof))
  
    if not p4.shape == (ndof, ndof):
        error_s.append("The external damping shape does not match the number of dof")
        error_s.append("pto damping shape: {} - dofs {}".format(p4.shape, ndof))
    if not p5.shape == (ndof, ndof):
        error_s.append("The external stiffness shape does not match the number of dof")
        error_s.append("pto damping shape: {} - dofs {}".format(p5.shape, ndof))
    
    if not p3.shape == p9.shape: 
        error_s.append("The scatter diagram shape does not match the power matrix one")
        error_s.append("Scatter diagram shape: {} - Power matrix shape {}".format(p9.shape, p3.shape))
    
    if not p9.shape == (len(p6), len(p7), len(p8)):
        error_s.append("The scatter diagram shape does not match the axis")
        error_s.append("Scatter diagram shape: {} - axis shape {}".format(p9.shape,  (len(p6), len(p7), len(p8))))

    elif not np.allclose(p9.sum(), 1., atol=1e-1):
        error_s.append("The probability of occurence of the different sea state does not sum to 1.")
        error_s.append("Verify the validity of the input data")
    return error_s
