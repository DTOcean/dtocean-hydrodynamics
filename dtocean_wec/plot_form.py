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
Created on Wed Jun 15 09:19:49 2016

.. moduleauthor:: Francesco Ferri <ff@civil.aau.dk>
"""

from PyQt4.QtCore import *
from PyQt4.QtGui import *

from plotter import Ui_Form as Ui_Plot
from submodule.visualise_motion import Visualiser

from submodule.utils.conversion_utils import *

from matplotlib.figure import Figure
from matplotlib.backends.backend_qt4agg import (
    FigureCanvasQTAgg as FigureCanvas,
    NavigationToolbar2QT as NavigationToolbar)

class Plotter(QWidget, Ui_Plot):
    trigger = pyqtSignal()
    def __init__(self, parent=None):
        QWidget.__init__(self, parent)
        self.setupUi(self)

        self.btn_plot_update.clicked.connect(self.plot_data)
        self.plotter = Visualiser()
        self.start_figure()
        
        QToolTip.setFont(QFont("SansSerif", 11))
        self.label_51.setToolTip("Wave direction expressed\nin the mesh coordinate system.\n0deg correspond to a wave travelling in the positive x-direction")
                                
        self.label_54.setToolTip("Wave direction expressed\nin the North-South/West-East coordinate system.\n0deg correspond to a wave travelling from North to South")
        
    def set_data(self, data):
        self._data = data
        if 'hyd' in data.keys() and not 'p_fit' in data.keys():
            self.enable_plot_area([1])
        elif 'p_fit' in data.keys():
            self.enable_plot_area([1,2])
        else:
            self.enable_plot_area()
            
    def start_figure(self):
        fig = Figure()
        self.addmpl(fig)
        self.enable_plot_area()
        
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
    
    def plot_excitation(self):
        fig = self.plotter.show_diffraction_problem(int(self.cb_dofi.currentIndex()),
                                                    float(self.cb_angle.currentIndex()))
        self.rmmpl()
        self.addmpl(fig)
# 
    def plot_radiation(self):
        fig = self.plotter.show_radiation_problem(int(self.cb_dofi.currentIndex()),
                                                  float(self.cb_dofj.currentIndex()))
        self.rmmpl()
        self.addmpl(fig)
    
    def plot_rao(self):
        fig = self.plotter.show_rao(int(self.cb_te.currentIndex()),
                                    int(self.cb_hm0.currentIndex()), 
                                    int(self.cb_wavedir.currentIndex()),
                                    int(self.cb_dofi.currentIndex()),
                                    int(self.cb_angle.currentIndex()))
        self.rmmpl()
        self.addmpl(fig)
    
    def plot_power_matrix(self):
        fig = self.plotter.show_power_matrix(int(self.cb_wavedir.currentIndex()))
        self.rmmpl()
        self.addmpl(fig)
    
    def plot_k_fit(self):
        fig = self.plotter.show_k_fit(int(self.cb_te.currentIndex()),
                                    int(self.cb_hm0.currentIndex()), 
                                    int(self.cb_wavedir.currentIndex()))
        self.rmmpl()
        self.addmpl(fig)
        
    def plot_c_fit(self):
        fig = self.plotter.show_c_fit(int(self.cb_te.currentIndex()),
                                    int(self.cb_hm0.currentIndex()), 
                                    int(self.cb_wavedir.currentIndex()))
        self.rmmpl()
        self.addmpl(fig)
    
    def plot_original_power_matrix(self):
        fig = self.plotter.show_original_power_matrix(int(self.cb_wavedir.currentIndex()))
        self.rmmpl()
        self.addmpl(fig)
    
    def plot_user_power_matrix(self):
        fig = self.plotter.show_user_power_matrix(int(self.cb_wavedir.currentIndex()))
        self.rmmpl()
        self.addmpl(fig)
        
    
    def plot_mass(self):
        fig = self.plotter.show_mass()
        self.rmmpl()
        self.addmpl(fig)
    
    def plot_hydrostatic(self):
        fig = self.plotter.show_hst()
        self.rmmpl()
        self.addmpl(fig)
        
    def enable_plot_area(self, index=[]):
        self.plotter.set_hydrodynamic_data(None)
        self.plotter.set_performance_fit_data(None)
        if index:
            self.groupBox_11.setEnabled(True)
            self.groupBox_9.setEnabled(True)
            self.btn_plot_update.setEnabled(True)
            for ind in index:
                if ind == 1:
                    self.enable_hydrodynamic(True)
                    self.set_hydrodynamic_dimensions(True)
                    self.plotter.set_hydrodynamic_data(self._data['hyd'])
                elif ind == 2:
                    self.enable_motion(True)
                    self.set_motion_dimensions(True)
                    self.plotter.set_performance_fit_data(self._data['p_fit'])
                else:
                    pass
        else:
            self.groupBox_11.setEnabled(False)
            self.groupBox_9.setEnabled(False)
            self.btn_plot_update.setEnabled(False)
            self.enable_hydrodynamic(False)
            self.enable_motion(False)
            

    def enable_hydrodynamic(self, choice):
        self.rb_excitation.setEnabled(choice)
        self.rb_radiation.setEnabled(choice)
        self.rb_mass.setEnabled(choice)
        self.rb_stiffness.setEnabled(choice)
        self.rb_radiation.setChecked(True)
        
    def enable_motion(self, choice):
        self.c_fit.setEnabled(choice)
        self.k_fit.setEnabled(choice)
        self.rb_userpowermat.setEnabled(choice)
        self.rb_origpowermat.setEnabled(choice)
        self.rb_rao.setEnabled(choice)
        self.rb_powermat.setEnabled(choice)
        
    def set_motion_dimensions(self, stat):
        if stat:
            te = self._data['p_fit']['te'].tolist()
            te = [str(x) for x in te]
            hm0 = self._data['p_fit']['hm0'].tolist()
            hm0 = [str(x) for x in hm0]
            wave_dir = self._data['p_fit']['wave_dir']
            
            wave_dir_ne = angle_wrap(-wave_dir-np.pi/2,'r2r')*180/np.pi
            wave_dir = [str(x) for x in wave_dir_ne.tolist()]

            
            self.cb_te.clear()
            self.cb_te.addItems(te)
            self.cb_hm0.clear()
            self.cb_hm0.addItems(hm0)
            self.cb_wavedir.clear()
            self.cb_wavedir.addItems(wave_dir)
        else:
            self.cb_te.clear()
            self.cb_hm0.clear()
            self.cb_wavedir.clear()
            
    def set_hydrodynamic_dimensions(self, stat):
        if stat:
            angle = self._data['hyd']['directions'].tolist()
            angle = [str(x) for x in angle]
            dofi = range(self._data['hyd']['m_m'].shape[0])
            dofi = [str(x) for x in dofi]
            self.cb_angle.clear()
            self.cb_angle.addItems(angle)
            self.cb_dofi.clear()
            self.cb_dofi.addItems(dofi)
            self.cb_dofj.clear()
            self.cb_dofj.addItems(dofi)
        else:
            self.cb_angle.clear()
            self.cb_dofi.clear()
            self.cb_dofj.clear()
        
    def plot_data(self):
        # check which tab is active first
        active_rb = [ix for ix, x in enumerate(self.groupBox_9.children()) if ix>0 and x.isChecked()]
        
        if active_rb[0] == 1:
            self.plot_radiation()
        elif active_rb[0] == 2:
            self.plot_excitation()
        elif active_rb[0] == 3:
            self.plot_mass()
        elif active_rb[0] == 4:
            self.plot_hydrostatic()
        elif active_rb[0] == 5:
            self.plot_c_fit()
        elif active_rb[0] == 6:
            self.plot_k_fit()
        elif active_rb[0] == 7:
            self.plot_original_power_matrix()
        elif active_rb[0] == 8:
            self.plot_user_power_matrix()
        elif active_rb[0] == 9:
            self.plot_power_matrix()
        elif active_rb[0] == 10:
            self.plot_rao()
        elif active_rb[0] == 11:
            pass
