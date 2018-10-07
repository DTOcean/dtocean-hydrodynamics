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
Created on Tue May 31 09:08:16 2016

.. moduleauthor:: Francesco Ferri <ff@civil.aau.dk>
"""

import numpy as np

from matplotlib.figure import Figure
from matplotlib.backends.backend_qt4agg import (
    FigureCanvasQTAgg as FigureCanvas,
    NavigationToolbar2QT as NavigationToolbar)


class VisualisePerformanceFit():
    def __init__(self):
        """Place Holder, for a future development
        """
        self._data_pf = None
    
    def set_data(self, data):
        self._data_pf = data
        
    def show_pto(self):
        if self._data_pf is None:
            print('Missing data to assess the radiation problem')
            return Figure()
            
        
        val = self._data_pf['machine']['c_pto']
        dofi = range(val.shape[0])
        
        return self.general_plot(val, dofi, dofi, 'dof i', 'dof j', 'Damping, [-]')

    def show_mooring(self):
        if self._data_pf is None:
            print('Missing data to assess the radiation problem')
            return Figure()
            
        
        val = self._data_pf['machine']['k_mooring']
        dofi = range(val.shape[0])
        
        return self.general_plot(val, dofi, dofi, 'dof i', 'dof j', 'Stiffness, [-]')
    
    def show_ext_d(self):
        if self._data_pf is None:
            print('Missing data to assess the radiation problem')
            return Figure()
            
        
        val = self._data_pf['machine']['ext_d']
        dofi = range(val.shape[0])
        
        return self.general_plot(val, dofi, dofi, 'dof i', 'dof j', 'Damping, [-]')
    
    def show_ext_k(self):
        if self._data_pf is None:
            print('Missing data to assess the radiation problem')
            return Figure()
            
        
        val = self._data_pf['machine']['ext_k']
        dofi = range(val.shape[0])
        
        return self.general_plot(val, dofi, dofi, 'dof i', 'dof j', 'Stiffness, [-]')
        
        
    def show_scatter_diagram(self, a_i):
        if self._data_pf is None:
            print('Missing data to assess the radiation problem')
            return Figure()
            
        
        val = self._data_pf['site']['probability_of_occurence'][:,:,a_i]
        te =  self._data_pf['site']['te']
        hm0 = self._data_pf['site']['hm0']
        
        return self.general_plot(val, te, hm0, 'Te, [s]', 'Hm0, [m]', 'Probability, [-]')
        
        
    def show_power_matrix(self, a_i):
        if self._data_pf is None:
            print('Missing data to assess the radiation problem')
            return Figure()
            
        
        val = self._data_pf['machine']['power_matrix'][:,:,a_i]
        te =  self._data_pf['site']['te']
        hm0 = self._data_pf['site']['hm0']
        
        return self.general_plot(val, te, hm0, 'Te, [s]', 'Hm0, [m]', 'Power, [W]')
        
    def general_plot(self, val, x, y, x_t, y_t, z_t):
        fig = Figure()
        ax = fig.add_subplot(111)
        sz_ticks=15
        sz_labels=25
        mtsh = ax.matshow(val.T, aspect="auto")
        for (i, j), z in np.ndenumerate(val.T):
            ax.text(j, i, '{:.0f}'.format(z), ha='center', va='center')
        ax.set_xticks(range(len(x)))
        ax.set_xticklabels(['{:.1f}'.format(tp) for tp in x], fontsize = sz_ticks)
        ax.set_xlabel(x_t, fontsize = sz_labels)
        ax.xaxis.set_label_position('top')
        ax.set_yticks(range(len(y)))
        ax.set_yticklabels(['{:.1f}'.format(hs) for hs in y], fontsize = sz_ticks)
        ax.set_ylabel(y_t, fontsize = sz_labels)
        cb = Figure.colorbar(fig, mtsh)
        cb.set_label(z_t, fontsize = sz_labels)
        ax.grid(False)
        return fig
            
   
