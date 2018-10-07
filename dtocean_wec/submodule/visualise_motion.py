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
Created on Tue May 31 09:08:16 2016

.. moduleauthor:: Francesco Ferri <ff@civil.aau.dk>
.. moduleauthor:: Mathew Topper <mathew.topper@dataonlygreater.com>
"""

import numpy as np

from dtocean_wave.utils.StrDyn import MotionFreq

from matplotlib.figure import Figure
from matplotlib.backends.backend_qt4agg import (
                                        FigureCanvasQTAgg as FigureCanvas)


class Visualiser():
    def __init__(self):
        """Place Holder, for a future development
        """
        self.__rao = None
        self._data_h = None
        self._data_p = None
    
    def set_hydrodynamic_data(self, data):
        status = self._check_data_hydro(data)
        if status[0]:
            self._data_h = data
            return (True, )
        else:
            return status
        
    def set_performance_fit_data(self, data):
        status = self._check_data_pfit(data)
        if status[0]:
            self._data_p = data
            return (True, )
        else:
            return status
        
    def _check_data_hydro(self, data):
        # add a check on the dictionaries keys and dimensions
        return (True, )
    
    def _check_data_pfit(self, data):
        # add a check on the dictionaries keys and dimensions
        return (True, )
        
    def show_k_fit(self, te_i, hs_i, ang_i):
        if self._data_p is None:
            print('Missing data to show the fitted stiffness')
            return Figure()
            
        fig = Figure()
        ax = fig.add_subplot(111)
        sz_ticks=15
        sz_labels=25
        dat = self._data_p['k_fit'][te_i, hs_i, ang_i,:,:]
        dofi, dofj = dat.shape
        
        mtsh = ax.matshow(dat.T, aspect="auto")
        for (i, j), z in np.ndenumerate(dat.T):
            ax.text(j, i, '{:.0f}'.format(z), ha='center', va='center')
        ax.set_xticks(range(dofi))
        ax.set_xticklabels(['{:.1f}'.format(tp) for tp in range(dofi)], fontsize = sz_ticks)
        ax.set_xlabel('dof i ', fontsize = sz_labels)
        ax.xaxis.set_label_position('top')
        ax.set_yticks(range(dofj))
        ax.set_yticklabels(['{:.1f}'.format(hs) for hs in range(dofj)], fontsize = sz_ticks)
        ax.set_ylabel('dof j', fontsize = sz_labels)
        cb = Figure.colorbar(fig, mtsh)
        cb.set_label('Fitted Stiffness (N/m [Nm/rad]) ', fontsize = sz_labels)
        ax.grid(False)
        return fig
    
    def show_c_fit(self, te_i, hs_i, ang_i):
        if self._data_p is None:
            print('Missing data to show the fitted stiffness')
            return Figure()
            
        fig = Figure()
        ax = fig.add_subplot(111)
        sz_ticks=15
        sz_labels=25
        dat = self._data_p['c_fit'][te_i, hs_i, ang_i,:,:]
        dofi, dofj = dat.shape
        
        mtsh = ax.matshow(dat.T, aspect="auto")
        for (i, j), z in np.ndenumerate(dat.T):
            ax.text(j, i, '{:.0f}'.format(z), ha='center', va='center')
        ax.set_xticks(range(dofi))
        ax.set_xticklabels(['{:.1f}'.format(tp) for tp in range(dofi)], fontsize = sz_ticks)
        ax.set_xlabel('dof i ', fontsize = sz_labels)
        ax.xaxis.set_label_position('top')
        ax.set_yticks(range(dofj))
        ax.set_yticklabels(['{:.1f}'.format(hs) for hs in range(dofj)], fontsize = sz_ticks)
        ax.set_ylabel('dof j', fontsize = sz_labels)
        cb = Figure.colorbar(fig, mtsh)
        cb.set_label('Fitted Damping (N/(m/s) [Nm/(rad/s)]) ', fontsize = sz_labels)
        ax.grid(False)
        return fig
        
    def show_mass(self):
        if self._data_h is None:
            print('Missing data to show the mass matrix')
            return Figure()
            
        fig = Figure()
        ax = fig.add_subplot(111)
        sz_ticks=15
        sz_labels=25
        mass = self._data_h['m_m']
        dofi, dofj = mass.shape
        
        mtsh = ax.matshow(mass.T, aspect="auto")
        for (i, j), z in np.ndenumerate(mass.T):
            ax.text(j, i, '{:.0f}'.format(z), ha='center', va='center')
        ax.set_xticks(range(dofi))
        ax.set_xticklabels(['{:.1f}'.format(tp) for tp in range(dofi)], fontsize = sz_ticks)
        ax.set_xlabel('dof i ', fontsize = sz_labels)
        ax.xaxis.set_label_position('top')
        ax.set_yticks(range(dofj))
        ax.set_yticklabels(['{:.1f}'.format(hs) for hs in range(dofj)], fontsize = sz_ticks)
        ax.set_ylabel('dof j', fontsize = sz_labels)
        cb = Figure.colorbar(fig, mtsh)
        cb.set_label('Mass (kg) ', fontsize = sz_labels)
        ax.grid(False)
        return fig
            
    def show_hst(self):
        if self._data_h is None:
            print('Missing data to show the hydrostatic matrix')
            return Figure()
        
        fig = Figure()
        ax = fig.add_subplot(111)
        sz_ticks=15
        sz_labels=25
        hst = self._data_h['k_hst']
        dofi, dofj = hst.shape
        
        mtsh = ax.matshow(hst.T, aspect="auto")
        for (i, j), z in np.ndenumerate(hst.T):
            ax.text(j, i, '{:.0f}'.format(z), ha='center', va='center')
        ax.set_xticks(range(dofi))
        ax.set_xticklabels(['{:.1f}'.format(tp) for tp in range(dofi)], fontsize = sz_ticks)
        ax.set_xlabel('dof i ', fontsize = sz_labels)
        ax.xaxis.set_label_position('top')
        ax.set_yticks(range(dofj))
        ax.set_yticklabels(['{:.1f}'.format(hs) for hs in range(dofj)], fontsize = sz_ticks)
        ax.set_ylabel('dof j', fontsize = sz_labels)
        cb = Figure.colorbar(fig, mtsh)
        cb.set_label('Hydrostatic Stiffness (N/m [Nm/rad])', fontsize = sz_labels)
        ax.grid(False)
        return fig
        
    def show_radiation_problem(self, d_i, d_j):
        if self._data_h is None:
            print('Missing data to show the radiation coefficients')
            return Figure()
		
        d_i = int(d_i)
        d_j = int(d_j)
		
        fig = Figure()
        host = fig.add_subplot(111)
		
        par1 = host.twinx()

        p1, = host.plot(1/self._data_h['periods'],
                        self._data_h['m_add'][:,int(d_i), int(d_j)],
                        "b-o", label="added mass")
        p2, = par1.plot(1/self._data_h['periods'],
                        self._data_h['c_rad'][:,int(d_i), int(d_j)],
                        "g-o", label="radiation damping")


        host.set_xlabel("frequency, [Hz]")
        host.set_ylabel("added mass, [(N/(m/s^2))/(Nm/(rad/s^2))]")
        par1.set_ylabel("radiation damping, [(N/(m/s))/(Nm/(rad/s))]")

        host.yaxis.label.set_color(p1.get_color())
        par1.yaxis.label.set_color(p2.get_color())

        tkw = dict(size=4, width=1.5)
        host.tick_params(axis='y', colors=p1.get_color(), **tkw)
        par1.tick_params(axis='y', colors=p2.get_color(), **tkw)
        host.tick_params(axis='x', **tkw)

        lines = [p1, p2]

        host.legend(lines, [l.get_label() for l in lines])

        return fig

    def show_diffraction_problem(self, d_i, dir_i):
        if self._data_h is None:
            print('Missing data to show the diffraction coefficients')
            return Figure()
			
        fig = Figure()
        host = fig.add_subplot(111)
        
        par1 = host.twinx()
        p1, = host.plot(1/self._data_h['periods'],
                        np.abs(self._data_h['f_ex'][:,int(dir_i), int(d_i)]),
                        "k-o", label="manitude")
        
        p2, = par1.plot(1/self._data_h['periods'],
                       np.angle(self._data_h['f_ex'][:,int(dir_i), int(d_i)]),
                       "r-o", label="phase")
        
        host.set_xlabel("frequency, [Hz]")
        host.set_ylabel("magnitude, [(N/m)/(Nm/m)]")
        par1.set_ylabel("phase, [rad/m]")

        host.yaxis.label.set_color(p1.get_color())
        par1.yaxis.label.set_color(p2.get_color())

        tkw = dict(size=4, width=1.5)
        host.tick_params(axis='y', colors=p1.get_color(), **tkw)
        par1.tick_params(axis='y', colors=p2.get_color(), **tkw)
        host.tick_params(axis='x', **tkw)

        lines = [p1, p2]

        host.legend(lines, [l.get_label() for l in lines])

        return fig        
        
    def show_rao(self, te_i, hs_i, ang_i, d_i, dir_i):
        if self._data_h is None or self._data_p is None:
            print('Missing data to show the rao')
            return Figure()
        rao = self.__eval_rao(int(te_i), int(hs_i), int(ang_i))
        fig = Figure()
        host = fig.add_subplot(111)

        par1 = host.twinx()

        p1, = host.plot(1/self._data_h['periods'],
                        np.abs(rao[:,int(dir_i), int(d_i)]),
                        "k-o", label="manitude")
        p2, = par1.plot(1/self._data_h['periods'],
                        np.angle(rao[:,int(dir_i), int(d_i)]),
                        "r-o", label="phase")

        host.set_xlabel("frequency, [Hz]")
        host.set_ylabel("rao manitude, [(m/m),(rad/m)]")
        par1.set_ylabel("rao phase, [rad/m]")

        host.yaxis.label.set_color(p1.get_color())
        par1.yaxis.label.set_color(p2.get_color())

        tkw = dict(size=4, width=1.5)
        host.tick_params(axis='y', colors=p1.get_color(), **tkw)
        par1.tick_params(axis='y', colors=p2.get_color(), **tkw)
        host.tick_params(axis='x', **tkw)

        lines = [p1, p2]

        host.legend(lines, [l.get_label() for l in lines])

        return fig
    
    def show_original_power_matrix(self, ang_i):
        if self._data_h is None or self._data_p is None:
            print('Missing data to show the power matrix')
            return Figure()
            
        fig = Figure()
        ax = fig.add_subplot(111)
        sz_ticks=15
        sz_labels=25
        power_matrix = self._data_p['original_power_matrix'][:,:,ang_i]
        dofi, dofj = power_matrix.shape
        
        mtsh = ax.matshow(power_matrix.T, aspect="auto")
        for (i, j), z in np.ndenumerate(power_matrix.T):
            ax.text(j, i, '{:.0f}'.format(z), ha='center', va='center')
        ax.set_xticks(range(dofi))
        ax.set_xticklabels(['{:.1f}'.format(te) for te in self._data_p['te']], fontsize = sz_ticks)
        ax.set_xlabel('Te ', fontsize = sz_labels)
        ax.xaxis.set_label_position('top')
        ax.set_yticks(range(dofj))
        ax.set_yticklabels(['{:.1f}'.format(hs) for hs in self._data_p['hm0']], fontsize = sz_ticks)
        ax.set_ylabel('Hm0', fontsize = sz_labels)
        cb = Figure.colorbar(fig, mtsh)
        cb.set_label('Power (W) - wave dir {}deg'.format(self._data_p['wave_dir'][ang_i]*180/np.pi), fontsize = sz_labels)
        ax.grid(False)
        return fig
    
    def show_user_power_matrix(self, ang_i):
        if self._data_h is None or self._data_p is None:
            print('Missing data to show the power matrix')
            return Figure()
            
        fig = Figure()
        ax = fig.add_subplot(111)
        sz_ticks=15
        sz_labels=25
        power_matrix = self._data_p['user_power_matrix'][:,:,ang_i]
        dofi, dofj = power_matrix.shape
        
        mtsh = ax.matshow(power_matrix.T, aspect="auto")
        for (i, j), z in np.ndenumerate(power_matrix.T):
            ax.text(j, i, '{:.0f}'.format(z), ha='center', va='center')
        ax.set_xticks(range(dofi))
        ax.set_xticklabels(['{:.1f}'.format(te) for te in self._data_p['te']], fontsize = sz_ticks)
        ax.set_xlabel('Te ', fontsize = sz_labels)
        ax.xaxis.set_label_position('top')
        ax.set_yticks(range(dofj))
        ax.set_yticklabels(['{:.1f}'.format(hs) for hs in self._data_p['hm0']], fontsize = sz_ticks)
        ax.set_ylabel('Hm0', fontsize = sz_labels)
        cb = Figure.colorbar(fig, mtsh)
        cb.set_label('Power (W) - wave dir {}deg'.format(self._data_p['wave_dir'][ang_i]*180/np.pi), fontsize = sz_labels)
        ax.grid(False)
        return fig
        
    def show_power_matrix(self, ang_i):
        if self._data_h is None or self._data_p is None:
            print('Missing data to show the power matrix')
            return Figure()
            
        fig = Figure()
        ax = fig.add_subplot(111)
        sz_ticks=15
        sz_labels=25
        power_matrix = self._data_p['power_matrix'][:,:,ang_i]
        dofi, dofj = power_matrix.shape
        
        mtsh = ax.matshow(power_matrix.T, aspect="auto")
        for (i, j), z in np.ndenumerate(power_matrix.T):
            ax.text(j, i, '{:.0f}'.format(z), ha='center', va='center')
        ax.set_xticks(range(dofi))
        ax.set_xticklabels(['{:.1f}'.format(te) for te in self._data_p['te']], fontsize = sz_ticks)
        ax.set_xlabel('Te ', fontsize = sz_labels)
        ax.xaxis.set_label_position('top')
        ax.set_yticks(range(dofj))
        ax.set_yticklabels(['{:.1f}'.format(hs) for hs in self._data_p['hm0']], fontsize = sz_ticks)
        ax.set_ylabel('Hm0', fontsize = sz_labels)
        cb = Figure.colorbar(fig, mtsh)
        cb.set_label('Power (W) - wave dir {}deg'.format(self._data_p['wave_dir'][ang_i]*180/np.pi), fontsize = sz_labels)
        ax.grid(False)
        return fig

    def show_power_function(self):
        pass

    def show_mesh(self):
        pass

    def __eval_rao(self, te, hs, ang):
        print("Evaluating the rao")
        omeg = 2.*np.pi/self._data_h['periods']
        nte, nhs, nang, ndof, ndofj = self._data_p['c_pto'].shape
        nfr = self._data_h['m_add'].shape[0]
        ndir = self._data_h['f_ex'].shape[1]

        rao = np.zeros((nfr, ndir, ndof),dtype='complex')
        for fr in range(nfr):
            for dr in range(ndir):
                rao[fr, dr, :] = MotionFreq(self._data_h['m_m'] + self._data_h['m_add'][fr,:,:],
                                      self._data_h['c_rad'][fr, :, :] + self._data_p['c_pto'][te, hs, ang, :, :] + self._data_p['c_fit'][te, hs, ang, :, :] + self._data_p['c_ext'][te, hs, ang, :, :],
                                      self._data_h['k_hst'] + self._data_p['k_mooring'][te, hs, ang, :, :] + self._data_p['k_fit'][te, hs, ang, :, :] + self._data_p['k_ext'][te, hs, ang, :, :],
                                      self._data_h['f_ex'][fr, dr].T,
                                      omeg[fr])

        return rao

if __name__ == "__main__":
    vi = Visualiser()
    data = {'periods': np.array([0,1]), 'f_ex': np.array([[[0,1]]], dtype=complex).T,
            'c_rad': np.array([[[0,1]]], 'f').T, 'm_add': np.array([[[0,1]]], 'f').T}
    figu = vi.show_diffraction_problem(0, 0)
    vi.set_hydrodynamic_data(data)
    figu = vi.show_radiation_problem(0,0)
    canvas = FigureCanvas(figu)
    canvas.draw()
    canvas.show()
    