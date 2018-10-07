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
import shutil

import file_utilities as f_util


def write_to_raw_results(project_folder):
    temp_folder = os.path.join(project_folder, 'raw_results')
    if os.path.isdir(temp_folder):
        if os.listdir(temp_folder):
            f_util.clean_prj_folder(temp_folder)
    else:
        os.mkdir(temp_folder)
    
    for el in os.listdir(project_folder):
        el_abs = os.path.join(project_folder, el)
        if os.path.isfile(el_abs):
            shutil.copy(el_abs, temp_folder)
    f_util.copy_result_to_project(project_folder, temp_folder)
    f_util.clean_prj_folder(project_folder, exept='raw_results')
    
#    os.mkdir(os.path.join(project_folder, 'hydrodynamic'))
#    os.mkdir(os.path.join(project_folder, 'hydrostatic'))
#    os.mkdir(os.path.join(project_folder, 'performance_fit'))
#    os.mkdir(os.path.join(project_folder, 'hydrodynamic','results'))
#    os.mkdir(os.path.join(project_folder, 'hydrodynamic','mesh'))
#    os.mkdir(os.path.join(project_folder, 'hydrostatic','body0'))
#    os.mkdir(os.path.join(project_folder, 'hydrostatic','body0','mesh'))
    
def write_cm(project_folder, m_add, periods):
    (nper, ndof, ndofj) = m_add.shape
    fid = open(os.path.join(project_folder,'hydrodynamic','results','CM.dat'),'w')
    fid.write('Number of periods:     {}\n'.format(int(nper)))
    for i_p, per in enumerate(periods): 
        fid.write(' {:.4f}\n'.format(per))
        for dof in range(ndof):
            line_range = 6
            for dofj in range(ndof):
                fid.write('  {:.6E}\t'.format(m_add[i_p, dof, dofj]))
                if dofj >= line_range-1:
                    fid.write('\n')
                    line_range += 6   
    fid.close()
    
def write_ca(project_folder, c_rad, periods):
    (nper, ndof, ndofj) = c_rad.shape
    fid = open(os.path.join(project_folder,'hydrodynamic','results','CA.dat'),'w')
    fid.write('Number of periods:     {}\n'.format(int(nper)))
    for i_p, per in enumerate(periods): 
        fid.write(' {:.4f}\n'.format(per))
        for dof in range(ndof):
            line_range = 6
            for dofj in range(ndof):
                fid.write('  {:.6E}\t'.format(c_rad[i_p, dof, dofj]))
                if dofj >= line_range-1:
                    fid.write('\n')
                    line_range += 6   
    fid.close()

def write_cyl_surface(project_folder, phi_d, phi_r):
    pass
def write_f_ex(project_folder, f_ex, periods, directions):
    (nper, ndir, ndof) = f_ex.shape
    fid = open(os.path.join(project_folder,'hydrodynamic','results','ExcitationForce.tec'),'w')
    fid.write('VARIABLES="w (rad/s)"\n')
    for dof in range(ndof):
        fid.write("abs(F   1   {})" "angle(F   1   {})\n".format(dof))
    
    fid.close()
def write_hst(project_folder, k_hst):
    pass