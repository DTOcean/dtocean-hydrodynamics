# -*- coding: utf-8 -*-

#    Copyright (C) 2016 Pau Mercadez Ruiz
#    Copyright (C) 2017-2022 Mathew Topper
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
This module contains the class used to generate the numerical model of the WEC
using the Nemoh software.

.. module:: hyd_Nemoh
   :platform: Windows
   :synopsis: Numerical model of WEC builder

.. moduleauthor:: Pau Mercadez Ruiz <pmr@civil.aau.dk>
.. moduleauthor:: Mathew Topper <mathew.topper@dataonlygreater.com>
"""

import os
import logging
import subprocess

import numpy as np
from numpy import linalg as LA

from hydrostatics import Hydrostatics_Nemohcal
from utils.mesh import MeshBem
from utils.multibody_analysis import MultiBodyAnalysis

# Start logging
module_logger = logging.getLogger(__name__)


class NemohExecute():
    """
    BEM_Nemoh: the class is used to generate the necessary files to run the BEM solver Nemoh,
                and to collect the results into a frequency domain model of the machine

    Args:
        iHydro (Hydro_pkg class): the object collects the information of the sea states, frequencies and directions
                                    to be analysed
        iInput (WP2input class): the object collects the WP2 package inputs

    Optional args:
        c_Nz (int): number of points used to discretise the cylinder circumscribing the WEC in the vertical direction
        c_Nth (int): number of points used to discretise the cylinder circumscribing the WEC in the angular direction
        debug (boolean): debug flag

    Attributes:
        _debug (boolean): debug flag
        water_depth (float): average water depth at the site
        Cyl_Nzeta (int):  number of points used to discretise the cylinder circumscribing the WEC in the vertical direction
        Cyl_Ntheta (int): number of points used to discretise the cylinder circumscribing the WEC in the angular direction
        reRunNemoh (boolean): flag used to decide whether to re-run Nemoh in case a previous calculation is available
        lCS (numpy.ndarray): position of the local coordinate system in the x, y, z space. The space is defined by the mesh file
        frequency_def (numpy.ndarray): min, max and number of wave frequency [rad/s] to be analysed
        angle_def  (numpy.ndarray): min, max and number of wave angles [rad] to be analysed
        in_path (str): path name where the WECmodes.wp2 is stored
        in_f_nm (str): file name of the WECmodes.wp2 file
        HOME (str): path name of the current folder
        path_hyd (str): folder name where the hydrodynamic solver will store the results
        totDOF (float): total number of degree of freedom of the system
        ptoDOF (numpy.ndarray): ID of the dofs connected to the power take off system
        mooDOF (numpy.ndarray): ID of the dofs connected to the mooring system
        MassM (numpy.ndarray): mass matrix for the machine, the dimension needs to agree with the number of totDOF
        sharedDOF (list): shared dofs given with respect to the mesh 0,0,0 point
        nDOF_shared (int): number of shared dof, these dofs describes the normal dof of the WEC
        extendedDOF (list): extended dofs given with respect to the mesh 0,0,0 point
        nDOF_extended (int): number of extended dof, used to build the generalised dof space for the multibody WEC
        n_bodies (int): number of body composing the WEC
        LocalAxes (list): local coordinate system position for each body composing the WEC.
                            The local coordiante system position is given with respect to the mesh 0,0,0 point.
        header (list): header string for each body composing the WEC, used in the generation of the nemoh.cal
        mesh_file (list): mesh name for each body composing the WEC
        nDOF_b (list): number of dofs for each body composing the WEC
        connectivity (list): definition of the connectivity between bodies for a multibody WEC
        nV (list): list containing the number of vertex for each body mesh
        nP (list): list ontaining the number of panels for each body mesh
        hyd_dir_name (str): name of the folder where the calculation is stored
        hydro_dir_path (str): defining the location of the hydrodynamic folder containing the BEM solver solution
        hydro_dir_path_res (str): defining the location of the result folder
        hydro_dir_path_mesh (str): path name defining the location of the mesh files to be used in the BEM software
        hydro_dir_path_hst (str): path name defining the location of the hydrostatic calculation related files
        hydro_dir_path_hst_bID (list): list of strings with the path name of the #ID body location used for the hydrostatic
                                        calculation
        hydro_dir_path_hst_bID_mesh (list): list of strings with the path name of the #ID body mesh used  for the hydrostatic
                                        calculation
        Outputs (tuple): contains the output of the class:
                           Scattering,
                           Radiation_m,
                           Fex_m,
                           Cm_m,
                           Ca_m,
                           Khyd,
                           self.MassM,
                           position,
                           self.Cyl_radius,
                           Th_cyl,
                           Z_cyl,
                           water_depth,
                           directions,
                           period,
                           dof

    Returns:
        Outputs (tuple): contains the output of the class:
                           Scattering,
                           Radiation_m,
                           Fex_m,
                           Cm_m,
                           Ca_m,
                           Khyd,
                           self.MassM,
                           position,
                           self.Cyl_radius,
                           Th_cyl,
                           Z_cyl,
                           water_depth,
                           directions,
                           period,
                           dof

    Note: the overall procedure to generate a numerical model if simplified to 5 steps:
        1) read inputs
        2) generate Nemoh folder structure
        3) generate Nemoh files
        4) run Nemoh
            a) hydrostatic
            b) hydrodynamic
        5) read results and generate the numerical model

    """
    def __init__(self,
                 prj_folder,
                 general_inputs,
                 body_inputs, debug=False,
                 get_array_mat=True): #wdepth, l_cs, c_Nz=50, c_Nth=60,
                
        self.rho = 1025
        self.g = 9.82
        self._debug = debug
        self.prj_folder = prj_folder
        self.__get_array_mat = get_array_mat
        self.g_d = general_inputs
        self.b_d = body_inputs

        self.water_depth = self.g_d['water_depth']
        self.Cyl_Nzeta = self.g_d['cyl_nzeta']
        self.Cyl_Ntheta = self.g_d['cyl_ntheta']
        self.tot_dof = self.g_d['ndof']
        self.pto_dof = self.g_d['pto_dof']
        self.moor_dof = self.g_d['mooring_dof']
        self.point_application = self.b_d['local_cs']
        self.frequency_def = self.g_d['frequency_def']
        self.angle_def = self.g_d['angle_def']
        
        self.shared_dof_binary = self.b_d['shared_dof']
        
        bodies = self.b_d['body'].values()  # list of bodies
        self.bodies = bodies
        self.n_bodies = len(bodies)
        
        
        self.vel_map = None
        self.bodies_dofs = None
        
        self.mass_m = None
        self.input_f_n = None
        self.LocalAxes = None
        self.header = None
        self.mesh_file = None
        self.nDOF_b = None
        self.connectivity = None
        self.sharedDOF = None
        self.nDOF_shared = None
        self.extendedDOF = None
        self.nDOF_extended = None
        self.Cyl_radius = None

        self.mesh_obj = []
        for nb in range(self.n_bodies):
            path = self.bodies[nb]['mesh']
            msh_obj = MeshBem( os.path.split(path)[1],  os.path.split(path)[0])

            self.bodies[nb]['nV'] = msh_obj.nV
            self.bodies[nb]['nP'] = msh_obj.nP
            self.mesh_obj += [msh_obj]

        # identify the size of the circumscribing cylinder
        if (self.Cyl_Nzeta == 0 or
                self.Cyl_Ntheta == 0 or
                    not self.__get_array_mat):
            self.Cyl_radius = 0
        else:
            self.Cyl_radius = _get_cylinder_radius(self.mesh_obj)
    
    def load_folder_tree(self):
        self.path_prj_hdy = os.path.join(self.prj_folder,'hydrodynamic')
        self.path_prj_dyn_res = os.path.join(self.path_prj_hdy,'results')
        self.path_prj_dyn_mesh = os.path.join(self.path_prj_hdy,'mesh')

    def gen_path(self):
        """
        gen_Path: generate the folder structure where the results are stored

        Optional args:
            hyd_dir (str):  name of the folder where the calculation is stored
            generate (boolean): trigger the generation on or off
        """
        self.load_folder_tree()
        os.makedirs(self.path_prj_hdy)
        os.makedirs(self.path_prj_dyn_res)
        os.makedirs(self.path_prj_dyn_mesh)


    def gen_mesh_files(self, show_norm=False):
        """
        gen_mesh_Files: generates the meshes to be used in the BEM solver

        Optional args:
            generate (boolean): triggers the generation on or off
            show_norm (boolean): triggers the visualisation of the mesh norms
        """
        # since Nemoh evaluate the hydrostatic matrix at the CG
        # it is possible to use that information rather then the translation of the body
        # the rotation is still needed though, to achieve a correct projection of the force components
        for nb in range(self.n_bodies):
            msh_obj = self.mesh_obj[nb]
            if self._debug: msh_obj.visualise_mesh()
            msh_obj.mesh_generation("nemoh",
                                    output_path=self.path_prj_dyn_mesh)
    
    def gen_multibody_structure(self):
        mb_obj = MultiBodyAnalysis(self.bodies, self.shared_dof_binary,  self.point_application)
        mb_obj.set_multibody_data()
        mb_obj.get_multibody_topology()
        mb_obj.set_joints_position_orientation()
        mb_obj.eval_velocity_transformation_matrix()
        self.bodies = mb_obj.bodies  # apply the new bodies numbering to the nemoh_obj
        self.vel_map = mb_obj.b_mat
        self.bodies_dofs = mb_obj.get_dofs_floating_bodies()
            
    def gen_hdyn_files(self):
        """
        gen_hdyn_Files: generates the files used to generate the hydrodynamic solution
    
        Optional args:
            generate (boolean): triggers the generation on or off
        """
        self.water_depth
        n_sh = len(self.bodies_dofs[0])
        
        f = open(os.path.join(self.path_prj_hdy,'ID.dat'), 'w')
        f.write('1\n')
        f.write('.')
        f.close()
    
        #input.txt
        f = open(os.path.join(self.path_prj_hdy,'input.txt'), 'w')
        f.write('--- Calculation parameters ------------------------------------------------------------------------------------\n')
        f.write('0				! Indiq_solver		! - 		! Solver (0) Direct Gauss (1) GMRES (2) GMRES with FMM acceleration (2 not implemented yet)\n')
        f.write('20				! IRES			! - 		! Restart parameter for GMRES\n')
        f.write('5.E-07				! TOL_GMRES		! -		! Stopping criterion for GMRES\n')
        f.write('100				! MAXIT			! - 		! Maximum iterations for GMRES\n')
        f.write('1				! Sav_potential		! -		! Save potential for visualization')
        f.close()
    
        #nemoh.cal
        f = open(os.path.join(self.path_prj_hdy,'nemoh.cal'), 'w')
        f.write('--- Environment ------------------------------------------------------------------------------------------------------------------\n')
        f.write('{:.2f}				! RHO 			! KG/M**3 	! Fluid specific volume \n'.format(self.rho))
        f.write('{:.2f}				! G			! M/S**2	! Gravity\n'.format(self.g))
        f.write('{}				! DEPTH			! M		! Water depth\n'.format(abs(self.water_depth)))#the water depth cannot be infinity
        f.write('0.	0.			! XEFF YEFF		! M		! Wave measurement point\n')
        f.write('--- Description of floating bodies -----------------------------------------------------------------------------------------------\n')
        f.write('{}				! Number of bodies\n'.format(self.n_bodies))
    
        for nb in range(self.n_bodies):
            mesh_fn = os.path.split(self.bodies[nb]['mesh'])[-1]
            f.write("----- Body {} ----------------------------------------------\n".format(self.bodies[nb]['ID']))
            f.write(".\mesh\{}dat			! Name of mesh file\n".format(mesh_fn[:-3]))
            f.write('{}\t{}			! Number of points and number of panels \n'.format(int(self.bodies[nb]['nV']),int(self.bodies[nb]['nP'])))
            f.write('{}				! Number of degrees of freedom\n'.format(n_sh))
            for dof in self.bodies_dofs[nb]:
                f.write('{:g} {} {} {} {} {} {}		! DOF type(1x1) direction(3x1) position(3x1)\n'.format(*dof))
    
            f.write('{}				! Number of generalised forces\n'.format(n_sh))
    
            for dof in self.bodies_dofs[nb]:
                f.write('{:g} {} {} {} {} {} {}		! DOF type(1x1) direction(3x1) position(3x1)\n'.format(*dof))
    
            f.write('0                                                    ! Number of lines of additional information\n')
    
        f.write('--- Load cases to be solved -------------------------------------------------------------------------------------------------------\n')
        f.write('{} {} {}		! Number of wave frequencies, Min, and Max (rad/s)\n'.format(self.frequency_def[0],self.frequency_def[1],self.frequency_def[2]))
        f.write('{} {} {}		! Number of wave angles, Min, and Max (deg)\n'.format(self.angle_def[0],self.angle_def[1],self.angle_def[2]))
        f.write('--- Post processing ---------------------------------------------------------------------------------------------------------------\n')
        f.write('0	0.01	10.		! IRF 				! IRF calculation (0 for no calculation), time step and duration\n')
        f.write('0				! Show pressure\n')
        f.write('0	0.	180.		! Kochin function 		! Number of directions of calculation (0 for no calculations), Min and Max (degrees)\n')
        f.write('0	2	1000.	2.	! Free surface elevation 	! Number of points in x direction (0 for no calcutions) and y direction and dimensions of domain in x and y direction\n')
        f.write('{}     {}     {}                                ! Cylindrical surface     ! Radius (0 for no calc), Number of points in theta direction, Number of points in z direction\n'.format(self.Cyl_radius,self.Cyl_Ntheta,self.Cyl_Nzeta))
        f.write(' --- End of file -------------------------------------------------------------------------------------------------------------------\n')
        f.write('\n')
        f.write('\n')
        f.close()

    def run_hst(self, start=True):
        """
        run_hst: run the hydrostatic solvers

        Optional args:
            nemoh_folder (str): location of the executable files
            start (boolean): triggers the calculation on or off
            generate (boolean): triggers the generation on or off
        """
        if start:
            directory = os.path.join(self.path_prj_hdy)
            Khst = Hydrostatics_Nemohcal(directory)
        
        return Khst

    def run_nemoh(self, nemoh_folder, start=True):
        """
        run_nemoh: run the hydrodynamic solvers

        Optional args:
            nemoh_folder (str): location of the executable files
            start (boolean): triggers the calculation on or off
            generate (boolean): triggers the generation on or off
        """
        
        if start:
            actual_dir = os.getcwd()
            os.chdir(self.path_prj_hdy)
            output_pre = execute(
                            os.path.join(nemoh_folder, 'preProcessor.exe')
                            )
            output_solver = execute(
                            os.path.join(nemoh_folder, 'solver.exe')
                            )
            output_post = execute(
                            os.path.join(nemoh_folder, 'postprocessor.exe')
                            )
            os.chdir(actual_dir)


def _get_cylinder_radius(meshes):
    
    coord = np.empty((0, 3))
    d = None
    
    for mesh in meshes:
        
        coord = np.vstack((coord, mesh.vertices))
        
        for panel in mesh.connectivity:
            
            d_panel = max([LA.norm(mesh.vertices[panel[0]] -
                                   mesh.vertices[panel[s]]) for s in [1, -1]])
            
            if d is None or d_panel > d:
                d = d_panel
    
    if d is None:
        raise RuntimeError("Can not determine mesh extents")
    
    centroid = coord.mean(0)
    
    if not np.sqrt(np.sum(centroid[:2] ** 2)) / coord.max() < 1e-3:
        module_logger.warning(
            "WARNING: the centroid of the mesh file is not centered "
            "at the origin of the mesh coordinate system.\n"
            "Consider regenerating a mesh to satisfy this condition.")
    
    return (np.sqrt(coord[:, 0] ** 2 + coord[:, 1] ** 2)).max() + d


def rot_matrix(ang):
    Rz = lambda angle: np.array([[np.cos(angle), -np.sin(angle), 0], [np.sin(angle), np.cos(angle), 0], [0, 0, 1]])
    Ry = lambda angle: np.array([[np.cos(angle), 0, np.sin(angle)], [0, 1, 0], [-np.sin(angle), 0, np.cos(angle)]])
    Rx = lambda angle: np.array([[1, 0, 0], [0, np.cos(angle), -np.sin(angle)], [0, np.sin(angle), np.cos(angle)]])
    
    return Rz(ang[2])*Ry(ang[1])*Rx(ang[0])

def execute(command):
    """
    unused TBD
    """
    import sys

    process = subprocess.Popen(command, stdout=subprocess.PIPE)
    for c in iter(lambda: process.stdout.read(1), ''):
        sys.stdout.write(c)
#    process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
#    output = ''
#
#    # Poll process for new output until finished
#    for line in iter(process.stdout.readline, ""):
#        module_logger.info(line),
#        output += line
#
#
#    process.wait()
#    exitCode = process.returncode
#
#    if (exitCode == 0):
#        return output
#    else:
#        raise Exception(command, exitCode, output)

if __name__ == "__main__":
    import sys
    sys.path.append(r"C:\Users\francesco\Desktop\test_gui\utils")
    from data_interface import DataStructure
    import dtocean_wave.utils.hdf5_interface as h5i
    data = h5i.load_dict_from_hdf5(r"C:\Users\francesco\Desktop\test_gui\test_prj\test_prj_data_collection.hdf5")
    
    dataobj = DataStructure(data)
    dataobj.body_inputs['body']['body0']['mesh'] = os.path.join(u'C:\\Users\\francesco\\Desktop\\test_gui', dataobj.body_inputs['body']['body0']['mesh'])
    dataobj.body_inputs['body']['body1']['mesh'] = os.path.join(u'C:\\Users\\francesco\\Desktop\\test_gui', dataobj.body_inputs['body']['body1']['mesh'])
    bem_obj = NemohExecute(dataobj.project_folder,
                                   dataobj.general_inputs,
                                   dataobj.body_inputs,
                                   get_array_mat=False,
                                   debug=False)
    bem_obj.gen_path()
    bem_obj.gen_mesh_files()
    bem_obj.gen_multibody_structure()
    bem_obj.gen_hdyn_files()
    bem_obj.run_nemoh()
    bem_obj.run_hst()
    