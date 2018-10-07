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
This module contains the classes used map the loose dofs into the generalised dof space

.. module:: input
   :platform: Windows
   :synopsis: Input module to DTOcean WP2

.. moduleauthor:: Francesco Ferri <ff@civil.aau.dk>
"""

import logging

import numpy as np
import utils.eom_multibody as eom_mb

# Start module_logger
module_logger = logging.getLogger(__name__)


class MultibodyMapping():
    """
    Multibody_mapping: the class is used to map from the loose degree of freedom to the generalised degree of freedom
                        for the case of a multibody WEC.

    Args:
        NemohObj (BEM_Nemoh class): the object contains the specification of the hydrodynamic model and bodies connectivity

    Optional args:
        J (numpy.ndarray): Jacobian matrix (linearised) used in the hydrostatic mapping

    Attributes:
        dofList (list): list of shared and extended degree of freedom
        dof (int): total number of generalised dofs
        Nb (int): number of bodies in the WEC
        nshared (int): number of shared dofs
        J (numpy.ndarray): Jacobian matrix (linearised) used in the hydrostatic mapping, implicit application of mechanical constraints
        RR (numpy.ndarray): composite rotation matrix
        c (numpy.ndarray): composite hydrostatic stiffness matrix built from the bodies matrix
        LocalAxes (list): list of local coordinate system
        lCS (numpy.ndarray): local coordinate system of the platform wrt the 0,0,0 of the mesh file
        connectivity (list): connectivity between bodies in the multibody system
        looseDOF (int): number of degree of freedom without considering constraints
        Khyd_mapped (numpy.ndarray): hydrostatic matrix defined in the generalised dof space

    Returns:
        Khyd_mapped (numpy.ndarray): hydrostatic matrix defined in the generalised dof space
    """
    def __init__(self, NemohObj, J=np.empty(0)):
        self.mapping =  NemohObj._mapping
        self.bodies = NemohObj.bodies
        self.shared_dof = NemohObj.sharedDOF
        self.l_cs = NemohObj.lCS
        self.inertia = [el[2] for el in NemohObj.mass_property]
        self.mass = [el[1] for el in NemohObj.mass_property]
        
        
    def _eval_reduced_mass_matrix(self):
        lCS = self.l_cs
        shared_dof = self.shared_dof
        bodies = self.bodies
        m = self.mass
        inertia = self.inertia
        
        # restrucure the inputs
        bodies, lCS, ndof, q, dq, x, t, shared_dof = eom_mb.set_data_for_mb_analysis(bodies, lCS, shared_dof)
        platform_tr, ind = eom_mb.set_platform_translational_dof(shared_dof, q)
        Ir0P1, kn_r_PnPn1 = eom_mb.set_global_poa_and_local_cs(bodies, platform_tr, lCS)
        KnomegaIKn, dof_rot4mat, rot_matrix = eom_mb.set_local_angular_velocity(bodies, shared_dof, q, dq, ind)         
                       
        # global positions of the dofs
        Ir0Pn, Iomega0Pn = eom_mb.set_global_state_of_application_points(bodies, rot_matrix, dof_rot4mat, Ir0P1)
        Ir0Sn, Iomega0Sn = eom_mb.set_global_state_of_cog(bodies, rot_matrix, Ir0Pn, Iomega0Pn)
        IvSn = eom_mb.set_global_velocity_of_cog(Ir0Sn, q, dq, t)
        
        T = eom_mb.evaluate_kinetic_energy(IvSn, KnomegaIKn, m, inertia)    
        
        M = eom_mb.lagrangian_equation(T, q, dq)
        
        Mlin = eom_mb.linearisation(M, q, x)
 
        return Mlin
    
    def _eval_hydrodynamic(self,Khst,Ca,Cm,Fex,Radiation):
        """
        _eval_hydrodynamic: evaluate the mapping for the excitation force, radiation force and partial waves

        Args:
            Ca (numpy.ndarray): Radiaton damping matrix function of the dofs and frequencies
            Cm (numpy.ndarray): Added mass matrix function of the dofs and frequencies
            Fex (numpy.ndarray): Excitation force matrix function of the dofs, directions and frequencies
            Radiation (numpy.ndarray): Radiaton partial waves function of the dofs and frequencies

        Returns:
            Ca_r (numpy.ndarray): Mapped radiaton damping matrix function of the dofs and frequencies
            Cm_r (numpy.ndarray): Mapped added mass matrix function of the dofs and frequencies
            Fex_r (numpy.ndarray): Mapped excitation force matrix function of the dofs, directions and frequencies
            Radiation_r (numpy.ndarray): Mapped radiaton partial waves function of the dofs and frequencies
        """
        Map = self.mapping
        
        Khst_r = np.dot(Map,np.dot(Khst,Map.T))
        Ca_r = np.dot(Map,np.dot(Ca,Map.T).transpose(2,1,0)).transpose(2,1,0)
        Cm_r = np.dot(Map,np.dot(Cm,Map.T).transpose(2,1,0)).transpose(2,1,0)
        Fex_r = np.dot(Fex,Map.T)
        Radiation_r = np.dot(Radiation.transpose(0,2,3,1),Map.T).transpose(0,3,1,2)
        
        return Khst_r, Ca_r, Cm_r, Fex_r, Radiation_r

    

if __name__ == "__main__":
    class emptyClass():
        def __init__(self):
            #print('Empty class created')
            return
    
    case = 5
    khy = [np.array([[0.0000000E+00,  0.0000000E+00,  0.0000000E+00,  0.0000000E+00,  0.0000000E+00,  0.0000000E+00],
              [0.0000000E+00,  0.0000000E+00,  0.0000000E+00,  0.0000000E+00,  0.0000000E+00,  0.0000000E+00],
              [0.0000000E+00,  0.0000000E+00,  0.7892646E+06,  0.2148438E-01,  0.5993652E-01,  0.0000000E+00],
              [0.0000000E+00,  0.0000000E+00,  0.2148438E-01, -0.3444551E+08,  0.0000000E+00,  0.0000000E+00],
              [0.0000000E+00,  0.0000000E+00,  0.5993652E-01,  0.0000000E+00, -0.3444972E+08,  0.0000000E+00],
              [0.0000000E+00,  0.0000000E+00,  0.0000000E+00,  0.0000000E+00,  0.0000000E+00,  0.0000000E+00]],'f')]


    if case==1:
        #print('Test 1')
        #print('Three body snake, the reference body is the central body')
        #print('+'*100)
        Ne = emptyClass()
        Ne.nDOF_b  =np.array([[ 3.,  4., 4.]])
        Ne.totDOF = 5
        Ne.N_bodies = 3
        Ne.connectivity = [0,0,0]

        Ne.Khst = khy*Ne.N_bodies


        Ne.lCS = np.array([0,0,0.])

        Ne.LocalAxes = [np.array([0,0,0,0.]),
                       np.array([0,10,0.,0]),
                       np.array([0,-10,0.,0])]


        Ne.sharedDOF = [np.array([1, 1, 0, 0, 0, 0, 0]),
                       np.array([1, 0, 0, 1, 0, 0, 0]),
                       np.array([2, 0, 1., 0, 0, 0, 0])]
        Ne.extendedDOF = [np.array([2 ,0, 1., 0, 5., 0., 0 ,1]),
                         np.array([2 ,0, 1., 0, -5., 0., 0 ,2])]

    elif case==2:
        #print('Test 2')
        #print('Three body snake, the reference body is head of the system')
        #print('+'*100)
        Ne = emptyClass()
        Ne.nDOF_b  =np.array([[ 3.,  4., 4.]])
        Ne.totDOF = 5
        Ne.N_bodies = 3
        Ne.connectivity = [0,0,1]

        Ne.Khst = khy*Ne.N_bodies


        Ne.lCS = np.array([-10,0,0.])

        Ne.LocalAxes = [np.array([0,-10,0,0.]),
                       np.array([0,0,0.,0]),
                       np.array([0,10,0.,0])]


        Ne.sharedDOF = [np.array([1, 1, 0, 0, 0, 0, 0]),
                       np.array([1, 0, 0, 1, 0, 0, 0]),
                       np.array([2, 0, 1., 0, 0, 0, 0])]
        Ne.extendedDOF = [np.array([2 ,0, 1., 0, -5., 0., 0 ,1]),
                         np.array([2 ,0, 1., 0, 5., 0., 0 ,2])]

    elif case==3:

        #print('Test 3')
        #print('Three body snake with last body with inclined axes (45°)')
        #print('+'*100)
        Ne = emptyClass()
        Ne.nDOF_b  =np.array([[ 3.,  4., 4.]])
        Ne.totDOF = 5
        Ne.N_bodies = 3
        Ne.connectivity = [0,0,1]

        Ne.Khst = khy*Ne.N_bodies



        Ne.lCS = np.array([0,0,0.])

        Ne.LocalAxes = [np.array([0,0,0,0.]),
                       np.array([0,10,0.,0]),
                       np.array([np.pi/4,10.+10./np.sqrt(2),10.+10./np.sqrt(2),0])]


        Ne.sharedDOF = [np.array([1, 1, 0, 0, 0, 0, 0]),
                       np.array([1, 0, 0, 1, 0, 0, 0]),
                       np.array([2, 0, 1., 0, 0, 0, 0])]
        Ne.extendedDOF = [np.array([2 ,0, 1., 0, 5., 0., 0 ,1]),
                         np.array([2 ,-1./np.sqrt(2), 1./np.sqrt(2), 0, 10+5./np.sqrt(2), 10+5./np.sqrt(2), 0 ,2])]

    elif case==4:

        #print('Test 4')
        #print('Three body snake with last two bodies with inclined axes (45°)')
        #print('+'*100)

        Ne = emptyClass()
        Ne.nDOF_b  =np.array([[ 3.,  4., 4.]])
        Ne.totDOF = 5
        Ne.N_bodies = 3
        Ne.connectivity = [0,0,1]

        khy = [np.array([[0.0000000E+00,  0.0000000E+00,  0.0000000E+00,  0.0000000E+00,  0.0000000E+00,  0.0000000E+00],
               [0.0000000E+00,  0.0000000E+00,  0.0000000E+00,  0.0000000E+00,  0.0000000E+00,  0.0000000E+00],
               [0.0000000E+00,  0.0000000E+00,  0.7892646E+06,  0.2148438E-01,  0.5993652E-01,  0.0000000E+00],
               [0.0000000E+00,  0.0000000E+00,  0.2148438E-01, -0.3444551E+08,  0.0000000E+00,  0.0000000E+00],
               [0.0000000E+00,  0.0000000E+00,  0.5993652E-01,  0.0000000E+00, -0.3444972E+08,  0.0000000E+00],
               [0.0000000E+00,  0.0000000E+00,  0.0000000E+00,  0.0000000E+00,  0.0000000E+00,  0.0000000E+00]],'f')]*Ne.N_bodies
        Ne.Khst = khy

        # position of the local coordinate system of the platform with respect to the
        # center of buoyancy
        Ne.lCS = np.array([0,0,0.])

        # position and orientation of the bodies'axis with respect to the platform
        Ne.LocalAxes = [np.array([0,0,0,0.]),
                        np.array([np.pi/4.,5.+5./np.sqrt(2),5./np.sqrt(2),0]),
                        np.array([np.pi/4.,5.+15./np.sqrt(2),15./np.sqrt(2),0])]


        Ne.sharedDOF = [np.array([1, 1, 0, 0, 0, 0, 0]),
                        np.array([1, 0, 0, 1, 0, 0, 0]),
                        np.array([2, 0, 1., 0, 0, 0, 0])]
        Ne.extendedDOF = [np.array([2 ,-1./np.sqrt(2), 1./np.sqrt(2),0, 5., 0., 0 ,1]),
                          np.array([2 ,-1./np.sqrt(2), 1./np.sqrt(2),0, 5.+10/np.sqrt(2), 10./np.sqrt(2), 0 ,2])]

        
    elif case==5:

        #print('Test 5')
        #print('Two body snake')
        #print('+'*100)

        Ne = emptyClass()
        Ne.nDOF_b  =np.array([[ 3.,  4.]])
        Ne.totDOF = 4
        Ne.N_bodies = 2
        Ne.connectivity = [0,0]

        khy = [np.array([[0.0000000E+00,  0.0000000E+00,  0.0000000E+00,  0.0000000E+00,  0.0000000E+00,  0.0000000E+00],
               [0.0000000E+00,  0.0000000E+00,  0.0000000E+00,  0.0000000E+00,  0.0000000E+00,  0.0000000E+00],
               [0.0000000E+00,  0.0000000E+00,  0.7892646E+06,  0.2148438E-01,  0.5993652E-01,  0.0000000E+00],
               [0.0000000E+00,  0.0000000E+00,  0.2148438E-01, -0.3444551E+08,  0.0000000E+00,  0.0000000E+00],
               [0.0000000E+00,  0.0000000E+00,  0.5993652E-01,  0.0000000E+00, -0.3444972E+08,  0.0000000E+00],
               [0.0000000E+00,  0.0000000E+00,  0.0000000E+00,  0.0000000E+00,  0.0000000E+00,  0.0000000E+00]],'f')]*Ne.N_bodies
        Ne.Khst = khy

        # position of the local coordinate system of the platform with respect to the
        # center of buoyancy
        Ne.lCS = np.array([-4,0,0.])

        # position and orientation of the bodies'axis with respect to the platform
        Ne.LocalAxes = [np.array([0, -4, 0, 0.]),
                        np.array([0, 6, 0, 0.])]

        Ne.sharedDOF = [np.array([1, 1, 0, 0, 0, 0, 0]),
                        np.array([1, 0, 0, 1, 0, 0, 0]),
                        np.array([2, 0, 1., 0, 0, 0, 0])]
        Ne.extendedDOF = [np.array([2, 0., 1., 0, 1., 0., 0, 2.])]


    MB = Multibody_mapping(Ne)
    MB._eval_hydrostatic(Ne.Khst)

    #print('Jacobian Matrix')
    #print('*'*100)
    #print(MB.J)
    #print('*'*100)
    #print(' ')
    #print('Mapped hydrostatic stiffness')
    #print('*'*100)
    #print(MB.Khyd_mapped)
    #print('*'*100)
    

