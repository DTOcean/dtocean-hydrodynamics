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

# Start module_logger
import logging
module_logger = logging.getLogger(__name__)

from scipy.linalg import block_diag,expm3,norm
import numpy as np

class Multibody_mapping():
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
        self.dofList =  list(NemohObj.sharedDOF)+list(NemohObj.extendedDOF)
        self.dof = NemohObj.totDOF
        self.Nb = NemohObj.N_bodies
        self.nshared = len(list(NemohObj.sharedDOF))
        
        if J.size==0:
            pass
        else:
            pass
        
        self.J = np.zeros((self.Nb*6,self.dof))
        self.RR = np.empty(0)
        self.c = np.empty(0)
        self.LocalAxes = NemohObj.LocalAxes
        self.lCS = NemohObj.lCS
        #build the connectivity from the platfomr to the end effector of each body
        self.connectivity = []
        Body_connectivity = ()
        for b in range(self.Nb):
            if NemohObj.connectivity[b]==0:
                Body_connectivity = (NemohObj.connectivity[b])
            else:
                connection = NemohObj.connectivity[b]
                Body_connectivity = ()
                while connection:
                    Body_connectivity += (connection,)
                    connection = NemohObj.connectivity[connection]
                Body_connectivity += (0,)
                Body_connectivity = Body_connectivity[::-1]
                    
                    
            self.connectivity += [Body_connectivity]
            

        self.looseDOF = NemohObj.nDOF_b.sum()
    
    def _eval_hydrodynamic(self,Ca,Cm,Fex,Radiation):
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
        totDOF = self.dof
        looseDOF = self.looseDOF
        nshared = self.nshared
        
        Map = np.eye(totDOF,nshared)
        for exDOF in range(nshared,totDOF):
            map_ex = np.zeros((totDOF,nshared+1))
            map_ex[:,:nshared] = np.eye(totDOF,nshared)
            map_ex[exDOF,-1] = 1
            Map = np.hstack((Map,map_ex))
        
        Ca_r = np.dot(Map,np.dot(Ca,Map.T).transpose(2,1,0)).transpose(2,1,0)
        Cm_r = np.dot(Map,np.dot(Cm,Map.T).transpose(2,1,0)).transpose(2,1,0)
        Fex_r = np.dot(Fex,Map.T)
        Radiation_r = np.dot(Radiation.transpose(0,2,3,1),Map.T).transpose(0,3,1,2)
        
        return Ca_r, Cm_r, Fex_r, Radiation_r

    def _eval_hydrostatic(self, Khyd):
        """
        _eval_hydrostatic: evaluate the mapping for the hydrostatic matrix

        Args:
            Khyd (list): list of hydrostatic matrix for each body

        Return:
            Khyd_r (numpy.ndarray): mapped hydrostatic matrix
        """
        
        def M(axis,theta): return expm3(np.cross(np.eye(3),axis/norm(axis)*theta))
        #Explanation of the M matrix.
        #if a is an unit vector of axis such that a = axis/nomr(axis)
        #and A = Ia is teh skew-symmetric matrix associated to a
        #then M = exp(theta, A) is the rotation matrix
        #expm3 computes the taylor series of the exponential: \sum_{k=0}^{2theta} \frac{1}{k!} (theta A)^k
        
        ex_dof = np.array(self.dofList[self.nshared:])        
        
        #create the block diagonal matrix for rotation and hydrostatic stiffness
        for b in range(self.Nb):
            R = block_diag(M([0,0,1],-self.LocalAxes[b][0]),M([0,0,1],-self.LocalAxes[b][0]))
            if b==0:
                self.RR = R
                self.c = Khyd[b]
            else:
                self.RR = block_diag(self.RR,R)
                self.c = block_diag(self.c,Khyd[b])
                
        Z = np.zeros((3))
        Jplatform = np.zeros((6,self.dof))
        
        #first solve the platform case
        module_logger.info('Calculating the Jacobian for the platform')
        for dof_i in range(self.dof):
            if len(self.dofList[dof_i])==7:
                k = self.dofList[dof_i][1:4]
                if self.dofList[dof_i][0]==1:
                    Jplatform[:6,dof_i] = np.r_[k,Z]
                else:    
                    r = self.LocalAxes[0][1:]-self.lCS
                    Jplatform[:6,dof_i] = np.r_[np.cross(k,r),k]
        
        self.J[:6,:] = Jplatform
        #print(Jplatform)
        
        
        for b in range(1,self.Nb):
            Jbody = np.zeros((6,self.dof))
            module_logger.info('Calculating the Jacobian for body #{}'.format(b))
            for dof_i in range(self.dof):
                if len(self.dofList[dof_i])==7:#common degree of freedom
                    #k is the dof axis
                    k = self.dofList[dof_i][1:4]
                    #differentiate between translation (1) and rotations (2)
                    if self.dofList[dof_i][0]==1:
                        Jbody[0:6,dof_i] = np.r_[k,Z]
                    else:    
                        r = self.LocalAxes[b][1:]-self.lCS
                        Jbody[0:6,dof_i] = np.r_[np.cross(k,r),k]
                        
                        
                        
                else:#extended degree of freedoms
                    if b==self.dofList[dof_i][7]-1:
                        #check if the body is connected to the platform
                        #by looking at the connectivity list   
                        if self.connectivity[b]==0:
                            k = self.dofList[dof_i][1:4]
                            pq = self.dofList[dof_i][4:7]-self.lCS
                            if self.dofList[dof_i][0]==1:
                                Jbody[0:6,dof_i] = np.r_[k,Z]
                            else:    
                                r = self.LocalAxes[b][1:]-self.lCS
                                Jbody[0:6,dof_i] = np.r_[np.cross(k,r-pq),k]
                        else:
                            for el in self.connectivity[b]:
                                
                                if el==0:
                                    k = self.dofList[dof_i][1:4]
                                    pq = self.dofList[dof_i][4:7]-self.lCS
                                    if self.dofList[dof_i][0]==1:
                                        Jbody[0:6,dof_i] = np.r_[k,Z]
                                    else:    
                                        r = self.LocalAxes[b][1:]-self.lCS
                                        Jbody[0:6,dof_i] = np.r_[np.cross(k,r-pq),k]
                                        
                                else:
                                    
                                    dof = ex_dof[:,7] == el
                                    ex_dof_i = [i for i,v in enumerate(dof) if v][0]+self.nshared
                                     
                                    k = ex_dof[dof,:][0][1:4]
                                    pq = ex_dof[dof,:][0][4:7]-self.lCS
                                    
                                    if ex_dof[dof,:][0][0]==1:
                                        Jbody[0:6,ex_dof_i] = np.r_[k,Z]
                                    else:    
                                        r = self.LocalAxes[b][1:]-self.lCS
                                        Jbody[0:6,ex_dof_i] = np.r_[np.cross(k,r-pq),k]
                            
            
            self.J[b*6:(b+1)*6,:] = Jbody
            #print(Jbody)
        self.Khyd_mapped = np.dot(self.J.T,np.dot(self.RR.T,np.dot(self.c,np.dot(self.RR,self.J))))
        
        return self.Khyd_mapped

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
    

