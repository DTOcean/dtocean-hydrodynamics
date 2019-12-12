# -*- coding: utf-8 -*-

#    Copyright (C) 2016 Francesco Ferri, Pau Mercadez Ruiz
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
This module contains the classes used to collect and modify the WP2 inputs

.. module:: input
   :platform: Windows
   :synopsis: Input module to DTOcean WP2

.. moduleauthor:: Francesco Ferri <ff@civil.aau.dk>
.. moduleauthor:: Pau Mercadez Ruiz <pmr@civil.aau.dk>
.. moduleauthor:: Mathew Topper <mathew.topper@dataonlygreater.com>
"""

from math import pi
from itertools import izip

import numpy as np
from numpy import (array,
                   zeros,
                   ones,
                   transpose,
                   dot,
                   sqrt,
                   arctan2,
                   exp,
                   cos,
                   sin,
                   linspace,
                   eye,
                   reshape,
                   meshgrid,
                   real,
                   imag)
from numpy.linalg import solve
from scipy.special import jv, yv

from dtocean_hydro.output import ReducedOutput

from .utils.WatWaves import len2
from .utils.StrDyn import EnergyProduction
from .utils.spec_class import wave_spec

# Start logging
import logging
module_logger = logging.getLogger(__name__)


class MultiBody(object):
    """
    MultiBody: the class assess the hydrodynamic model for an array of WECs, given the a WEC object in input

    Args:
        iHydro (Hydro_pkg class): object containing the specification of the sea states, frequencies and directions to be analysed
        iwec (wec class): object containing the description of the single body hydrodynamic

    Optional args:
        cylamplitude (boolean): used to save or not the partial wave, associated with the array

    Attributes:
        coord (numpy.ndarray) [m]: UTM coordinates (Easting-Northing) describing the array layout.
        iHydro (Hydro_pkg class): object containing the specification of the sea states, frequencies and directions to be analysed
        cylamplitude (boolean): used to save or not the partial wave, associated with the array
        depth (float) [m]: average water depth at the site
        Hs (numpy.ndarray) [m]: vector containing the wave heights associated to the scatter diagram
        Tp (numpy.ndarray) [s]: vector containing the wave periods associated to the scatter diagram
        B (numpy.ndarray) [rad]: vector containing the wave directions associated to the scatter diagram
        ScatDiag (numpy.ndarray): normalised probability of occurence of the sea states
        ScatDiag_spec (tuple): the tuple gathers the spectrum type, the frequency and the directional spreading parameters
        period (numpy.ndarray) [s]: vector containing the wave periods to be analysed
        wnumber (numpy.ndarray) [rad/m]: vector containing the wave numbers to be analysed
        dir (numpy.ndarray) [rad]: vector containing the wave directions to be analysed
        energyproduction (float): average energy production for the given sea states
        Fex (numpy.ndarray): excitation force in function of the dof, frequency and direction
        Madd (numpy.ndarray): added mass in function of the dof and frequency
        Crad (numpy.ndarray): radiation damping in function of the dof and frequency
        Khyd (numpy.ndarray): hydrostatic matrix function of the dofs
        M (numpy.ndarray): mass matrix function of the dofs
        CPTO (numpy.ndarray): PTO damping function of the sea states and dof for the array
        Kmooring (numpy.ndarray): mooring stiffness function of the sea states and dof for the array
        Cfit (numpy.ndarray): damping fitting, calculated to reduce the error between the cetified and calculated power matrixes
        Kfit (numpy.ndarray): stiffness fitting, calculated to reduce the error between the cetified and calculated power matrixes
        power_prod_perD_perS (numpy.ndarray) [W]: power production per device per sea states
        power_prod_perD (numpy.ndarray) [W]: power production per device
        AEP_perD_perS (numpy.ndarray) [Wh]: annual energy production per device per sea states
        AEP_perD (numpy.ndarray) [Wh]: annual energy production per device
        AEP_array (float) [Wh]: annual energy production of the array
        q_perD (numpy.ndarray): q-factor per device
        q_array (numpy.ndarray): average q-factor for the array
        Nbodies (int): number of bodies in the array
        Device_Position (dic) [m]: gathers the UTM location of each machine, the dictionary keys are "machineID", where ID is the
                                machine ith number.
        Resource_reduction (float) [-]: ratio of outgoing energy over the incoming energy, for the worst case only
        TI (-): unused
    """
    def __init__(self, ihydro, iwec, cylamplitude=True, debug=False):

        # coordinates of WECs
        self.coord = None
        self.debug = debug
        # cylamplitude = True ==> will save cylindrical amplitude coefficients
        self.cylamplitude = cylamplitude
        self.depth = ihydro.depth
        # Sea States(Hs, Tp, B)
        self.Hs = ihydro.Hs
        self.Tp = ihydro.Tp
        self.B = ihydro.B
        # Probability of Occurrence of each sea state
        self.ScatDiag = ihydro.ScatDiag
        self.ScatDiag_spec = ihydro.specType
        # Analysis of regular waves or
        # Discretization of sea states for irregula waves analysis
        self.period = ihydro.period
        self.wnumber = ihydro.wnumber
        self.dir = ihydro.dir
        self.iwec = iwec
        (self.cpto,
         self.cfit,
         self.kmooring,
         self.kfit) = iwec.matrix_zoh_interp(ihydro)
        
        self._mooring_fex_index = self.estimate_design_wave()
    
    def estimate_design_wave(self):
        """
        estimate_design_wave: asses the index of the design wave to be passed to the mooring and 
        fundation module.
        This apporach is a quick-fix for an input which as not be requested
        to the user by the mooring and fundation module, but it is far to be considered accurate or exaustive.
        
        WARNING:
        Since the hydrodynamic module has been built around operational conditions, 
        this method, used to estimate the desing wave, is extremely inaccurate because the scatter diagram will
        not cover the extreme sea states.
        In the opinion of the module author a better apporach will be to ask the user to specify the design wave himself/herself.
        """
        ihs = np.argmax(self.Hs)
        itp = np.argmax(self.Tp)
        
        return (ihs, itp, 0)

    def energy(self, coord):
        """
        energy: calculates the energy production of the given array layout

        Args:
            coordinates (numpy.ndarray):
                UTM coordinates (Easting-Northing) describing the array layout.

        Returns:
            res (ReducedOutput class):
                contains the minimum set of information required by the
                WP2Output class.
        """
        self.coord = coord

        # Handy names
        NBo = len2(self.coord)

        # loop around orientations, O, with the corresponding
        # wave directions, Os, to analyse.
        # initialize
        Pyradd = zeros((NBo,) + self.ScatDiag.shape, dtype=float)
        Padd = Pyradd.copy()
        PyraddWEC = zeros((1,) + self.ScatDiag.shape, dtype=float)
        PaddWEC = PyraddWEC.copy()

        Nm = array(range(-self.iwec.order.max(), self.iwec.order.max() + 1),
                   dtype=float)
        m, n = np.meshgrid(Nm, Nm, indexing='ij')
        Gp = np.transpose(self.iwec.G, axes = (0,2,1))
        wec_Fex = None
        wec_dir = None
                
        orients = self.dir[::2]
        theta_groups = self.dir[1::2]
                
        for orient, thetas in izip(orients, theta_groups):
            
            ths = array([th for th in thetas] , dtype=float)

            d_rot = self.iwec.D*np.exp(-1j*(n-m)*orient)
            g_rot = np.transpose(Gp*np.exp(1j*Nm*orient), axes=(0,2,1))
            ar_rot = self.iwec.AR*np.exp(1j*Nm*orient)
            
            # solve for Madd, Crad and Fex
            self.Hydrodynamics(self.iwec, ths, d_rot, g_rot, ar_rot)
            
            Pyr, P = EnergyProduction(NBo,
                                      self.B,
                                      self.Hs,
                                      self.Tp,
                                      thetas,
                                      self.period,
                                      [self.ScatDiag, self.ScatDiag_spec],
                                      self.iwec.M,
                                      self.Madd,
                                      self.cpto,
                                      self.Crad,
                                      self.kmooring,
                                      self.iwec.Khyd,
                                      self.Fex,
                                      self.kfit,
                                      self.cfit,
                                      self.iwec.rated_power)

            # Excitation Force isolated WEC
            (mths,mode) = meshgrid(ths, Nm, indexing='ij')
            AP = exp(-1j*mode*(pi/2+mths))
            Fex_iso = zeros((self.Fex.shape[:2] + (self.iwec.Fex.shape[2],)),
                            dtype=np.complex64)
            for f in range(len2(self.period)):
                Fex_iso[f,:,:]= dot(AP,g_rot[f])
            
            if wec_Fex is None:
                wec_Fex = Fex_iso[:,0,:]
                wec_dir = np.array([ths[0]-orient])
            
            PyrWEC, PWEC = EnergyProduction(1,
                                            self.B,
                                            self.Hs,
                                            self.Tp,
                                            thetas,
                                            self.period,
                                            [self.ScatDiag,
                                             self.ScatDiag_spec],
                                            self.iwec.M,
                                            self.iwec.Madd,
                                            self.cpto,
                                            self.iwec.Crad,
                                            self.kmooring,
                                            self.iwec.Khyd,
                                            Fex_iso,
                                            self.kfit,
                                            self.cfit,
                                            self.iwec.rated_power)

            Pyradd += Pyr
            Padd += P

            PyraddWEC += PyrWEC
            PaddWEC += PWEC

        self.energyproduction = Pyradd.sum((1,2,3))
        self.power = Padd

        AEP_single = ((PyraddWEC)*365*24).sum()

        # new Device_Model format as agreed with Sam in Cork
        # Mode matrix (e.g. [[1, 0, 0],[0, 0, 0],[0, 0, 0]] if only surge has been analysed)
        # Wave frequencies (list)
        # Wave heights (list)
        # Excitation force (LxMxN, for L frequencies, M directions and N degrees of freedom

        tNb, tNt, tNh, tNd = Padd.shape
        Psd = zeros((tNt*tNh*tNd,tNb)) 
        indexFlat = -1
        for idi in range(tNd):
            for ih in range(tNh):
                for it in range(tNt):
                    indexFlat += 1
                    Psd[indexFlat,:] = Padd[:,it,ih,idi]
                    
        energy_balance = self.EnergyBalance(Psd.T)
        aep_dev = Pyradd.sum((1,2,3))*365*24
        translational_modes = [0,0,0]
        Fex_for_mooring = [[],[],[]]
        
        for mode_ind, mode_val in enumerate(self.iwec.modes):
            if mode_val[1:4].tolist() == [1,0,0]:
                Fex_for_mooring[0] = wec_Fex[:,mode_ind].tolist()
                translational_modes[0] = 1
            elif mode_val[1:4].tolist() == [0,1,0]:
                Fex_for_mooring[1] = wec_Fex[:,mode_ind].tolist()
                translational_modes[1] = 1
            elif mode_val[1:4].tolist() == [0,0,1]:
                Fex_for_mooring[2] = wec_Fex[:,mode_ind].tolist()
                translational_modes[2] = 1
                
        device_model = {'wave_fr': 1 / self.period,
                        'mode_def': translational_modes,
                        'wave_dir': wec_dir,
                        'fex':Fex_for_mooring}
        
        machine = {}
        for jj in range(NBo):
            strn = 'Device%d'%jj
            machine[strn] = (self.coord[jj, 0], self.coord[jj, 1])

        if self.debug:  # add attributes for inspection
            self.power_prod_perD_perS= Psd.T
            self.power_prod_perD = Pyradd.sum((1, 2, 3))
            self.AEP_perD_perS = Pyradd*365*24
            self.AEP_perD = Pyradd.sum((1,2,3))*365*24
            self.AEP_array = aep_dev.sum()
            self.q_perD = aep_dev/AEP_single
            self.q_array = self.AEP_array/(AEP_single*NBo)
            self.Nbodies = NBo
            self.Device_Model = device_model
            self.Device_Position = machine
            self.Resource_reduction = energy_balance
            self.TI = []

        res = ReducedOutput(aep_dev.sum(),
                            aep_dev,
                            aep_dev.sum() / (AEP_single * NBo),
                            aep_dev / AEP_single,
                            Pyradd.sum((1, 2, 3)),
                            Psd.T,
                            NBo,
                            machine,
                            energy_balance,
                            None,
                            device_model,
                            PaddWEC[0,:,:,:],
                            None)
        
        return res

    def Hydrodynamics(self, iWEC, direction, D_rot, G_rot, AR_rot):
        """
        Hydrodynamics: Computes hydrodynamic interactions between the bodies.

        Args:
            iWEC (WEC class): It is an instance of WEC class.
            direction (numpy.ndarray) [rad]: It contains wave-directions considered for
                                                diffraction problems.
            order (int): Truncation order for the power series solution of cylindric waves (Kagemoto).
        """
        # Check point {start}"
        TolCheck = 1e-3
        requires = [abs(self.period-iWEC.period) < TolCheck, 
                    [abs(self.depth-iWEC.depth) < TolCheck]]
        
        requires = [all(requires[ind]) for ind in range(len(requires))]
                    
        if not all(requires):
#            print 'period {} , depth {} '.format(requires[0],requires[1])
            raise IOError('Mismatch between iHydro and retreieved sea state '
                          'data. iWECpickler may be stale.')
        
        k0 = self.wnumber
        cfreq = 2*pi/self.period
        
        # Interaction theory from Kagemoto
        (TS,
         TR,
         DS_array,
         DR_array,
         GS_array,
         GR_array,
         MS_inter,
         MR_inter) = self.Interaction(k0,
                                      D_rot,
                                      G_rot,
                                      iWEC.truncorder)
    
        # Wave scattering between bodies
        self.Scattering(k0,
                        direction,
                        TS,
                        DS_array,
                        GS_array,
                        MS_inter,
                        iWEC.order[0])
        
        # Radiated waves and scattering between bodies
        self.Radiation(k0,
                       cfreq,
                       TR,
                       DR_array,
                       GR_array,
                       MR_inter,
                       AR_rot,
                       iWEC.Madd,
                       iWEC.Crad,
                       iWEC.order.max())
        
        return
        
    def Interaction(self, k0, D, G, TruncOrder):
        """
        Interaction: The interaction theory yields a system of equations to be solved for amplitude scattered
        wave coefficients. This method provides the matrix of the system already inverted for the 
        given spatial distribution of the bodies.

        Args::
            k0 (numpy.ndarray): It contains wave-numbers considered for diffraction and/or
                                radiation problems.
            D (numpy.ndarray): It is the diffraction transfer matrix of the isolated body at
                        different wave-periods. Shape: (period,dim,dim), where dim is 2*Nm+1 and Nm is the order.
            G (numpy.ndarray): It is the force transfer matrix of the isolated body at
                        different wave-periods. Shape: (period,dof,dim), where dof is the total number of degrees
                        of freedom of a single body.
            TruncOrder (int): order of the partial waves used to solve the system of equation

        Returns:
            TS (list): It contains the transformation matrix, T, shrinked down to the truncation order for 
                       diffraction problems. There is a different truncation order for each wave-period so that
                       len(TS) = period.
            TR (list): Same as TS but for radiation problems.
            DS_array (list): It contains the diffraction transfer matrix, D, for the entire array and shrinked down to the truncation order for 
                             diffraction problems. There is a different truncation order for each wave-period so that
                             len(DS) = period.
            DR_array (list): Same as DS_array but for radiation problems.
            GS_array (list): It contains the force transfer matrix, G, for the entire array and shrinked down to the truncation order for 
                             diffraction problems. There is a different truncation order for each wave-period so that
                             len(GS) = period.
            GR_array (list): Same as GS_array but for radiation problems.
            MS_inter (list): It contains the matrix of the linear system of equations, already inverted, to solve 
                             the interaction in diffraction problems. len(MS_inter) = period.
            MR_inter (list): Same as MS_inter but for radiation problems.

        Notes:
            k0 should be in accordance with D and G, i.e. it should be obtained from period (see above).
        """
        dof = len2(G[0,0])
        Nm = (len2(D[0])-1)/2
        # Geometric properties of the array: "
        # Calculating the distances between devices as Dx and Dy "
        Nb= len2(self.coord)
        col=range(Nb)*ones((Nb,Nb))
        row=transpose(col)
        Dx=transpose(self.coord[:,0]*ones((Nb,Nb)))-self.coord[:,0]*ones((Nb,Nb))
        Dy=transpose(self.coord[:,1]*ones((Nb,Nb)))-self.coord[:,1]*ones((Nb,Nb))
        # """ Return the inter-distance between devices, L=[L[0,1],L[0,2],...,L[Nb-1,Nb]] and
        # the azimuth, alf=[alf[0,1],alf[0,2],...,alf[Nb-1,Nb]]
        # and the body-to-body indexes, i.e. [i,j]"""
        L= sqrt(Dx**2+Dy**2)[row<col]
        alf= arctan2(Dy,Dx)[row<col]
        ij0= row[row<col]
        ij1= col[row<col]
        # " Transformation matrix: "
        # " Calculating the transformation matrix "
        dim= 2*Nm+1 # it will be used recurrently
        # """ Calculating Hankel function for all the positive orders and inter-distances
        # H[i,:,:]=[[H[0,k0[i]*L[0]],...,H[0,k0[i]*L[end]]],...,[H[2*Nm,k0[i]*L[0]],...,H[2*Nm,k0[i]*L[end]]]] """
        # " indexing to match jv and yv output to the aforementioned H[i,:,;] "
        matk0,mode,matL= meshgrid(k0,array(range(dim)),L,indexing='ij') # k0 first index, modes are rows and L are columns
        H= jv(mode,matk0*matL)-1j*yv(mode,matk0*matL)
        # """ Assembling the matrix of indexes, ind, to distribuite H throughout Tij[i,:,:]
        # ind= [[0,1,...,dim],[-1,0,...,dim-1],...,[-dim,-dim+1,...,0]].
        # Notice that we will work all the time considering the inner 2D matrix Tij[i,:,:]
        # and for different i we will get Tij for k0[i]. Therefore, we generate an extra
        # matrix of indexes, indp, to distribuite H throughout Tij[:,:,:]. """
        col= range(dim)*ones((dim,dim),dtype=int)
        ind= col-transpose(col)
        indp= array([i*ones((dim,dim)) for i in range(len2(k0))],dtype=int)
        # " From ind take (-1)**ind, since we have just calculated Hankel for positive orders "
        indpp= col-transpose(col) # to prevent dividing by zero in si, which used to be defined si= (ind/abs(ind))**ind
        indpp[indpp == 0]= 1
        si= (indpp/abs(indpp))**ind
        # """ Assembling the transformation matrix
        # T=[[0 T[1,2],...,T[1,Nb]],....[T[Nb,1],...,T[Nb,Nb-1],0]] """
        T = zeros((len2(k0), Nb*dim, Nb*dim), dtype=np.complex64)
                
        length_of_L = len2(L)
                
        for i in xrange(length_of_L): # We will construct T[i,j] and T[j,i] at the same time so we just need len2(L)
            
            row = int(dim * ij0[i])
            col = int(dim * ij1[i])
            
            # " Assembling Tij and Tji "
            Tij = si * H[:,:,i][indp, abs(ind)] * exp(ind * alf[i] * 1j)
            
            row_end = int(row + dim)
            col_end = int(col + dim)
            
            T[:, row:row_end, col:col_end] = Tij[:,:,:] # Tij
            T[:, col:col_end, row:row_end] = Tij[:,:,:] * exp(ind * pi * 1j) # Tji
            
        # """ Transpose for later usage since we finally  decided
        # to work with transposed(T). D and G for the isolated device
        # have been already obtained transposed as well as
        # AP, aS, AR and aR """
        T= transpose(T,(0,2,1))
        # " Interaction: "
        # " Get the matrix of the system of equations inverted and construct D_array and G_array "
        D_array= zeros((len2(k0),Nb*dim,Nb*dim), dtype=np.complex64)
        G_array= zeros((len2(k0),Nb*dim,Nb*dof), dtype=np.complex64) # G is already transposed, numb columns= Nb*dof ok!
        
        for i in xrange(Nb):
            
            start_dim = dim * i
            end_dim = dim * (i + 1)
            
            start_dof = dof * i
            end_dof = dof * (i + 1)
            
            D_array[:, start_dim:end_dim, start_dim:end_dim] = D[:, :, :] # diffraction transfer matrix for the isolated device
            G_array[:, start_dim:end_dim, start_dof:end_dof] = G[:, :, :]
            
        dimp = 2 * TruncOrder + 1
        MS_inter, MR_inter, TS, TR = [], [], [], []
        DS_array, DR_array, GS_array, GR_array = [], [], [], []
        
        length_of_k0 = len2(k0)
            
        for i in xrange(length_of_k0):
            
            rowcolS = (array([range(dimp[i, 0])] * Nb).transpose() +
                       array(range(0, dim * Nb, dim))).transpose().reshape(-1)
            
            rowcolR = (array([range(dimp[i, 1])]*Nb).transpose() +
                       array(range(0, dim * Nb, dim))).transpose().reshape(-1)
            
            rowS, colS = meshgrid(rowcolS, rowcolS, indexing='ij')
            rowR, colR = meshgrid(rowcolR, rowcolR, indexing='ij')
            
            TS.append(T[i, rowS, colS])
            TR.append(T[i, rowR, colR])
        
        # Free memory
        del(T)
        
        for i in xrange(length_of_k0):
            
            rowcolS = (array([range(dimp[i, 0])] * Nb).transpose() +
                       array(range(0, dim * Nb, dim))).transpose().reshape(-1)
            
            rowcolR = (array([range(dimp[i, 1])]*Nb).transpose() +
                       array(range(0, dim * Nb, dim))).transpose().reshape(-1)
            
            rowS, colS = meshgrid(rowcolS, rowcolS, indexing='ij')
            rowR, colR = meshgrid(rowcolR, rowcolR, indexing='ij')
            
            iS = Nm - TruncOrder[i, 0]
            iR = Nm - TruncOrder[i, 1]
            
            DS_array.append(D_array[i, rowS + iS, colS + iS]) 
            DR_array.append(D_array[i, rowR + iR, colR + iR])
        
        # Free memory
        del(D_array)
                
        for i in xrange(length_of_k0):
            
            rowcolS = (array([range(dimp[i, 0])] * Nb).transpose() +
                       array(range(0, dim * Nb, dim))).transpose().reshape(-1)
            
            rowcolR = (array([range(dimp[i, 1])]*Nb).transpose() +
                       array(range(0, dim * Nb, dim))).transpose().reshape(-1)
            
            rowS, colS = meshgrid(rowcolS, rowcolS, indexing='ij')
            rowR, colR = meshgrid(rowcolR, rowcolR, indexing='ij')
            
            iS = Nm - TruncOrder[i, 0]
            iR = Nm - TruncOrder[i, 1]
            
            GS_array.append(G_array[i, rowcolS + iS, :]) 
            GR_array.append(G_array[i, rowcolR + iR, :])
          
        # Free memory
        del(G_array)
            
        for i in xrange(length_of_k0):
            
            rowcolS = (array([range(dimp[i, 0])] * Nb).transpose() +
                       array(range(0, dim * Nb, dim))).transpose().reshape(-1)
            
            rowcolR = (array([range(dimp[i, 1])]*Nb).transpose() +
                       array(range(0, dim * Nb, dim))).transpose().reshape(-1)
            
            rowS, colS = meshgrid(rowcolS, rowcolS, indexing='ij')
            rowR, colR = meshgrid(rowcolR, rowcolR, indexing='ij')

            ms = solve(eye(Nb * dimp[i, 0]) - dot(TS[i], DS_array[i]),
                       eye(Nb * dimp[i, 0]))
            
            mr = solve(eye(Nb * dimp[i, 1]) - dot(TR[i], DR_array[i]),
                       eye(Nb * dimp[i, 1]))
            
            MS_inter.append(ms)
            MR_inter.append(mr)
            
        return (TS,
                TR,
                DS_array,
                DR_array,
                GS_array,
                GR_array,
                MS_inter,
                MR_inter)
        
    def Scattering(self, k0,
                         direction,
                         T,
                         D_array,
                         G_array,
                         M_interaction,
                         BaseOrder):
        """
        Scattering: Computes the excitation force at different wave-periods and directions.

        Args:
            k0 (numpy.ndarray): It contains wave-numbers considered for diffraction problems.
            direction (numpy.ndarray)[rad]: It contains wave-directions considered for
                                diffraction problems.
            T (numpy.ndarray): It is the transformation matrix at different wave-periods. Shape:
                             (period,dim*Nb,dim*Nb), where dim is 2*Nm+1, Nm is the order and Nb the number of bodies.
            D_array (numpy.ndarray): Diffraction transfer matrix for the multibody system. Shape:
                             (period,dim*Nb,dim*Nb).
            G_array (numpy.ndarray): Force transfer matrix for the multibody system. Shape:
                             (period,dof*Nb,dim*Nb), where dof is the total number of degrees of freedom of a single body.
            M_interaction (numpy.ndarray): Inverted matrix of the system of equations to be
                            solved at different wave-periods and directions. Shape: (period,dim*Nb,dim*Nb).
        """
        Nb = len2(self.coord)
        Nm = BaseOrder
        dim = 2 * Nm + 1
        # " Get the ambient planar wave, AP, for the scattering problem for the entire array "
        # " a little bit of indexing to order the results into (k0,dir,bodycoord,mode) "
        matk0, matdir, X, mode = meshgrid(k0,
                                          direction,self.coord[:,0],
                                          array(range(-Nm,Nm+1)),
                                          indexing='ij')
        Y = meshgrid(k0,
                     direction,
                     self.coord[:,1],
                     array(range(-Nm,Nm+1)),
                     indexing='ij')[2]
        
        # " Calculate AP using (k0,dir,bodycoord,mode) "
        AP= exp(-1j*(matk0*(X*cos(matdir)+Y*sin(matdir))+mode*(pi/2+matdir)))
        # """ Arrange AP from 2D to 1D array with respect to (bodycoord,mode)
        # so that [body[0][-Nm:Nm],...,body[Nb-1][-Nm:Nm]] """
        AP= reshape(AP,(len2(k0),len2(direction),(2*Nm+1)*Nb))
        Fex= zeros((len2(k0),len2(direction),len2(G_array[0][0])),
                   dtype=np.complex64)
        # " Save amplitude coefficients "
        if self.cylamplitude:
            self.aS= []
            self.AP= []
        for i in range(len2(k0)):
            # """ Get the scatter cylindric wave for the scattering problem for
            # the entire array (i.e., solve system of equations, i.e. direct matrix method!!!) """
            Nmi= (len2(T[i])/Nb-1)/2
            dimi= 2*Nmi+1
            col= (array([range(dimi)]*Nb).transpose()+\
            array(range(0,dim*Nb,dim))).transpose().reshape(-1)+(Nm-Nmi)
            aS=dot(dot(AP[i][:,col],D_array[i]),M_interaction[i]) # aS= (Id-D*T)\D*AP with M_Interaction=(Id-D*T)**-1. However, all have been transposed for sake of convenience
            # " Save amplitude coefficients "
            if self.cylamplitude:
                self.aS.append(aS)
                self.AP.append(AP[i][:,col])
            # " Calculation of the excitation force for the entire array "
            Fex[i,:,:] = dot(AP[i][:,col]+dot(aS,T[i]),G_array[i]) # Fex= G*aI with aI=AP+T*aS so the overall incident wave for the scattering problem
        
        self.Fex = Fex

    def Radiation(self, k0,
                        freq,
                        T,
                        D_array,
                        G_array,
                        M_interaction,
                        AR_iso,
                        Madd_iso,
                        Crad_iso,
                        BaseOrder):
        """
        Radiation: Computes the radiation forces at different wave-periods for all the degrees of freedom.

        Args:
            k0 (numpy.ndarray): It contains wave-numbers considered for diffraction problems.
            direction (numpy.ndarray)[rad]: It contains wave-directions considered for
                                diffraction problems.
            T (numpy.ndarray): It is the transformation matrix at different wave-periods. Shape:
                             (period,dim*Nb,dim*Nb), where dim is 2*Nm+1, Nm is the order and Nb the number of bodies.
            D_array (numpy.ndarray): Diffraction transfer matrix for the multibody system. Shape:
                             (period,dim*Nb,dim*Nb).
            G_array (numpy.ndarray): Force transfer matrix for the multibody system. Shape:
                             (period,dof*Nb,dim*Nb), where dof is the total number of degrees of freedom of a single body.
            M_interaction (numpy.ndarray): Inverted matrix of the system of equations to be
                            solved at different wave-periods and directions. Shape: (period,dim*Nb,dim*Nb).
            freq (numpy.ndarray)[rad/s]. It contains wave-frequencies considered for
                    radiation problems.
            AR_iso (numpy.ndarray): Amplitude radiated cylindrical wave coefficients at different
                    wave-periods for all the degrees of freedom of the isolated body. Shape: (period,dof,dim).
            Madd_iso (numpy.ndarray): Added mass of the isolated body at different wave-periods.
                        Shape: (period,dof,dof).
            Crad_iso (numpy.ndarray): Radiation damping of the isolated body at different
                    wave-periods. Shape: (period,dof,dof).
        """
        Nb= len2(self.coord)
        Nm= BaseOrder
        dof= len2(Madd_iso[0])
        # " From now on we work for each wave-frequency "
        Madd= zeros((len2(k0),Nb*dof,Nb*dof))
        Crad= zeros((len2(k0),Nb*dof,Nb*dof))
        # " For later usage D, G, T and M_inter are transposed "
        # " Save amplitude coefficients "
        if self.cylamplitude:
            self.aR = []
            self.AR = []
        for k in range(len2(k0)):
            # """ Get the ambient radiated wave, AR, for the radiation problems (dof*Nb problems), i.e
            # (dof_0 induced to all Nb, ..., dof_f induced to all Nb), for the entire array """
            # " cols=[col_0,col_dim,col_2*dim,....col_(Nb-1)*dim,col_1,col_(dim+1),col_(2*dim+1),...] "
            Nmk = (len2(T[k])/Nb-1)/2
            dimk = 2*Nmk+1
            colT = array(range(0,dimk*Nb,dimk)*dimk,dtype=int) +\
                        linspace(0,dimk,dimk*Nb,endpoint=False, dtype=int)
            colAR = (array([range(dimk)])+(Nm-Nmk)).reshape(-1)
            # " Reshape T matrix to get T_re=[T:,0,T:,1,...,T:,Nb-1] so that we append the columns from T:,j "
            T_re = reshape(T[k][colT],(dimk,dimk*Nb**2)) # remember that T is already transposed so we are already getting T columns
            # " AR=[AR[dof_0]*T_re,...,AR[dof_f]*T_re] "
            AR = dot(AR_iso[k][:,colAR],T_re) # each row of aR is constant dof.
            # """ Row_0 of AR will be the ambient radiated wave when bodies drive dof_0, and so on
            # for each Row_i it should be solved Nb subproblems, i.e. for body_0 undergoing dof_i (AR[i,0:dim*Nb])
            # up to body_Nb-1 undergoing dof_i (AR[i,dim*Nb*(Nb-1):dim*Nb*Nb]
            # so it is better to reshape AR as (dof,body which undergoes dof,Nb*dim) where Nb*dim is
            # used to solve aR_(when body j moves in dof_i)= (Id-D*T)\D*AR[i,j,:] """
            AR = reshape(AR,(dof,Nb,Nb*dimk))
            # """ Get the radiated cylindric wave for the radiation problems for
            # the entire array (i.e., solve system of equations, i.e. direct matrix method!!!) """
            Fex_rad = zeros((dof,Nb,dof*Nb), dtype=np.complex64)
            # " Save amplitude coefficients "
            if self.cylamplitude:
                    self.AR.append(AR)
                    aRaux = zeros(AR.shape, dtype=np.complex64)
            for i in range(dof):
                # """ ¤^t means transpose(¤)
                # aR^t=(M_Interaction*D*AR)^t= (AR^t*D^t)*M_interaction^t and AR is already "transposed" since
                # for each dof_i, the columns for each constant row (body_j) are the (2*Nm+1)*Nb ambient
                # wave coefficients required for solving 1 radiation problem (i.e. the radiation problem
                # corresponding to body_j undergoing dof_i) """
                aR = dot(dot(AR[i,:,:],D_array[k]),M_interaction[k])
                # " Save amplitude coefficients "
                if self.cylamplitude:
                    aRaux[i] = aR
                # """ Fex_rad^t= (G*aI)^t=aI^t*G^t and aI^t=AR^t+aR^t*T^t and AR and aR are already "transposed"
                # since for each dof_i, the columns for each constant row (body_j) are the (2*Nm+1)*Nb ambient
                # and scattered wave coefficients  required for solving "excitation force" due to an ambient wave
                # AR[i,j,:], Fex_rad[i,j,k]= excitation force on body_dof (k) due to radiation of body (j) in
                # motion dof (i). body_dof (k) is a global position within the row Fex_rad[i,j,:] so that
                # [0,....,dof,dof+1,...,2*dof+1,2*dof+2,...,3*dof+2,....k=x*dof+y,.....,dof*Nb] and (k)
                # means body (y) in motion dof (x*dof) """
                Fex_rad[i,:,:] = dot(AR[i,:,:] + dot(aR,T[k]), G_array[k])
            # " Save amplitude coefficients "
            if self.cylamplitude:
                self.aR.append(aRaux)
            # " Calculation of the radiation force for the entire array "
            # """ To Fex_rad we should add the hydrodynamics of the device generating the wave that causes
            # Fex_rad[i,j,:] since this just accounts as it was a regular diffraction problem. This is achieved
            # through Madd_iso and Crad_iso. Madd_iso[:,i] and Crad_iso[:,i] gives the radiated force to
            # be added to Fex_rad[i,j,dof*j:dof*(j+1)] for j=range(Nb). """
            Madd_t = Madd_iso[k,:,:].transpose().repeat(Nb,axis=0).reshape(
                                                                (dof,1,dof*Nb))
            Crad_t = Crad_iso[k,:,:].transpose().repeat(Nb,axis=0).reshape(
                                                                (dof,1,dof*Nb)) 
            # """ Crad_t[i,0,:]= [Crad[:,i],...,Crad[:,i]] "Nb times",
            # so column i from Crad is repeated Nb times and distributed along (i,0,:) """
            Ones = eye(Nb,dtype=int).repeat(dof,axis=1) # will be used to distribute Crad_t over eye() type
            Frad_k = -(-freq[k]**2*Ones*Madd_t+1j*freq[k]*Ones*Crad_t)+Fex_rad
            # " Redistribution to a conventional way "
            Frad_k = Frad_k.reshape((Nb*dof,Nb*dof), order= 'F').T
            Madd[k,:,:] = 1/freq[k]**2*real(Frad_k) 
            Crad[k,:,:] = -1/freq[k]*imag(Frad_k)
            
        self.Madd = Madd
        self.Crad = Crad

    def EnergyBalance(self, power_prod_perD_perS):
        """
        EnergyBalance: assess the ratio between outgoing and incoming energy
        Args:
            Nt: unused

        Returns:
            (float): max normalised resource reduction
        """

        # Handy names
        NB = len2(self.B)
        NHs = len2(self.Hs)
        NTp = len2(self.Tp)
        
        # Initialization
        e = zeros((NTp*NHs*NB), dtype = float)
        eP = zeros(e.shape, dtype = float)        
        
        Spec_ = wave_spec(1./self.period,1.,1.)
        
        Spec_.s = self.ScatDiag_spec[2]
        Spec_.gamma = self.ScatDiag_spec[1]
        Spec_.spec_type = self.ScatDiag_spec[0]
        Spec_.add_spectrum()
        self.Spectrum = Spec_
        # loop over sea states 
        ind_SS = -1
        for i0 in range(NB):
            for i1 in range(NHs):
                for i2 in range(NTp):
                    ind_SS +=1
                    
                    eP[ind_SS] = (1025. * 9.81 ** 2 / 64. / np.pi *
                                 self.Hs[i1] ** 2 *
                                 self.Tp[i2] *
                                 np.sqrt((max(self.coord[:, 0]) -
                                          min(self.coord[:, 0])) ** 2 + 
                                         (max(self.coord[:, 1]) -
                                          min(self.coord[:, 1])) ** 2))
                    e[ind_SS] = eP[ind_SS]-power_prod_perD_perS[:,ind_SS].sum()
                    
        self.e = e
        self.eP = eP
                   
        EnergyWarray = e # energy per sea state
        EnergyWOarray = eP # energy per sea state
            
        Balance = EnergyWOarray-EnergyWarray # balance per sea state
        
        return np.max(Balance/EnergyWOarray)
