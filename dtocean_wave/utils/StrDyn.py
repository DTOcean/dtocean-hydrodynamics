# -*- coding: utf-8 -*-

#    Copyright (C) 2016 Pau Mercadez Ruiz
#    Copyright (C) 2017-2019 Mathew Topper
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
This module provides several methods to solve structure dynamics problems

.. module:: StrDyn
   :platform: Windows
   :synopsis: Structural Dynamic Solver

.. moduleauthor:: Pau Mercadez Ruiz <pmr@civil.aau.dk>
.. moduleauthor:: Mathew Topper <mathew.topper@dataonlygreater.com>
"""

import numpy as np
from numpy.linalg import solve
from scipy.linalg import block_diag

from WatWaves import len2
from spec_class import wave_spec


def MotionFreq(Mass, Damping, Stiffness, Force, freq):
    """
    MotionFreq: Solves the equation of motion in frequency domain for a
    particular frequency and for different forces provided as column-vectors:
        
    Mass*Acceleration+Damping*Velocity+Stiffness*Displacement = Force
    
    where Force=[Force_0,Force_1,...,Force_N]

    Args:
        Mass (numpy.ndarray): Scaling factor for the acceleration.
        Damping (numpy.ndarray): Scaling factor for the velocity.
        Stiffness (numpy.ndarray): Scaling factor for the displacement.
        Force (numpy.ndarray): It contains all the forces (column-vectors).
        freq (float): It is considered a periodic time variation and the
                      periodicity is provided through freq (see in comments).
    Outputs:
        (numpy.ndarray): response amplitude operator for the given frequency

    Notes:
        Considering a periodic time variation such that:
            (displacement=real(X*exp(i*freq*time)))
            
        then the equation of motion in frequency domain is:
            
        X= H\Force, where H=-freq**2*Mass+i*freq*Damping+Stiffness
    """
        
    H = -freq ** 2 * Mass + 1j * freq * Damping + Stiffness
    
    return solve(H, Force)


def scatterdiagram_threshold(p, threshold=0.001):
    """
    scatterdiagram_threshold: apply threshold to the scatter diagram to avoid
    calculation for negligible probability of occurence

    Args:
        p (numpy.ndarray): probability of occurrence of the different sea
                           states

    Optional args:
        threshold (float): cutting level

    Returns:
        p_th (numpy.ndarray): re-scaled probability of occurrence of the
                              different sea states
    """

    p_shape = p.shape

    if not len(p_shape) == 3:
        raise IOError('The number of dimensions of the scatter diagram needs '
                      'to be 3.')

    if not np.allclose(p.sum(), 1.):
        raise IOError('The sum of the probability of occurence of the '
                      'different sea states is not 100%.')

    p_out = p.copy()
    for i0 in range(p_shape[0]):
        for i1 in range(p_shape[1]):
            for i2 in range(p_shape[2]):
                if p_out[i0, i1, i2] < threshold:
                    p_out[i0, i1, i2] = 0.

    return p_out / p_out.sum()


def EnergyProduction(NBodies,
                     Dirs,
                     Hs,
                     Tp,
                     dirs,
                     period,
                     ScatDiag,
                     M,
                     Madd,
                     CPTO,
                     Crad,
                     Kmoor,
                     Khyd,
                     Fex,
                     Kfit,
                     Cfit,
                     RatedPower=None):
    
    """
    EnergyProduction: calculates the energy production for the given sea
    states, based on the given numerical model

    Args:
        NBodies (int):
            number of bodies
        Dirs (numpy.ndarray) [rad]:
            wave directions vector
        Hs (numpy.ndarray) [rad]:
            significant wave heights vector
        Tp (numpy.ndarray) [rad]:
            wave periods vector
        dirs (tuple) [rad]:
            tuple of wave directions associated with the current orientation
        period (numpy.ndarray) [s]:
            wave period used in the definition of the numerical model
        ScatDiag (list):
            list containing the probability of occurrence of the sea states and
            a tuple with the type of wave spectrum (Spec name, gamma,
            directional spreading)
        M (numpy.ndarray):
            Cumulative mass matrix of the isolated WEC
        Madd (numpy.ndarray):
            Cumulative added mass matrix of the array, function of the dof and
            wave frequencies
        Crad (numpy.ndarray):
            Cumulative radiation damping matrix of the array, function of the
            dof and wave frequencies
        CPTO (numpy.ndarray):
            Cumulative PTO damping matrix of the isolated WEC
        Cfit (numpy.ndarray):
            Cumulative fitting damping matrix of the isolated WEC
        Kmoor (numpy.ndarray):
            Cumulative mooring stiffness matrix of the isolated WEC
        Kfit (numpy.ndarray):
            Cumulative fitting stiffness matrix of the isolated WEC
        Khyd (numpy.ndarray):
            Cumulative hydrostatic stiffness matrix of the isolated WEC
        Fex (numpy.ndarray):
            Cumulative excitation force matrix of the array, function of the
            dofs, directions and wave frequencies
        RatedPower (float, optional) [W]:
            Rated power of the device. Any values recorded above this are
            truncated.

    Returns:
        Pyr (numpy.ndarray):
            power production per device per sea states normalised by the
            probability of occurrence of the sea states
        P_dev (numpy.ndarray):
            power production per device per sea states
    
    """
    
    NDir = len2(Dirs)
    NHs = len2(Hs)
    NTp = len2(Tp)
    Np = len2(period)
    ndof = len(Khyd[0, :])
    if RatedPower is None: RatedPower = np.inf
    
    # convert the tuple to numpy.ndarray to ease calculations
    dirs = np.array(dirs)
    prob_occ = scatterdiagram_threshold(ScatDiag[0])
    
    # initialize spectrum
    df = np.abs(1. / period[1:] - 1. / period[:-1])
    fr = 1. / period
    
    Spec_ = wave_spec(fr, 1, 1)
    Spec_.s = ScatDiag[1][2]
    
    # is s=0 or s=30 there is no need for directional spreading.
    if Spec_.s <= 0 or Spec_.s > 30:
        Nd_subset = 1
    else:
        Nd_subset = 3
        
    Spec_.gamma = ScatDiag[1][1]
    Spec_.spec_type = ScatDiag[1][0]
    
    Spec_.add_spectrum()
    
    # initialize output
    P_dev = np.zeros((NBodies, NTp, NHs, NDir), dtype=float)
    Pyr = np.zeros((NBodies, NTp, NHs, NDir), dtype=float)
    
    # time saving for MotionFreq call
    block1 = block_diag(*[M] * NBodies)
    
    for i_Dir in range(NDir):
        
        # define a dirs subset to account for directional spreading
        search_region = range(len(dirs) / Nd_subset)
        
        # pivoting angle index. The angles in the search_region are elements of
        # the original B vector, therefore it is possible to search for a true
        # 0 rather then a minimum difference index
        dir_subset_ind = np.where(
                            np.abs(dirs[search_region] - Dirs[i_Dir]) == 0)[0]
        
        if not dir_subset_ind.size: continue
    
        i_dir = [dir_subset_ind[0] + el * len(search_region)
                                                for el in range(Nd_subset)]
        
        for i_Hs in range(NHs):
            
            for i_Tp in range(NTp):
                
                ## Compute power function (Device dependent)
                powfun = np.zeros((NBodies, Nd_subset, Np), dtype=float)
                
#                if np.isclose(prob_occ[i_Tp, i_Hs, i_Dir], 0): continue
            
                # time saving for MotionFreq call
                block2 = block_diag(*[CPTO[i_Tp, i_Hs, i_Dir] +
                                      Cfit[i_Tp, i_Hs, i_Dir]] * NBodies)
            
                block3 = block_diag(*[Kmoor[i_Tp, i_Hs, i_Dir] +
                                      Khyd +
                                      Kfit[i_Tp, i_Hs, i_Dir]] * NBodies)
            
                for i_fr in range(Np):
                    
                    velo = MotionFreq(block1 + Madd[i_fr],
                                      block2 + Crad[i_fr],
                                      block3,
                                      Fex[i_fr, i_dir].T,
                                      2. * np.pi * fr[i_fr])

                    velo *= 1j * 2. * np.pi * fr[i_fr]
                    
                    for bdy in range(NBodies):
                        
                        for ind in range(Nd_subset):
                            
                            start_idx = ndof * bdy
                            end_idx = ndof * (bdy + 1)
                            
                            Vr = velo[start_idx:end_idx, ind].real
                            Vi = velo[start_idx:end_idx, ind].imag
                            
                            Vrdot = np.dot(CPTO[i_Tp, i_Hs, i_Dir], Vr)
                            Vidot = np.dot(CPTO[i_Tp, i_Hs, i_Dir], Vi)
                            
                            Vrpowfun = 0.5 * np.dot(Vr.T, Vrdot)
                            Vipowfun = 0.5 * np.dot(Vi.T, Vidot)
                            
                            powfun[bdy, ind, i_fr] = Vrpowfun + Vipowfun
                
                ## Compute power matrix
                Spec_.rm_spectrum()
                
                Spec_.Hs = Hs[i_Hs]
                Spec_.fp = 1. / Tp[i_Tp]
                Spec_.t_mean = Dirs[i_Dir]
                Spec_.t = dirs[i_dir]
                
                Spec_.add_spectrum()
                
                # integrate over frequencies
                I = 2. * Spec_.specs[0][2].T * powfun  # a**2 = 2*Spec_val*df
                I = .5 * ((I[:, :, 1:] + I[:, :, :-1]) * df).sum(axis=-1)
                
                # integrate over directions
                if Nd_subset > 1:
                    I = .5 * ((I[:, 1:] + I[:, :-1]) * Spec_.dth).sum(axis=-1)
                
                free_power = I.reshape(-1)
                
                if (free_power > RatedPower).any():
                    free_power = np.clip(free_power, None, RatedPower)
                    
                P_dev[:, i_Tp, i_Hs, i_Dir] = free_power
                
                ## Compute yearly power production (Site dependent)
                powprob = prob_occ[i_Tp, i_Hs, i_Dir]
                powdev = P_dev[:, i_Tp, i_Hs, i_Dir]
                
                Pyr[:, i_Tp, i_Hs, i_Dir] = powprob * powdev
    
    return Pyr, P_dev

# if __name__ == "__main__":
#
#    NBo=1
#    B= np.array([0,180.])/180.*np.pi
#    Hs = np.array([0.5,1,1.5,2.])
#    Tp = np.array([3.5,4.5,5.5,6.5,7.5,8.5,9.5])
#    wdir = np.linspace(0,360,30,endpoint=False)/180.*np.pi
#    period = np.linspace(2,15,50)
#    ScatDiag = (np.ones((7,4,2))/(7*4*2),('Jonswap',3.3,0))
#    M = np.eye(3)
#    Madd = np.array([np.eye(3)]*50)
#    Cpto = np.array([[[np.eye(3)]*2]*4]*7)
#    Crad = np.array([np.eye(3)]*50)
#    Khyd = np.eye(3)
#    Fex = np.ones((50,30,3))
#    Kfit = np.array([[[np.eye(3)]*2]*4]*7)
#    Cfit = np.array([[[np.eye(3)]*2]*4]*7)
#
#
#    Pyr, P_dev = EnergyProduction(NBo, B, Hs, Tp, wdir, period, ScatDiag,
#                                  M, Madd, Cpto, Crad, Cpto, Khyd, Fex,
#                                  Kfit, Cfit)
#
#    print(Pyr)
#    print(P_dev)
