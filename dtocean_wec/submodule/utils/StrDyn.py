#!/usr/bin/python2.7
# encoding: utf-8
"""
This module provides several tools for structure dynamics problems

.. module:: input
   :platform: Windows
   :synopsis: Input module to DTOcean WP2

.. moduleauthor:: Pau Mercadez Ruiz <pmr@civil.aau.dk>
"""
from numpy.linalg import solve
from scipy.linalg import block_diag
import numpy as np
from WatWaves import len2
from spec_class import wave_spec


def MotionFreq(Mass, Damping, Stiffness, Force, freq):
    """
    MotionFreq: Solves the equation of motion in frequency domain for a particular
    frequency and for different forces provided as column-vectors:
    Mass*Acceleration+Damping*Velocity+Stiffness*Displacement=Force=[Force_0,Force_1,...,Force_N]

    Args:
        Mass (numpy.ndarray): Scaling factor for the acceleration.
        Damping (numpy.ndarray): Scaling factor for the velocity.
        Stiffness (numpy.ndarray): Scaling factor for the displacement.
        Force (numpy.ndarray): It contains all the forces (column-vectors).
        freq (float): It is considered a periodic time variation and the periodicity
                is provided through freq (see in comments).
    Outputs:
        (numpy.ndarray): response amplitude operator for the given frequency

    Notes:
        Considering a periodic time variation (displacement=real(X*exp(i*freq*time)))
        the equation of motion in frequency domain turns out:
        X= H\Force, where H=-freq**2*Mass+i*freq*Damping+Stiffness """
    H = -freq**2*Mass+1j*freq*Damping+Stiffness
    return solve(H, Force)


def scatterdiagram_threshold(p, threshold=0.001):
    """
    scatterdiagram_threshold: apply threshold to the scatter diagram to avoid calculation for negligible probability of occurence

    Args:
        p (numpy.ndarray): probability of occurrence of the different sea states

    Optional args:
        threshold (float): cutting level

    Returns:
        p_th (numpy.ndarray): re-scaled probability of occurrence of the different sea states
    """

    p_shape = p.shape

    if not len(p_shape) == 3:
        raise IOError('The number of dimensions of the scatter diagram needs to be 3.')

    if not np.allclose(p.sum(), 1.):
        raise IOError('The sum of the probability of occurence of the different sea states'
                      '\n is not 100%.')

    p_out = p.copy()
    for i0 in range(p_shape[0]):
        for i1 in range(p_shape[1]):
            for i2 in range(p_shape[2]):
                if p_out[i0, i1, i2] < threshold:
                    p_out[i0, i1, i2] = 0.

    return p_out/p_out.sum()


def EnergyProduction(NBodies, Dirs, Hs, Tp, dirs, period, ScatDiag,
                      M, Madd, CPTO, Crad, Kmoor, Khyd, Fex,
                      Kfit, Cfit):
    """
    EnergyProduction: calculates the energy production for the given sea states, based on the given numerical model

    Args:
        NBodies (int): number of bodies
        Dirs (numpy.ndarray) [rad]: wave directions vector
        Hs (numpy.ndarray) [rad]: significant wave heights vector
        Tp (numpy.ndarray) [rad]: wave periods vector
        dirs (tuple) [rad]: tuple of wave directions associated with the current orientation
        period (numpy.ndarray) [s]: wave period used in the definition of the numerical model
        ScatDiag (list): list containing the probability of occurrence of the sea states and a tuple with the type of wave
                            spectrum (Spec name, gamma, directional spreading)
        M (numpy.ndarray): Cumulative mass matrix of the isolated WEC
        Madd (numpy.ndarray): Cumulative added mass matrix of the array, function of the dof and wave frequencies
        Crad (numpy.ndarray): Cumulative radiation damping matrix of the array, function of the dof and wave frequencies
        CPTO (numpy.ndarray): Cumulative PTO damping matrix of the isolated WEC
        Cfit (numpy.ndarray): Cumulative fitting damping matrix of the isolated WEC
        Kmoor (numpy.ndarray): Cumulative mooring stiffness matrix of the isolated WEC
        Kfit (numpy.ndarray): Cumulative fitting stiffness matrix of the isolated WEC
        Khyd (numpy.ndarray): Cumulative hydrostatic stiffness matrix of the isolated WEC
        Fex (numpy.ndarray): Cumulative excitation force matrix of the array, function of the dofs, directions and wave frequencies

    Returns:
        Pyr (numpy.ndarray): power production per device per sea states normalised by the probability of occurrence of the
                                sea states
        P_dev (numpy.ndarray): power production per device per sea states
    """
    NDir, NHs, NTp, Np, ndof = len2(Dirs), len2(Hs), len2(Tp), len2(period), len(Khyd[0,:])
    dirs = np.array(dirs)  #convert the tuple to numpy.ndarray to ease calculations
    prob_occ = scatterdiagram_threshold(ScatDiag[0])
    # initialize spectrum
    df = np.abs(1./period[1:]-1./period[:-1])
    fr = 1/period
    Spec_ = wave_spec(fr, 1, 1)
    Spec_.s = ScatDiag[1][2]
    if Spec_.s <= 0 or Spec_.s > 30:  # is s=0 or s=30 there is no need for directional spreading.
        Nd_subset = 1
    else:
        Nd_subset = 3
    Spec_.gamma = ScatDiag[1][1]
    Spec_.spec_type = ScatDiag[1][0]
    Spec_.add_spectrum()
    # initialize output
    P_dev = np.zeros((NBodies, NTp, NHs, NDir), dtype=float)
    Pyr = np.zeros((NBodies, NTp, NHs, NDir), dtype=float)
    for i_Dir in range(NDir):
        search_region = range(len(dirs)/Nd_subset) # define a dirs subset to account for directional spreading
        dir_subset_ind = np.where(np.abs(dirs[search_region]-Dirs[i_Dir])==0)[0] # pivoting angle index. The angles in the search_region are elements of the original B vector, therefore it is possible to search for a true 0 rather then a minimum difference index
        if not dir_subset_ind.size: continue
        i_dir = [dir_subset_ind[0]+el*len(search_region) for el in range(Nd_subset)]
        for i_Hs in range(NHs):
            for i_Tp in range(NTp):
                ## Compute power function (Device dependent)
                powfun = np.zeros((NBodies, Nd_subset, Np), dtype=float)
                if prob_occ[i_Tp, i_Hs, i_Dir] == 0.: continue
                for i_fr in range(Np):
                    velo = MotionFreq(block_diag(*[M]*NBodies) + Madd[i_fr],
                                     block_diag(*[CPTO[i_Tp, i_Hs, i_Dir]+Cfit[i_Tp, i_Hs, i_Dir]]*NBodies) + Crad[i_fr],
                                     block_diag(*[Kmoor[i_Tp, i_Hs, i_Dir] + Khyd + Kfit[i_Tp, i_Hs, i_Dir]]*NBodies),
                                     Fex[i_fr, i_dir].T,
                                     2.*np.pi*fr[i_fr])

                    velo *= 1j*2.*np.pi*fr[i_fr]
                    for bdy in range(NBodies):
                        for ind in range(Nd_subset):
                            Vr = velo[ndof*bdy:ndof*(bdy+1), ind].real
                            Vi = velo[ndof*bdy:ndof*(bdy+1), ind].imag
                            powfun[bdy, ind, i_fr] = 0.5*np.dot(Vr.T, np.dot(CPTO[i_Tp, i_Hs, i_Dir], Vr))
                            powfun[bdy, ind, i_fr] += 0.5*np.dot(Vi.T, np.dot(CPTO[i_Tp, i_Hs, i_Dir], Vi))
                ## Compute power matrix
                Spec_.rm_spectrum()
                Spec_.Hs, Spec_.fp, Spec_.t_mean, Spec_.t = Hs[i_Hs], 1/Tp[i_Tp], Dirs[i_Dir], dirs[i_dir]
                Spec_.add_spectrum()
                I = 2*Spec_.specs[0][2].T*powfun # a**2 = 2*Spec_val*df
                I = .5*((I[:,:,1:]+I[:,:,:-1])*df).sum(axis=-1) # integrate over frequencies
                if Nd_subset > 1 :
                    I = .5*((I[:,1:]+I[:,:-1])*Spec_.dth).sum(axis=-1) # integrate over directions
                P_dev[:, i_Tp, i_Hs, i_Dir] = I.reshape(-1)
                ## Compute yearly power production (Site dependent)
                Pyr[:, i_Tp, i_Hs, i_Dir] = prob_occ[i_Tp, i_Hs, i_Dir]*P_dev[:, i_Tp, i_Hs, i_Dir]
    return Pyr, P_dev

if __name__ == "__main__":

    NBo=1
    B= np.array([0,180.])/180.*np.pi
    Hs = np.array([0.5,1,1.5,2.])
    Tp = np.array([3.5,4.5,5.5,6.5,7.5,8.5,9.5])
    wdir = np.linspace(0,360,30,endpoint=False)/180.*np.pi
    period = np.linspace(2,15,50)
    ScatDiag = (np.ones((7,4,2))/(7*4*2),('Jonswap',3.3,0))
    M = np.eye(3)
    Madd = np.array([np.eye(3)]*50)
    Cpto = np.array([[[np.eye(3)]*2]*4]*7)
    Crad = np.array([np.eye(3)]*50)
    Khyd = np.eye(3)
    Fex = np.ones((50,30,3))
    Kfit = np.array([[[np.eye(3)]*2]*4]*7)
    Cfit = np.array([[[np.eye(3)]*2]*4]*7)


    Pyr, P_dev = EnergyProduction(NBo, B, Hs, Tp, wdir, period, ScatDiag,
                                                M, Madd, Cpto, Crad, Cpto, Khyd, Fex,
                                                Kfit, Cfit)

    print(Pyr)
    print(P_dev)
