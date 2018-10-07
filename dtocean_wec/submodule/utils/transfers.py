# -*- coding: utf-8 -*-

#    Copyright (C) 2016 Pau Mercadez Ruiz
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
This module contains the methods to derive the diffraction and force
transfer matrices, which are required input of the direct matrix
method (Kagemoto and Yue, 1986).

.. module:: transfers
   :platform: Windows
   :synopsis: Numerical model of WEC builder

.. moduleauthor:: Pau Mercadez Ruiz <pmr@civil.aau.dk>
"""

from math import log

import numpy as np
from scipy.special import jv, yv

from dtocean_wave.utils.WatWaves import WNumber


def transfers(water_depth,
              directions,
              periods,
              discrete_cyl,
              vpot_scat,
              vpot_rad,
              fex):
    """ Computes cylindrical amplitude coefficients from the velocity
    potential values on a cylinder.

    :param water_depth: water depth (m)
    :type water_depth: float
    :param directions: wave heading angles (rad)
    :type directions: list or 1D numpy array
    :param periods: wave periods (s)
    :type periods: list or 1D numpy array
    :param discrete_cyl: (radius_cyl, azimuth_cyl, axial_cyl)
    :type discrete_cyl: tuple
    :param radius_cyl: radius (m) of the cylinder
    :type radius_cyl: float
    :param azimuth_cyl: azimuthal discretization (rad) of the cylinder. An equispaced
                        discretization is assumed. Each element of the array is
                        an azimuth. Same axial discretization is assumed for each
                        azimuth
    :type azimuth_cyl: 1D numpy array
    :param axial_cyl: axial discretization (m) of the cylinder. Each element of the array is
                      the z-axial coordinate of a discrete point of the cylinder. Same
                      azimuthal discretization is assumed for each z-coordinate
    :type axial_cyl: 1D numpy array
    :param vpot_scat: complex amplitude of the velocity potential (m**2/s) for each discrete point
                  of the cylinder and for the inputted water depth, wave periods and wave,
                  directions for the diffraction problems. shape (number wave frequencies,
                  number wave directions, number of axial discretization, number of azimuthal
                  discretization)
    :type vpot_scat: 4D numpy array
    :param vpot_rad: complex amplitude of the velocity potential (m**2/s) for each discrete point
                  of the cylinder and for the inputted water depth, wave periods and degrees
                  of freedom, for the diffraction problems. shape (number wave frequencies,
                  number of degrees of freedom, number of axial discretization, number of
                  azimuthal discretization)
    :type vpot_rad: 4D numpy array
    :param fex: complex amplitude of the diffraction excitation force () for the inputted
                water depth, wave periods, wave directions and degrees of freedom.
                shape (number wave frequencies, number of wave directions, number of
                 degrees of freedom)
    :type fex: 3D numpy array
    """
    targ_order = int((len(directions)-1)/2)
    act_order = np.zeros((len(periods), 2), dtype=int)
    decimals = np.zeros((len(periods), 2), dtype=int)
    diffmat = np.zeros((len(periods), 2*targ_order+1, 2*targ_order+1), dtype=complex)
    frcmat = np.zeros((len(periods), 2*targ_order+1, fex.shape[-1]), dtype=complex)
    a_s_rad = np.zeros((len(periods), fex.shape[-1], 2*targ_order+1), dtype=complex)
    dirs, modes = np.meshgrid(directions, range(-targ_order, targ_order+1),
                              indexing='ij', sparse=True)
    for ind, per in enumerate(periods):
        wave_cond = (water_depth, 2.*np.pi/per, WNumber(per, water_depth))
        a_s_scat = bem2cyl(wave_cond, discrete_cyl, vpot_scat[ind], targ_order)
        a_s_rad[ind] = bem2cyl(wave_cond, discrete_cyl, vpot_rad[ind], targ_order)
        a_i_plane = np.exp(-1j*modes*(np.pi/2.+dirs))
        diffmat[ind] = np.linalg.lstsq(a_i_plane, a_s_scat)[0]
        frcmat[ind] = np.linalg.lstsq(a_i_plane, fex[ind])[0]
        # find maximum truncation order
        act_order[ind, 0], decimals[ind, 0] = max_trunc_order(a_s_scat, targ_order, 1e-6)
        act_order[ind, 1], decimals[ind, 1] = max_trunc_order(a_s_rad[ind], targ_order, 1e-6)
    # Shrink G, D and AR according to the truncation order Nm
    ini = targ_order-act_order.max()
    fin = ini+2*act_order.max()+1
    return (diffmat[:, ini:fin, ini:fin].round(decimals.max()),
            frcmat[:, ini:fin, :].round(decimals.max()),
            a_s_rad[:, :, ini:fin].round(decimals.max()),
            act_order.max(axis=0),
            act_order)

def max_trunc_order(a_prob,
                    targ_order,
                    tol):
    """ Selects the biggest wave mode fulfiling a_prob_mode > abs(a_prob).max()*tol

    :param a_prob: cylindrical amplitude coefficients. First dimension is taken for
                   wave directions, if diffraction problem, or for degrees of freedom,
                   if radiation problem.
    :type a_prob: 2D numpy array
    :param targ_order: maximum possible truncation order
    :type targ_order: int
    :param tol: 1e-(number of significant decimals)
    :type tol: float
    """
    tol = max([tol*abs(a_prob).max(), 1e-99])
    decimals = int(abs(log(tol, 10)))
    act_order = 0
    for n_mode in (range(targ_order)):
        if any(abs(a_prob[:, n_mode]) > tol):
            act_order = targ_order-n_mode
            break
    return (act_order, decimals)

def bem2cyl(wave_cond,
            discrete_cyl,
            vpot_cyl,
            trunc_ord):
    """ Computes cylindrical amplitude coefficients from the velocity
    potential values on a cylinder.

    :param wave_cond: (water_depth, cfreq, wnum)
    :type wave_cond: tuple
    :param water_depth: water depth (m)
    :type water_depth: float
    :param cfreq: cyclic wave frequency (rad/s)
    :type cfreq: float
    :param wnum: wave number (rad/m) associated to that wave frequency and water depth
    :type wnum: float
    :param discrete_cyl: (radius_cyl, azimuth_cyl, axial_cyl)
    :type discrete_cyl: tuple
    :param radius_cyl: radius (m) of the cylinder
    :type radius_cyl: float
    :param azimuth_cyl: azimuthal discretization (rad) of the cylinder. An equispaced
                        discretization is assumed. Each element of the array is
                        an azimuth. Same axial discretization is assumed for each
                        azimuth
    :type azimuth_cyl: 1D numpy array
    :param axial_cyl: axial discretization (m) of the cylinder. Each element of the array is
                      the z-axial coordinate of a discrete point of the cylinder. Same
                      azimuthal discretization is assumed for each z-coordinate
    :type axial_cyl: 1D numpy array
    :param vpot_cyl: complex amplitude of the velocity potential (m**2/s) for each discrete point
                     of the cylinder and for the inputted water depth and wave frequency. shape
                     (extra axis, number of axial discretization, number of azimuthal
                     discretization)
    :type vpot_cyl: 3D numpy array
    :param trunc_ord: truncation order for the number of wave modes included.
                      Total number of wave modes is 2*trunc_ord+1
    :type trunc_ord: int
    """
    (water_depth, cfreq, wnum) = wave_cond
    (radius_cyl, azimuth_cyl, axial_cyl) = discrete_cyl
    dz = axial_cyl[1:]-axial_cyl[:-1]
    dth = azimuth_cyl[1]-azimuth_cyl[0] # equispaced is assumed
    rightz = all(dz > 0)
    rightth = dth > 0
    rightth2 = azimuth_cyl[-1] == 2.*np.pi + azimuth_cyl[0] or azimuth_cyl[-1] == azimuth_cyl[0]
    # Build integration domain
    (z_cyl, th_cyl) = np.meshgrid(axial_cyl, azimuth_cyl, indexing='ij')
    # Initialize
    a_s = np.zeros((vpot_cyl.shape[0], 2*trunc_ord+1), dtype=complex)
    for n_mode, mode in enumerate(range(-trunc_ord, trunc_ord+1)):
        integrand = vpot_cyl*np.cosh(wnum*(z_cyl+water_depth))*np.exp(-1j*mode*th_cyl)
        # Integrate along th
        int_th = (integrand[:, :, 1:]+integrand[:, :, :-1]).sum(axis=2)*.5*dth
        if not rightth2: # add last paralepipede
            int_th += (integrand[:, :, 0]+integrand[:, :, -1])*.5*dth
        if not rightth:
            int_th *= -1
        # Integrate I_th along z
        int_th_z = ((int_th[:, 1:]+int_th[:, :-1])*dz).sum(axis=1)*.5
        if not rightz:
            int_th_z *= -1
        # Cm
        cntm = -1j*cfreq/(2*np.pi*9.809)
        cntm *= 2*np.cosh(wnum*water_depth)
        cntm /= water_depth*(1+np.sinh(2*wnum*water_depth)/(2*wnum*water_depth))
        cntm /= jv(mode, wnum*radius_cyl)-1j*yv(mode, wnum*radius_cyl)
        # amplitude coefficients
        a_s[:, n_mode] = cntm*int_th_z
    return a_s

if __name__ == "__main__":
    from nemoh_reader import NemohReader
    run_bem = False
    clear_prj = True
    folder = r'C:\Users\pmr\Documents\Python Scripts\wp2DTOceanExamples\inputs_wave\Cyl'
    prj_folder = folder + r"\Cylinder_out"
    data_folder = folder + r"\Cylinder"
    reader = NemohReader(prj_folder, data_folder, run_bem, clear_prj=clear_prj)
    reader.load_data()
    (diff,
     frc,
     arad,
     trunc_order_max,
     trunc_order) = transfers(reader.water_depth,
                              reader.directions*np.pi/180.,
                              reader.periods,
                              (reader.cyl_r, reader.cyl_t, reader.cyl_z),
                              reader.phi_s,
                              reader.phi_r,
                              reader.f_ex)
