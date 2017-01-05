#!/usr/bin/python2.7
# encoding: utf-8
"""
This module contains the main classes used to obtain the solution of the hydrodynamic problem for a given wave energy converter

.. module:: output
   :platform: Windows
   :synopsis: wec_external_module module to DTOcean

.. moduleauthor:: Francesco Ferri <ff@civil.aau.dk>
"""
# External Package
import numpy as np


def angle_wrap(angle,conversion=None):
    """
    anglewrap: wrap the angle in the range 0-360. The input and output format are
    decided based on the conversion type

    Args:
        angle (numpy.ndarray): vector containing the angles to be wrapped in the 0-360 range
        conversion (str): decribes the input and output types:
                            Cases:
                                None: maps from degrees to degrees
                                r2r: maps from radiants to radiants
                                r2d: maps from radiants to degrees
                                d2r: maps from degrees to radiants


    Returns:
        wrap (numpy.ndarray): wrapped angle n agreement with the chosen conversion
    """
    def f (x): return (x%360+360)%360
    if not conversion:
        wrap = f(angle)
    elif conversion=='r2r':
        angle *= 180./np.pi
        wrap = f(angle)
        wrap *= np.pi/180.
    elif conversion=='r2d':
        angle *= 180./np.pi
        wrap = f(angle)
    elif conversion=='d2r':
        wrap = f(angle)
        wrap *= np.pi/180.
    else:
        wrap = angle
        raise ValueError('Anglewrap method: Wrong data conversion type!')

    return wrap


def convert_angle(angle):
    """
    convertangle is used to move forth and back between the angle conventions for the wave case.
    This is done to avoid confusion when defining the different coordinate system.
    Hystorically the metocean condition for waves are given in the following convention
    Wave travelling North-South are tagged 0 deg,
                    West-East are tagged 270 deg.
    This convention does not fit with the selected axis East (x) North (y).

    Args:
        angle [rad]:
            angles to be converted from East-North to South-West or viceversa

    Return:
        same object with rotated coordinate system

    """

    return angle_wrap(-(angle+np.pi/2),'r2r')


def convert_te2tp(te, spec_type, gamma):
    coeff = np.array([[  1.22139232e+00],
                            [ -7.26257028e-02],
                            [  1.74397331e-02],
                            [ -2.19288663e-03],
                            [  1.07357912e-04]])
    # convert Te to Tp for the Metocean condition relative to the deployment site
    conversion_factor = 1.16450471
    if spec_type == 'Jonswap':
        if gamma > 7 or gamma < 1:
            print("warning: gamma value outside the range of confidence")
            #module_logger.warning('The gamma value of the JONSWAP spectrum in the metocean data specification is out of the confident range [1-7].')
        conversion_factor = coeff[0] + coeff[1] * gamma + coeff[2] * gamma**2 + coeff[3] * gamma**3 + coeff[4] * gamma**4
    tp = te*conversion_factor
    
    return tp

def convert_tp2te(tp, spec_type, gamma):
    coeff = np.array([[  1.22139232e+00],
                            [ -7.26257028e-02],
                            [  1.74397331e-02],
                            [ -2.19288663e-03],
                            [  1.07357912e-04]])
    # convert Te to Tp for the Metocean condition relative to the deployment site
    conversion_factor = 1.16450471
    if spec_type == 'Jonswap':
        if gamma > 7 or gamma < 1:
            print("warning: gamma value outside the range of confidence")
            #module_logger.warning('The gamma value of the JONSWAP spectrum in the metocean data specification is out of the confident range [1-7].')
        conversion_factor = coeff[0] + coeff[1] * gamma + coeff[2] * gamma**2 + coeff[3] * gamma**3 + coeff[4] * gamma**4
    te = tp/conversion_factor
    
    return te


def set_wdirs_with_yaw(B, s, mA, yA, debug=False):
    """
    set_wdirs_multibody: identifies the orientations of the devices and the relative angles based on the
                        wave directions to be analysed, the directional spreading, and the yaw angle span
                        of the machine

    Args:
        B (nunpy.ndarray): [rad] wave direction vector from scatter diagram
        s (float): [-] spreading parameter from scatter diagram
        mA (float): [rad] main angle
        yA (float): [rad] yawing semi-span

    Optional args:
        debug (boolean): debug flag

    Returns:
        B_new (list): list of orientations and associated angles to be analysed in the array.
                        The associated angles are placed into a tuple that contains at the minimum the orientation itself
                       The output is structured as follow:
                        [ orientation_i, (angles_j)_i ]
    """

    # if not all(x<y for x, y in zip(B, B[1:])):
    #    B.sort()
    if len(set([x for x in B if B.tolist().count(x) > 1])):
        errStr = ('Repeated elements in the wave angle vector')
        raise IOError(errStr)

    # s parameter from metocean
    # assuming s=1 distribution 90° and s=30 distribution 5°, a linear relation is used
    # dir_spr = s*0.0989-0.011636
    dir_spr = -0.0512*s+1.621
    if dir_spr > np.pi/2:
        dir_spr = np.pi
    if dir_spr < np.pi/36-0.001 or s <= 0:
        dir_spr = 0

    # rotate the B vector in agreement with the main angle
    rot_B = angle_wrap(B[:]-mA, 'r2r')
    # identify the set of feasibel/unfeasible orientations
    yaw_ind = np.logical_or(rot_B >= angle_wrap(-(yA+0.0001), 'r2r'), rot_B <= yA)
    yaw_angles = B[yaw_ind]
    notyaw_angles = B[np.logical_not(yaw_ind)]

    if not np.any(yaw_angles):
        tuple_angles = tuple(B)
        if dir_spr:
            tuple_angles += tuple([angle_wrap(el-dir_spr/3.*2, 'r2r') for el in tuple_angles]) + tuple([angle_wrap(el+dir_spr/3.*2, 'r2r') for el in tuple_angles])
        B_new = [mA, tuple_angles]
    else:
        if not notyaw_angles.shape[0] == 0:  # check if the machine can rotate in the whole circle
            # define the min and max angles, those are used as reference for all the angles outised the yaw range
            plusBound = np.argmin(np.abs(angle_wrap(yaw_angles-mA, 'r2r')-yA))
            minusBound = np.argmin(np.abs(angle_wrap(yaw_angles-mA, 'r2r')-angle_wrap(-yA, 'r2r')))

            # distribute the unfeasible orientations between the lower and higher bounds
            ang_distr = angle_wrap(notyaw_angles-(mA+np.pi), 'r2r')
            if yaw_angles.shape[0] == 1:  # one orientation case--> set the unfeasible orientations to the feasible one
                max_range = tuple([angle_wrap(el+mA+np.pi, 'r2r') for el in ang_distr])
            else:
                max_range = tuple([angle_wrap(el+mA+np.pi, 'r2r') for el in ang_distr if el >= np.pi])

            # upper bound
            tuple_angles = (yaw_angles[plusBound],)+max_range
            if dir_spr:
                tuple_angles += tuple([angle_wrap(el-dir_spr/3.*2, 'r2r') for el in tuple_angles]) + tuple([angle_wrap(el+dir_spr/3.*2, 'r2r') for el in tuple_angles])
            B_new = [yaw_angles[plusBound], tuple_angles]

            # lower bound
            if not yaw_angles.shape[0] == 1:
                min_range = tuple([angle_wrap(el+mA+np.pi, 'r2r') for el in ang_distr if el < np.pi])
                tuple_angles = (yaw_angles[minusBound],)+min_range
                if dir_spr:
                    tuple_angles += tuple([angle_wrap(el-dir_spr/3.*2, 'r2r') for el in tuple_angles]) + tuple([angle_wrap(el+dir_spr/3.*2, 'r2r') for el in tuple_angles])
                B_new += [yaw_angles[minusBound],tuple_angles]

            # inside elements
            mask = np.ones(yaw_angles.size, 'bool')
            mask[[plusBound, minusBound]] = False
        else:
            B_new = []
            mask = np.ones(yaw_angles.size, 'bool')

        for elA in yaw_angles[mask]:
            tuple_angles = (elA,)
            if dir_spr:
                tuple_angles += tuple([angle_wrap(el-dir_spr/3.*2,'r2r') for el in tuple_angles]) + tuple([angle_wrap(el+dir_spr/3.*2,'r2r') for el in tuple_angles])
            B_new += [elA, tuple_angles]

    return B_new