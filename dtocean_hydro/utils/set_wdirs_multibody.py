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
This module contains the function used to set the wave direction by taking into account the
possible self orientation feature of the device

.. module:: set_wdirs_multibody
   :platform: Windows
   :synopsis: Defines wave direction for the array layout

.. moduleauthor:: Francesco Ferri <ff@civil.aau.dk>
.. moduleauthor:: Mathew Topper <mathew.topper@dataonlygreater.com>
"""

import numpy as np

def anglewrap(angle,conversion=None):
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

def convertangle(angle):
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

    return anglewrap(90 - angle, 'd2r')


def set_wdirs_multibody(B, s, mA, yA, debug=False):
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
    rot_B = anglewrap(B[:]-mA, 'r2r')
    # identify the set of feasibel/unfeasible orientations
    yaw_ind = np.logical_or(rot_B >= anglewrap(-(yA+0.0001), 'r2r'), rot_B <= yA)
    yaw_angles = B[yaw_ind]
    notyaw_angles = B[np.logical_not(yaw_ind)]

    if not np.any(yaw_angles):
        tuple_angles = tuple(B)
        if dir_spr:
            tuple_angles += tuple([anglewrap(el-dir_spr/3.*2, 'r2r') for el in tuple_angles]) + tuple([anglewrap(el+dir_spr/3.*2, 'r2r') for el in tuple_angles])
        B_new = [mA, tuple_angles]
    else:
        if not notyaw_angles.shape[0] == 0:  # check if the machine can rotate in the whole circle
            # define the min and max angles, those are used as reference for all the angles outised the yaw range
            plusBound = np.argmin(np.abs(anglewrap(yaw_angles-mA, 'r2r')-yA))
            minusBound = np.argmin(np.abs(anglewrap(yaw_angles-mA, 'r2r')-anglewrap(-yA, 'r2r')))

            # distribute the unfeasible orientations between the lower and higher bounds
            ang_distr = anglewrap(notyaw_angles-(mA+np.pi), 'r2r')
            if yaw_angles.shape[0] == 1:  # one orientation case--> set the unfeasible orientations to the feasible one
                max_range = tuple([anglewrap(el+mA+np.pi, 'r2r') for el in ang_distr])
            else:
                max_range = tuple([anglewrap(el+mA+np.pi, 'r2r') for el in ang_distr if el >= np.pi])

            # upper bound
            tuple_angles = (yaw_angles[plusBound],)+max_range
            if dir_spr:
                tuple_angles += tuple([anglewrap(el-dir_spr/3.*2, 'r2r') for el in tuple_angles]) + tuple([anglewrap(el+dir_spr/3.*2, 'r2r') for el in tuple_angles])
            B_new = [yaw_angles[plusBound], tuple_angles]

            # lower bound
            if not yaw_angles.shape[0] == 1:
                min_range = tuple([anglewrap(el+mA+np.pi, 'r2r') for el in ang_distr if el < np.pi])
                tuple_angles = (yaw_angles[minusBound],)+min_range
                if dir_spr:
                    tuple_angles += tuple([anglewrap(el-dir_spr/3.*2, 'r2r') for el in tuple_angles]) + tuple([anglewrap(el+dir_spr/3.*2, 'r2r') for el in tuple_angles])
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
                tuple_angles += tuple([anglewrap(el-dir_spr/3.*2,'r2r') for el in tuple_angles]) + tuple([anglewrap(el+dir_spr/3.*2,'r2r') for el in tuple_angles])
            B_new += [elA, tuple_angles]

    if debug:
        import matplotlib as mpl
        import matplotlib.pyplot as plt
        import matplotlib.colors as colors
        import matplotlib
        from matplotlib.patches import Wedge
        from matplotlib.collections import PatchCollection

        col = colors.cnames
        f = plt.figure()
        ax = f.add_subplot(111)

        MA_x = [0, np.cos((mA))]
        MA_y = [0, np.sin((mA))]
        mA_x = [0, np.cos((mA+np.pi))]
        mA_y = [0, np.sin((mA+np.pi))]
        yA_x = [0, np.cos((mA-yA))]
        yA_y = [0, np.sin((mA-yA))]
        YA_x = [0, np.cos((mA+yA))]
        YA_y = [0, np.sin(mA+yA)]

        # add a wedge for the angle span
        patches = []
        t = mpl.transforms.Affine2D().rotate_deg(mA*180/np.pi) + ax.transData

        # add a wedge for the upper angle range
        patches += [Wedge([0, 0], 0.5, -yA*180/np.pi,
                               yA*180/np.pi, ec="none",
                                alpha=0.3)]

        if np.any(yaw_angles) and not notyaw_angles.shape[0] == 0:
            if yaw_angles.shape[0] == 1:
                patches +=[Wedge((0, 0), 1.2, (yaw_angles[plusBound]-mA)*180/np.pi,
                               (np.pi)*180/np.pi, width=0.10)]
            else:
                patches +=[Wedge((0, 0), 1.2, (yaw_angles[plusBound]-mA)*180/np.pi,
                               (np.pi)*180/np.pi, width=0.10),
                         Wedge((0, 0), 1.2, (np.pi)*180/np.pi,
                               (yaw_angles[minusBound]-mA)*180/np.pi, width=0.10)]

        colors = 100*np.random.rand(len(patches))
        p = PatchCollection(patches, cmap=matplotlib.cm.jet, alpha=0.4)
        p.set_array(np.array(colors))
        p.set_transform(t)
        ax.add_collection(p)

        ax.plot(mA_x, mA_y, '-.', color='#999999')
        ax.plot(MA_x, MA_y, color='#999999')
        ax.plot(yA_x, yA_y, '--', color='#aaaaaa')
        ax.plot(YA_x, YA_y, '--', color='#aaaaaa')

        for ang in B.copy():
            ax.plot([0.9*np.cos(ang),
                     1.*np.cos(ang)],
                    [0.9*np.sin(ang),
                     1.*np.sin(ang)],'k',lw=5)

        for i, e in enumerate(B_new):
            if not i%2:
                ax.plot([0.8*np.cos(e),1.1*np.cos(e)],
                         [0.8*np.sin(e),1.1*np.sin(e)],
                          color=col[col.keys()[i]],lw=3)
            else:
                for eA in e:
                    ax.plot([0.6*np.cos(eA),0.7*np.cos(eA)],
                             [0.6*np.sin(eA),0.7*np.sin(eA)],
                              color=col[col.keys()[i-1]],lw=2)

        plt.axis('equal')
        plt.show()

    return B_new

if __name__ == "__main__":
    mA = 50.*np.pi/180
    yA = 100*np.pi/180
    B = np.linspace(0, 360, 8, endpoint=False)
    B *= np.pi/180

    for s in [28]:
        a = set_wdirs_multibody(B, s, mA, yA, debug=True)

