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
This module contains the method used to read the WECmodes.wp2 file.

.. module:: readWECmodes
   :platform: Windows
   :synopsis: read WECmodes.wp2

.. moduleauthor:: Francesco Ferri <ff@civil.aau.dk>
"""

from glob import glob
from numpy import array, zeros, genfromtxt, exp, transpose, pi, conj, \
real, imag, ones, mean, sqrt
import os
import numpy as np


def readWECmodes(path, f_n="WECmodes.wp2"):
    """
    readWECmodes: reads the information stored in the specified file

    Args:
        path (str): location of the file to read. It needs to be a valid path name

    Optional args:
        f_n (str): name of the file to be read.

    Returns:
        totDOF (float): total number of degree of freedom of the system
        ptoDOF (numpy.ndarray): ID of the dofs connected to the power take off system
        mooDOF (numpy.ndarray): ID of the dofs connected to the mooring system
        freqs (numpy.ndarray): vector of wave frequencies describing the numerical model of the WEC
        direction (numpy.ndarray): vector of wave directions describing the numerical model of the WEC
    """

    with open(glob(os.path.join(path,f_n))[0],'r') as open_file: 
        lines = open_file.readlines()
        
    totDOF = array(lines[1].partition('!')[0],dtype=int)
    ptoDOF = array(lines[2].partition('!')[0].split(),dtype=int)
    mooDOF = array(lines[3].partition('!')[0].split(),dtype=int)
    if len(lines) > 5:  # Nemoh input file
        freqs = array(lines[5].partition('!')[0].split(),dtype=float)*(2.*np.pi)
        Ndir = array(lines[6].partition('!')[0],dtype=int)
        direction = np.linspace(0,2*np.pi,Ndir,endpoint=False)
    else:  # WAMIT case read pot file
        # Open, read pathname.cfg and check IPERIN "
        perORfreq = 1
        dof_gener = 0
        with open(glob(''.join((path,'\*.cfg')))[0],'r') as fcfg:
            for line in fcfg:
                assignment = line.partition('=')
                if assignment[0].split()[0] == 'IPERIN':
                    perORfreq = int(assignment[2])
                elif assignment[0].split()[0] == 'NEWMDS':
                    dof_gener = int(assignment[2])
                    
        # Open, read pathname.pot and get number of wave directions (Ndir) and wave periods (Nper) "
        fpot = open(glob(''.join((path,'\*.pot')))[0],'r').readlines()  # pot file from WAMIT
        wdepth = float(fpot[1].split()[0])
        # Get periods
        Nper= int(fpot[3].split()[0])  # Nper is found in the 4th line (3 for python) and 1st column (0 for python). split() is used for breaking characters into a list with sep=whitespaces
        rowORcolarray = array(fpot[3+1].split(),dtype=float)  # used to know whether they are arranged through a row (one line) or through column
        adjust = 1  # used to know where Ndir starts from
        if perORfreq == 2:
            if len(rowORcolarray) < Nper:
                period = 2*pi/array(fpot[3+1:3+1+Nper],dtype=float)
            elif len(rowORcolarray) == Nper:
                period = 2*pi/rowORcolarray
                adjust = 1./Nper
            else :
                raise IOError('Check on .POT file wave periods')
        else: 
            if len(rowORcolarray) < Nper:
                period = array(fpot[3+1:3+1+Nper],dtype=float)
            elif len(rowORcolarray) == Nper:
                period = rowORcolarray
                adjust = 1./Nper
            else :
                raise IOError('Check on .POT file wave periods')
                
        # Get directions
        start = 3+1+int(Nper*adjust)
        Ndir= int(fpot[start].split()[0])
        rowORcolarray = array(fpot[start+1].split(),dtype=float)
        adjust = 1
        if len(rowORcolarray) < Ndir:
                direction = array(fpot[start+1:start+1+Ndir],dtype=float)
        elif len(rowORcolarray) == Ndir:
            direction = rowORcolarray
            adjust = 1./Ndir
        else :
            raise IOError('Check on .POT file wave directions')
            
        direction *= np.pi/180 
        freqs = 2.*np.pi/period     

    return totDOF, ptoDOF, mooDOF, freqs, direction
