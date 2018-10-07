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
This module contains the class used to generate the numerical model of the WEC
using a given WAMIT solution.

.. module:: hyd_WAMIT
   :platform: Windows
   :synopsis: Numerical model of WEC builder

.. moduleauthor:: Pau Mercadez Ruiz <pmr@civil.aau.dk>
"""

import os

import numpy as np

from dtocean_wave.utils.WatWaves import len2, WNumber

import utils.file_utilities as f_util
from utils.transfers import transfers


class WamitReader():
    def __init__(self, data_folder, general_input, get_array_mat=True, debug=False):
        """
        WamitReader: reads output files from WAMIT and buld a frequency domain model of the WEC

        Args:
            pathname (str): pathname for WAMIT output files. E.g. C:\Desktop\WAMIT\case1(.6p, .pot, etc.)

        Returns:
            output (tuple): contains all the output of the module:
                       PhiD, velocity potential for diffraction problem
                       PhiP, velocity potential for income wave
                       PhiR, velocity potential for radiation problem
                       Fex, excitation force frequency function
                       Madd, added mass frequecy function
                       Crad, radiation damping frequency function
                       Khyd, hydrostatic matrix
                       M, mass matrix
                       position,
                       r, cylinder radius
                       t, angular discretisation of the cylinder
                       z, vertical discretisation of the cylinder
                       wdepth, water depth
                       direction, wave directions
                       period, wave periods
                       dof, degree of freedom
                       pto_dof, degree of freedom connected with the pto
                       mooring_dof, degree of freedom connected with the pto

        Notes:
            the folder specified in the pathname nees to contain some key files such as:
                WECmodes.wp2
                .cfg configuration file
                .pot potential file
                .mmx mass matrix file
                .hst hydrostatic matrix file
                .1 radiation result file
                .2 excitation result file
                .fpt multibody WEC decription
                .6p pressure on the circunscribing cylinder file
        """
        self.debug = debug
        self.g = 9.82
        self.rho = 1000.0
        self.data_folder = data_folder
        self.n_dof = general_input['ndof'][0]
        self.pto_dof = general_input['pto_dof']
        self.moor_dof = general_input['mooring_dof']
        self.__get_array_mat = get_array_mat


        self.diffraction_tr_mat = None
        self.force_tr_mat = None
        self.amplitude_coefficient_radiation = None
        self.order = None
        self.truncation_order = None
        self.f_ex = None
        self.m_add = None
        self.c_rad = None
        self.k_hst = None
        self.m_m = None
        self.cyl_radius = None
        self.cyl_z = None
        self.cyl_t = None
        self.water_depth = None
        self.directions = None
        self.periods = None
        self.modes = None


    def load_data(self):
        # read configuration file
        #(n_dof, self.pto_dof, self.moor_dof) = self.__read_gnr()
        self.modes = self.__read_mds()
        (self.__per_or_freq, self.__n_dof_gener) = self.__read_cfg()

        # read potential file
        (self.water_depth,
         self.periods,
         self.directions,
         self.l_cs,
         self.rigid_body_dof) = self.__read_pot()


        n_dof = int(self.rigid_body_dof.sum() + self.__n_dof_gener)
        if not self.n_dof == n_dof:
            raise ValueError("The number of dofs specified in the general input and in the pot file does not match.")
        if [el for el in self.pto_dof if el > n_dof]:
            raise ValueError("The pto dof number is larger than the total number of dofs of the machine.")
        if [el for el in self.moor_dof if el > n_dof]:
            raise ValueError("The mooring dof number is larger than the total number of dofs of the machine.")

        # read radiation problem solution
        (self.m_add, self.c_rad) = self.__read_1()

        # read diffraction problem solution
        self.f_ex = self.__read_2()

        # read mass matrix
        self.m_m = self.__read_mmx()

        # read hydrostatic matrix
        self.k_hst = self.__read_hst()

        # read the field point .fpt
        (x, y, self.cyl_z, self.cyl_radius, self.cyl_t) = self.__read_fpt()

        # read 6p field potential
        (phi_s, phi_r) = self.__read_6p(x, y)

        # calculate diffraction and force transfer matrix
        if self.__get_array_mat:
            (self.diffraction_tr_mat,
            self.force_tr_mat,
            self.amplitude_coefficient_radiation,
            self.order,
            self.truncation_order) = transfers(self.water_depth,
                                      self.directions*np.pi/180.,
                                      self.periods,
                                      (self.cyl_radius, self.cyl_t, self.cyl_z),
                                      phi_s,
                                      phi_r,
                                      self.f_ex)
        else:
            (self.diffraction_tr_mat,
            self.force_tr_mat,
            self.amplitude_coefficient_radiation,
            self.order,
            self.truncation_order) = np.array([0]), np.array([0]), np.array([0]), np.array([0]), np.array([0])
       
        del phi_s, phi_r
        del x, y

    @classmethod
    def __search_file_extension(WamitReader, path, estension):
        file_ls = os.listdir(path)
        filter_files = [el for el in file_ls if el.endswith(estension)]
        if not filter_files:
            return False
        return filter_files

    @classmethod
    def __load_ASCII(WamitReader, path, estension, POINTER=False):
        file_list = WamitReader.__search_file_extension(path, estension)
        if len(file_list) > 1:
            raise ValueError('The given data folder contains multiple entries of the {} file. Please remove the unused ones'.format(estension))
        fid = open(os.path.join(path,file_list[0]),'r') # .xxx file from WAMIT
        if POINTER:  # this option is a lazy reader for the big 6p file
            return fid
        lines = fid.readlines()
        fid.close()

        return lines

    def __read_mds(self, run_fl=True):
        file_name = [el for el in os.listdir(self.data_folder) if el.endswith('.pot')][0]
        f_mds = f_util.load_file(self.data_folder, file_name)[-1]
        modes = f_util.split_string(f_mds, int)
        full_list = [[1, 1, 0, 0, 0, 0, 0],
                     [1, 0, 1, 0, 0, 0, 0],
                     [1, 0, 0, 1, 0, 0, 0],
                     [2, 1, 0, 0, 0, 0, 0],
                     [2, 0, 1, 0, 0, 0, 0],
                     [2, 0, 0, 1, 0, 0, 0]]
            
        return [el for ind, el in enumerate(full_list) if modes[ind]==1]


    def __read_cfg(self):
        if self.debug: print("Read cfg file")
        fcfg = self.__load_ASCII(self.data_folder, '.cfg')
        per_or_freq = 1
        n_dof_gener = 0
        iperin = [int(float(el.split('=')[-1])) for el in fcfg if 'IPERIN' in el]
        newmds = [int(float(el.split('=')[-1])) for el in fcfg if 'NEWMDS' in el]
        if iperin:
            per_or_freq = iperin[0]
        if newmds:
            n_dof_gener = newmds[0]

        return (per_or_freq, n_dof_gener)

    def __read_pot(self):
        if self.debug: print("Read pot file")
        #" Open, read pathname.pot and get number of wave directions (Ndir) and wave periods (Nper) "
        fpot = self.__load_ASCII(self.data_folder, '.pot')
        i_s = 1  # water depth index
        wdepth = float(fpot[i_s])
        i_s += 2  # number of periods index
        n_per = int(float(fpot[i_s].split()[0]))
        i_s += 1  # period indexes start
        if len(fpot[i_s].split()) == n_per:  # data arranged in row
            periods = fpot[i_s].split()
        else:  # data arrange in column
            periods = fpot[i_s:i_s+n_per]
        periods = np.asarray(periods,'f')
        if self.__per_or_freq==2:
            periods = 2*np.pi/periods.copy()
        elif not self.__per_or_freq == 1:
            raise ValueError("The accepted IPERIN option is either 1 or 2.")

        i_s += n_per  # number of direction index
        n_dir = int(float(fpot[i_s].split()[0]))
        i_s += 1  # direction index start
        if len(fpot[i_s].split()) == n_dir:  # data arranged in row
            directions = fpot[i_s].split()
        else:  # data arrange in column
            directions = fpot[i_s:i_s+n_dir]
        directions = np.asarray(directions,'f')
        i_s += n_dir  # number of body index
        # NOTE: the number of bodies need to be 1
        n_bodies = int(float(fpot[i_s]))
        if n_bodies > 1:
            raise ValueError("Wrong number of bodies specified in the pot file.")
        i_s += 2  # coordinate system index
        l_cs = np.asarray(fpot[i_s].split()[:2],'f')
        rigid_body_dof = np.asarray(fpot[i_s+1].split(),'u8')

        return (wdepth, periods, directions, l_cs, rigid_body_dof)

    def __read_1(self):
        if self.debug: print("Read .1 file")
        f1 = self.__load_ASCII(self.data_folder, '.1')[1:]
        n_per = len(self.periods)
        m_add = []
        c_rad = []
        n_lines = len(f1)
        for i, line in enumerate(f1):# skip first line, i.e. header
            linep= line.split()
            m_add.append(float(linep[-2]))
            c_rad.append(float(linep[-1])) # Cradij component is given in the last column
        dofcheck = int(np.sqrt((n_lines+1)/n_per))
        if dofcheck != self.n_dof:
            raise IOError('Total number of dof from Madd does not match total number of dof from dof_solid + dof_gener*NBODY')
        m_add = self.rho*np.asarray(m_add,'f').reshape((n_per,self.n_dof,self.n_dof))  # dimentionalise the wamit output
        c_rad = np.asarray(c_rad, dtype=float).reshape(m_add.shape) # The reshape is easy the way Crad is given by WAMIT (as it was seen for the mass matrix)
        periodp = self.periods.repeat(self.n_dof**2).reshape(c_rad.shape)
        c_rad *= self.rho*2*np.pi/periodp  # dimentionalise the wamit output

        return (m_add, c_rad)

    def __read_2(self):
        if self.debug: print("Read .2 file")
        f2 = self.__load_ASCII(self.data_folder, '.2')[1:]
        n_per = len(self.periods)
        n_dir = len(self.directions)
        fex = np.zeros(n_per*n_dir*self.n_dof,dtype=complex)
        for i,line in enumerate(f2):  # skip first line, i.e. header
            linep = line.split()
            fex[i ]= complex(float(linep[-2]),float(linep[-1]))  # Fex is given as real(fex) and im(fex) in the last two columns
        fex = self.rho*self.g*fex.reshape((n_per,n_dir,self.n_dof))  # The reshape is easy the way Fex is given by WAMIT (as it was seen for the mass matrix)

        return fex

    def __read_mmx(self):
        if self.debug: print("Read .mmx file")
        fmmx = self.__load_ASCII(self.data_folder, '.mmx')[13:]
        mm = np.zeros((self.n_dof,self.n_dof),'f')
        dof_gener = self.__n_dof_gener
        for line in fmmx[:(6+dof_gener)**2]:
            linep = line.split()
            i_in = int(linep[0])-1
            j_in = int(linep[1])-1
            mm[i_in, j_in] = float(linep[2])  # M[row_i,column_j]. M is given [row1,...,rowdof] so then afterwards a reshape will work super!

        return mm

    def __read_hst(self):
        if self.debug: print("Read .hst file")
        fhst = self.__load_ASCII(self.data_folder, '.hst')[1:]
        khst = np.zeros((self.n_dof,self.n_dof),'f')
        for line in fhst:
            linep = line.split()
            i_in = int(linep[0])-1
            j_in = int(linep[1])-1
            khst[i_in, j_in] = float(linep[2])
        khst *= self.rho*self.g

        return khst

    def __read_fpt(self):
        if not self.__get_array_mat:
            return (0,0,0,0,0)
            
        if self.debug: print("Read .fpt file")
        ffpt = self.__load_ASCII(self.data_folder, '.fpt')[1:]  # skip first line, i.e. header
        fpt = []
        for line in ffpt:
            linep = line.split()  # split characters into a list where sep is whitespace
            fpt.append(np.asarray(linep,'f'))  #[linep[0],linep[1],linep[2]),linep[3]])
        fpt = np.asarray(fpt)  # fpt= [ID, x(ID), y(ID), z(ID)]
        r = np.sqrt(fpt[0,1]**2+fpt[0,2]**2)  # radius of the cylinder
        x = fpt[fpt[:,3] == fpt[0,3],1]
        y = fpt[fpt[:,3] == fpt[0,3],2]
        t = np.arctan2(y,x)  # theta of the cylinder
        z = fpt[fpt[:,1] == fpt[0,1],3]  # z of the cylinder

        return (x, y, z, r, t)

    def __read_6p(self, x, y):
        if not self.__get_array_mat:
            nper = len(self.periods)
            ndir = len(self.directions)
            return (1j*np.zeros((nper, ndir,1,1)), 1j*np.zeros((nper, ndir,1,1)))
            
        if self.debug: print("Read .6p file. This will take a while. Take a cup of coffe, sit down and relax!")
        z = self.cyl_z
        t = self.cyl_t
        period = self.periods
        direction = self.directions.copy()
        n_per = len(period)
        n_dir = len(direction)
        dof = self.n_dof
        dof_solid = int(self.rigid_body_dof.sum())
        # Total number of field points
        n_points= len2(z)*len2(t)
        # Due to the size of the 6p file the realines method cannot be used,
        # a lazy method is used instead
        f6p = self.__load_ASCII(self.data_folder, '.6p', POINTER=True)
        burnheader = f6p.readline(); del(burnheader)
        lineR = len2(f6p.readline()) # get the total number of characters within the 1st line
        if dof > dof_solid:
            lineR += len2(f6p.readline()) # if generalized modes are regarded a 2nd line should be read
        inbetween = f6p.read(lineR*(n_points-1)); del(inbetween)
        lineD = len2(f6p.readline())
        f6p.close()

        f6p = self.__load_ASCII(self.data_folder, '.6p', POINTER=True)
        burnheader = f6p.readline(); del(burnheader)
        PhiD = np.zeros(n_per*n_dir*len2(z)*len2(t),dtype=complex)
        PhiR = np.zeros((n_per*dof,len2(z)*len2(t)),dtype=complex)
        for i in range(n_per):
            if self.debug: print("reading period #{} over {}".format(i, n_per))
            # Radiation
            linep = np.asarray(f6p.read(lineR*n_points).split(),
                               'f').reshape((n_points,2+dof*2))

            wfreq = 2*np.pi/period[i]
            dimR = self.g/wfreq**2
            PhiR[dof*i:dof*(i+1),:] = dimR*np.transpose(linep[:,2:][:,range(0,
                                                                           dof*2,2)]+1j*linep[:,3:][:,range(0,dof*2,2)])
            dimD = 1j*self.g/wfreq

            # Diffraction
            for i2 in range(n_dir):
                linep = np.asarray(f6p.read(lineD*n_points).split(),
                                   'f').reshape((n_points,5))

                PhiD[n_points*n_dir*i:n_points*n_dir*(i+1)][n_points*i2:n_points*(i2+1)]= dimD*(linep[:,
                                                                                                -2]+1j*linep[:,-1])
        # the way the nodes of the cylinder are given and the way WAMIT provides PhiD facilitates the reshape
        PhiD = np.reshape(PhiD,(n_per,n_dir,len2(z),len2(t)))
        PhiR = np.reshape(PhiR,(n_per,dof,len2(z),len2(t)))

        # Radiation due to unit amplitude motion
        for f,fr in enumerate(2*np.pi/period) :
            PhiR[f] *= 1j*fr
        direction *= np.pi/180
        wfreq, wdir, fz, fx = np.meshgrid(2*np.pi/period, direction, z, x, indexing='ij', sparse=True)

        # check
        if abs(self.water_depth - max(abs(z))) > 1e-2:
            raise IOError('abs(wdepth - max(abs(z))) > 1e-2')
        wnumber = WNumber(period, self.water_depth)

        wnumber, wdir, fz, fy= np.meshgrid(wnumber, direction, z, y, indexing='ij', sparse=True)

        PhiP = 1j*self.g/wfreq*np.cosh(wnumber*(fz+
                        self.water_depth))/np.cosh(wnumber*self.water_depth)*np.exp(-1j*wnumber*(fx*
                                                                    np.cos(wdir)+fy*np.sin(wdir)))
        f6p.close()

        return (PhiD-PhiP, PhiR)


if __name__ == "__main__":
    #reader = WamitReader(r"C:\Users\francesco\Desktop\Pelamis_input", debug=True)
    #reader.load_data()

    import sys
    sys.path.append(r"C:\Users\francesco\Desktop\test_gui\utils")
    from data_interface import DataStructure
    import dtocean_wave.utils.hdf5_interface as h5i
    data = h5i.load_dict_from_hdf5(r"C:\Users\francesco\Desktop\test_gui\prj\prj_data_collection.hdf5")
    data_path = r"C:\Users\francesco\Desktop\Pelamis_input"
    dataobj = DataStructure(data)


    reader = WamitReader( data_path,
                          dataobj.general_inputs,
                          get_array_mat=False,
                          debug=False)
    #reader.check_inputs()
    reader.load_data()