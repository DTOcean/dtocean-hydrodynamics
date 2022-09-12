# -*- coding: utf-8 -*-

#    Copyright (C) 2017-2022 Mathew Topper
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
Created on Tue Sep 25 16:28:33 2018

.. moduleauthor:: Mathew Topper <mathew.topper@dataonlygreater.com>
"""

import os

import pytest
import numpy as np
from scipy.stats import multivariate_normal, norm
from PyQt4 import QtGui

from polite.paths import Directory
from dtocean_wec.main import MainWindow

this_dir = os.path.abspath(os.path.dirname(__file__))


@pytest.fixture
def tidalsite():

    MinDist = 20
    dx = MinDist*4
    dy = dx/4
    pos = []
    for i in range(5):
        for j in range(5):
            if not (i)%2:
                temp = [i*dx,j*dy-dy/2]
            else:
                temp = [i*dx,j*dy]
            
            pos.append(temp)
            
    pos = [item for item in pos if item[1] >= 0]
    
    # #x,y coordinate of the statistical analysis
    # # Statistical analysis generation
    # # --------------------------------           
    x = np.linspace(0., 1000., 100)
    y = np.linspace(0., 500., 50)
    
    # Lease Area
    leaseAreaVertexUTM = np.array([[50., 50.],
                                   [950., 50.],
                                   [950., 500.],
                                   [50., 500.]],dtype=float)
    
    # Nogo areas
    Nogoareas_wave = [] 
            
    nx = len(x)
    ny = len(y)
    
    # Tidal time series
    time_points = 1
    
    rv = norm()
    time_pdf = rv.pdf(np.linspace(-2, 2, time_points))
    time_scaled = time_pdf * (1. / np.amax(time_pdf))
    
    xgrid, ygrid = np.meshgrid(x,y)
    pos = np.dstack((xgrid, ygrid))
    
    rv = multivariate_normal([x.mean(), y.mean()],
                             [[max(x)*5., max(y)*2.],
                              [max(y)*2., max(x)*5.]])
    
    #u_max = 10.
    u_max = 5.
    v_max = 1.
    ssh_max = 1.
    
    grid_pdf = rv.pdf(pos)
    
    u_scaled = grid_pdf * (u_max / np.amax(grid_pdf))
    v_scaled = np.ones((ny, nx)) * v_max
    ssh_scaled = grid_pdf * (ssh_max / np.amax(grid_pdf))
    
    u_arrays = []
    v_arrays = []
    ssh_arrays = []
    
    for multiplier in time_scaled:
        
        u_arrays.append(u_scaled * multiplier)
        v_arrays.append(v_scaled * multiplier)
        ssh_arrays.append(ssh_scaled * multiplier)
    
    U = np.dstack(u_arrays)
    V = np.dstack(v_arrays)
    SSH = np.dstack(ssh_arrays)
    TI = np.array([0.1])
    p = np.ones(U.shape[-1])
    
    # END of Statistical analysis generation
    # --------------------------------------- 
    Meteocean = {'V':V,'U':U,'p':p,'TI':TI,'x':x,'y':y,'SSH':SSH} 
    MainDirection = None#np.array([1.,1.])
    #ang = np.pi*0.25
    #MainDirection = np.array([np.cos(ang),np.sin(ang)])
    
    #Temp check nogo areas
    #xb = np.linspace(0,100,10)
    #yb = np.linspace(0,50,5)
    
    Bathymetry = np.array([-60.])
    
    BR = 1.
    electrical_connection_point = (-1000.0, -4000.0)
    
    out = [leaseAreaVertexUTM,
           Nogoareas_wave,
           Meteocean,
           None,
           None,
           MainDirection,
           Bathymetry,
           BR,
           electrical_connection_point]

    return out


@pytest.fixture
def tidalsite_simple():
    
    MinDist = 20
    dx = MinDist*4
    dy = dx/4
    pos = []
    for i in range(5):
        for j in range(5):
            if not (i)%2:
                temp = [i*dx,j*dy-dy/2]
            else:
                temp = [i*dx,j*dy]
            
            pos.append(temp)
    
    pos = [item for item in pos if item[1] >= 0]
    
    # #x,y coordinate of the statistical analysis
    # # Statistical analysis generation
    # # --------------------------------
    x = np.linspace(0., 1000., 100)
    y = np.linspace(0., 300., 30)
    
    # Lease Area
    leaseAreaVertexUTM = np.array([[50., 50.],
                                   [950., 50.],
                                   [950., 250.],
                                   [50., 250.]],dtype=float)
    
    # Nogo areas
    Nogoareas_wave = [] 
    
    nx = len(x)
    ny = len(y)
    
    # Tidal time series
    time_points = 1
    
    xgrid, ygrid = np.meshgrid(x,y)
    pos = np.dstack((xgrid, ygrid))
    
    u_max = 0.
    v_max = -5.
    ssh_max = 0.
    
    u_scaled = np.ones((ny, nx)) * u_max
    v_scaled = np.ones((ny, nx)) * v_max
    ssh_scaled = np.ones((ny, nx)) * ssh_max
    
    u_arrays = []
    v_arrays = []
    ssh_arrays = []
    
    for _ in xrange(time_points):
        u_arrays.append(u_scaled)
        v_arrays.append(v_scaled)
        ssh_arrays.append(ssh_scaled)
    
    U = np.dstack(u_arrays)
    V = np.dstack(v_arrays)
    SSH = np.dstack(ssh_arrays)
    TI = np.array([0.1])
    p = np.ones(U.shape[-1])
    
    # END of Statistical analysis generation
    # --------------------------------------- 
    Meteocean = {'V':V,'U':U,'p':p,'TI':TI,'x':x,'y':y,'SSH':SSH} 
    MainDirection = None#np.array([1.,1.])
    #ang = np.pi*0.25
    #MainDirection = np.array([np.cos(ang),np.sin(ang)])
    
    #Temp check nogo areas
    #xb = np.linspace(0,100,10)
    #yb = np.linspace(0,50,5)
    
    Bathymetry = np.array([-60.])
    
    BR = 1.
    electrical_connection_point = (-1000.0, -4000.0)
    
    out = [leaseAreaVertexUTM,
           Nogoareas_wave,
           Meteocean,
           None,
           None,
           MainDirection,
           Bathymetry,
           BR,
           electrical_connection_point]

    return out


@pytest.fixture
def wavesite():

    leaseAreaVertexUTM = np.array([[0, 0],
                                   [1000., 0],
                                   [1000., 300.],
                                   [0, 300.]],
                                   dtype=float) - 1.

    Nogoareas_wave =  [np.array([[0, 0],[.1, 0],[.1, .1],[0, .1]])]
    
    "Statistical analysis generation"
    "--------------------------------"
    B=  np.array([0.,270.])/180*np.pi
    H=  np.array([1.])
    T=  np.array([4.])
    
    #this needs to be checked cos so far the definition was T,H,dirs while it
    # was changed later
    p= 1.0/len(B)/len(H)/len(T)* np.ones((len(T),len(H),len(B)))
    
    # Sea state discretization that will be added into the dictionary
    # Meteocean_wave
    #if the frequency are not in increasing order 
    #Nemoh will have a problem in setting up the case to solve
    #if nemoh is used and the the frequency are not uniformaly distributed 
    #there will be a mismatch between the input and the calculated values.
    
    SSH = 0.
    specType = ('Jonswap',3.3,0)
    "END of Statistical analysis generation"
    "---------------------------------------"
    
    Meteocean_wave = {'Te':T,'Hs':H,'B':B,'p':p,'specType':specType,'SSH':SSH}
    MainDirection = None
    
    x = np.linspace(0.,1000.,75)
    y = np.linspace(0.,300.,11) 
    
    X, Y = np.meshgrid(x,y)
    
    Z = -X*0-50.
    
    xyz = np.vstack((X.ravel(),Y.ravel(),Z.ravel())).T
    Bathymetry = xyz
    
    BR = np.empty(0)
    electrical_connection_point = (-1000.0, -4000.0)

    out = [leaseAreaVertexUTM,
           Nogoareas_wave,
           Meteocean_wave,
           None,
           None,
           MainDirection,
           Bathymetry,
           BR,
           electrical_connection_point]
    
    return out


@pytest.fixture
def wavesitebiggamma():

    leaseAreaVertexUTM = np.array([[0, 0],
                                   [1000., 0],
                                   [1000., 300.],
                                   [0, 300.]],
                                   dtype=float) - 1.

    Nogoareas_wave =  [np.array([[0, 0],[.1, 0],[.1, .1],[0, .1]])]
    
    "Statistical analysis generation"
    "--------------------------------"
    B=  np.array([0.,270.])/180*np.pi
    H=  np.array([1.])
    T=  np.array([4.])
    
    #this needs to be checked cos so far the definition was T,H,dirs while it 
    #was changed later
    p= 1.0/len(B)/len(H)/len(T)* np.ones((len(T),len(H),len(B)))
    
    # Sea state discretization that will be added into the dictionary
    # Meteocean_wave
    #if the frequency are not in increasing order 
    #Nemoh will have a problem in setting up the case to solve
    #if nemoh is used and the the frequency are not uniformaly distributed 
    #there will be a mismatch between the input and the calculated values.
    
    SSH = 0.
    specType = ('Jonswap', 10, 0)
    "END of Statistical analysis generation"
    "---------------------------------------"
    
    Meteocean_wave = {'Te':T,'Hs':H,'B':B,'p':p,'specType':specType,'SSH':SSH}
    MainDirection = None
    
    x = np.linspace(0.,1000.,75)
    y = np.linspace(0.,300.,11) 
    
    X, Y = np.meshgrid(x,y)
    
    Z = -X*0-50.
    
    xyz = np.vstack((X.ravel(),Y.ravel(),Z.ravel())).T
    Bathymetry = xyz
    
    BR = np.empty(0)
    electrical_connection_point = (-1000.0, -4000.0)

    out = [leaseAreaVertexUTM,
           Nogoareas_wave,
           Meteocean_wave,
           None,
           None,
           MainDirection,
           Bathymetry,
           BR,
           electrical_connection_point]
    
    return out


@pytest.fixture
def tidal():

    # Machine data
    Type = 'Tidal'
    lCS = np.array([0,0,30])
    Clen = (30,20)
    YawAngle = 0./180*np.pi#Need to be clarified
    Float_flag = False
    InstalDepth = [-np.inf,0]
    MinDist = (100,25)
    OpThreshold = 0
    
    dx = MinDist[0]
    dy = dx
    pos = []
    for i in range(5):
        for j in range(5):
            temp = [i*dx + 70, j*dy + 70]
            pos.append(temp)
            
    pos = [item for item in pos if item[1] >= 0]
    FixedArrayLayout = np.array(pos, dtype=float) + 0.1
    
    UserArray = {'Option': 2,
                 'Value': FixedArrayLayout}

    RatedPowerArray = 50
    RatedPowerDevice = 1

    out = [Type,
           lCS,
           Clen,
           YawAngle,
           Float_flag,
           InstalDepth,
           MinDist,
           OpThreshold,
           UserArray,
           RatedPowerArray,
           RatedPowerDevice]
    
    return out


@pytest.fixture
def tidal_kwargs():
    
    X = np.array([   0.        ,   0.1010101 ,   0.2020202 ,   0.3030303 ,
                     0.4040404 ,   0.50505051,   0.60606061,   0.70707071,
                     0.80808081,   0.90909091,   1.01010101,   1.11111111,
                     1.21212121,   1.31313131,   1.41414141,   1.51515152,
                     1.61616162,   1.71717172,   1.81818182,   1.91919192,
                     2.02020202,   2.12121212,   2.22222222,   2.32323232,
                     2.42424242,   2.52525253,   2.62626263,   2.72727273,
                     2.82828283,   2.92929293,   3.03030303,   3.13131313,
                     3.23232323,   3.33333333,   3.43434343,   3.53535354,
                     3.63636364,   3.73737374,   3.83838384,   3.93939394,
                     4.04040404,   4.14141414,   4.24242424,   4.34343434,
                     4.44444444,   4.54545455,   4.64646465,   4.74747475,
                     4.84848485,   4.94949495,   5.05050505,   5.15151515,
                     5.25252525,   5.35353535,   5.45454545,   5.55555556,
                     5.65656566,   5.75757576,   5.85858586,   5.95959596,
                     6.06060606,   6.16161616,   6.26262626,   6.36363636,
                     6.46464646,   6.56565657,   6.66666667,   6.76767677,
                     6.86868687,   6.96969697,   7.07070707,   7.17171717,
                     7.27272727,   7.37373737,   7.47474747,   7.57575758,
                     7.67676768,   7.77777778,   7.87878788,   7.97979798,
                     8.08080808,   8.18181818,   8.28282828,   8.38383838,
                     8.48484848,   8.58585859,   8.68686869,   8.78787879,
                     8.88888889,   8.98989899,   9.09090909,   9.19191919,
                     9.29292929,   9.39393939,   9.49494949,   9.5959596 ,
                     9.6969697 ,   9.7979798 ,   9.8989899 ,   10.          ])
    
    
    Cp = np.array([ 
                0.        ,  0.        ,  0.        ,  0.        ,  0.        ,
                0.        ,  0.        ,  0.        ,  0.        ,  0.        ,
                0.00248182,  0.0273    ,  0.05211818,  0.07693636,  0.10175455,
                0.12657273,  0.15139091,  0.17620909,  0.20102727,  0.22584545,
                0.25066364,  0.27548182,  0.3003    ,  0.32511818,  0.34993636,
                0.37475455,  0.39957273,  0.42439091,  0.44920909,  0.47402727,
                0.49884545,  0.52366364,  0.54848182,  0.5733    ,  0.59811818,
                0.62293636,  0.64775455,  0.67257273,  0.69739091,  0.72220909,
                0.74702727,  0.77184545,  0.79666364,  0.82148182,  0.8463    ,
                0.86      ,  0.86      ,  0.86      ,  0.86      ,  0.86      ,
                0.86      ,  0.86      ,  0.86      ,  0.86      ,  0.86      ,
                0.86      ,  0.86      ,  0.86      ,  0.86      ,  0.86      ,
                0.86      ,  0.86      ,  0.86      ,  0.86      ,  0.86      ,
                0.86      ,  0.86      ,  0.86      ,  0.86      ,  0.86      ,
                0.86      ,  0.86      ,  0.86      ,  0.86      ,  0.86      ,
                0.86      ,  0.86      ,  0.86      ,  0.86      ,  0.86      ,
                0.86      ,  0.86      ,  0.86      ,  0.86      ,  0.86      ,
                0.86      ,  0.86      ,  0.86      ,  0.86      ,  0.86      ,
                0.86      ,  0.86      ,  0.86      ,  0.86      ,  0.86      ,
                0.86      ,  0.86      ,  0.86      ,  0.86      ,  0.86      
                 ]) * 10.
            
            
    Ct = 0.4 * np.ones((100))
    Bidirection = False
    C_IO = np.array([1., 4.])
    
    out = {"tidal_power_curve": Cp,
           "tidal_thrust_curve": Ct,
           "tidal_bidirectional": Bidirection,
           "tidal_cutinout": C_IO,
           "tidal_velocity_curve": X}
    
    return out


@pytest.fixture
def wave():
    
    Type = 'Wave'
    lCS = np.array([0, 0., 0.])
    Clen = (1,)
    YawAngle = 0./180*np.pi
    
    Float_flag = True
    InstalDepth = [-np.inf,0]
    MinDist = (80,)
    OpThreshold = 0.8
    
    dx = MinDist[0] * 4
    dy = dx/4
    pos = []
    for i in range(3):
        for j in range(3):
            if not (i)%2:
                temp = [i*dx,j*dy-dy/2]
            else:
                temp = [i*dx,j*dy]
            
            pos.append(temp)
            
    pos = [item for item in pos if item[1] >= 0]
    
    FixedArrayLayout = np.array(pos, dtype=float) + 0.1
    
    UserArray = {'Option': 2,
                 'Value': FixedArrayLayout}
    
    RatedPowerArray = 200.
    RatedPowerDevice = 1.
    
    out = [Type,
           lCS,
           Clen,
           YawAngle,
           Float_flag,
           InstalDepth,
           MinDist,
           OpThreshold,
           UserArray,
           RatedPowerArray,
           RatedPowerDevice]
    
    return out

 
@pytest.fixture
def wave_data_folder():
    return os.path.join(this_dir, "..", "examples", 'inputs_wave')


@pytest.fixture
def install_lines():
    
    lines = [
        "[dtocean_tidal]",
        "share_path=dtocean_tidal_mock",
        "",
        "[dtocean_wec]",
        "share_path=dtocean_wec_mock",
        "",
        "[global]",
        "prefix=mock",
        "bin_path=bin_mock"
        ]
    
    return "\n".join(lines)


@pytest.fixture
def main_window(mocker, qtbot, tmpdir, install_lines):
    
    from dtocean_wec.main import QMessageBox
    
    exe_path = tmpdir / "python.exe"
    ini_file = tmpdir / "etc" / "dtocean-data" / "install.ini"
    ini_file.write(install_lines, ensure=True)
    
    mocker.patch('polite.paths.sys.executable', new=str(exe_path))
    mocker.patch('polite.paths.system', new='win32')
    mocker.patch('dtocean_hydro.configure.SiteDataDirectory',
                 return_value=Directory(str(tmpdir)))
    
    mocker.patch.object(QMessageBox,
                        'question',
                        return_value=QtGui.QMessageBox.Yes)
    
    window = MainWindow()
    window.show()
    qtbot.addWidget(window)
    
    return window


@pytest.fixture
def test_data_folder():
    return os.path.join(this_dir, "..", "test_data")
