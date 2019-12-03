# -*- coding: utf-8 -*-
"""
Example 1 tidal case
"""

import numpy as np

from scipy.stats import multivariate_normal, norm

from dtocean_hydro import start_logging
from dtocean_hydro.input import WP2_SiteData, WP2_MachineData, WP2input
from dtocean_hydro.main import WP2
from dtocean_hydro.utils.Coordinates import *

from convert import make_tide_statistics


# Start the logging system
start_logging()



FixedArrayLayout = np.c_[np.linspace(940, 60, 9), np.asarray([150]*9)]


# #x,y coordinate of the statistical analysis
# # Statistical analysis generation
# # --------------------------------           
x = np.linspace(0., 1000., 100)
y = np.linspace(0., 300., 30)

# Lease Area
leaseAreaVertexUTM = np.array([[50., 50.],[950., 50.],[950., 250.],[50., 250.]],dtype=float)

# Nogo areas
Nogoareas_wave = [] #[np.array([[50., 50.],[200., 50.],[200., 100.],[50., 100.]],dtype=float)]
        
nx = len(x)
ny = len(y)

# Tidal time series
n_bins = 6
time_points = 48
t = np.linspace(0, 1, time_points)

rv = norm()
time_sin = np.sin(np.linspace(0, 4*np.pi, time_points))
time_scaled = time_sin * (1. / np.amax(time_sin))

xgrid, ygrid = np.meshgrid(x,y)
pos = np.dstack((xgrid, ygrid))
pos = np.swapaxes(pos, 0, 1)

rv = multivariate_normal([x.mean(), y.mean()],
                         [[max(x)*5., max(y)*2.],
                          [max(y)*2., max(x)*5.]])

u_max = 0.
v_max = 6.
ssh_max = 1.
TI = 0.1

grid_pdf = rv.pdf(pos)

u_scaled = np.ones((nx, ny)) * u_max
v_scaled = np.ones((nx, ny)) * v_max
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
TI = np.ones(SSH.shape) * TI

xc = x[int(nx/2)]
yc = y[int(ny/2)]

tide_dict = {"U": U, 
             "V": V, 
             "SSH": SSH, 
             "TI": TI, 
             "x": x, 
             "y": y,
             "t": t,  
             "xc": xc, 
             "yc": yc,
             "ns": n_bins
             }

# END of Statistical analysis generation
# --------------------------------------- 
Meteocean = make_tide_statistics(tide_dict)
VelocityShear = np.array([7.])
MainDirection = None#np.array([1.,1.])
#ang = np.pi*0.25
#MainDirection = np.array([np.cos(ang),np.sin(ang)])

#Temp check nogo areas
#xb = np.linspace(0,100,10)
#yb = np.linspace(0,50,5)

X, Y = np.meshgrid(x,y)

Z = (-X*0.1-1)*0-60

xyz = np.vstack((X.ravel(),Y.ravel(),Z.ravel())).T
Bathymetry = np.array([-60.])

BR = 1.
electrical_connection_point = (-1000.0, -4000.0)

Site = WP2_SiteData(leaseAreaVertexUTM,
                    Nogoareas_wave,
                    Meteocean,
                    None,
                    None,
                    MainDirection,
                    Bathymetry,
                    Geophysics,
                    BR,
                    electrical_connection_point)         


# Machine data
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


Cp = np.array([ 0.        ,  0.        ,  0.        ,  0.        ,  0.        ,
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
             ])*10.
        
        
Ct = 0.4*np.ones((100))
y = np.empty((1))
z = np.empty((1))
Bidirection = False
C_IO = np.array([1.,4.])

Type = 'Tidal'
lCS = np.array([0,0,30])
Clen = (8,)
YawAngle = 0.
Float_flag = False
InstalDepth = [-np.inf,0]
MinDist = (10,)
OpThreshold = 0
#UserArray = {'Option':1,'Value':'rectangular'}
UserArray = {'Option':2,'Value':FixedArrayLayout}
#UserArray = {'Option':1,'Value':'full'}
#UserArray = {'Option':3,'Value':np.random.rand(20,2)*200}
RatedPowerArray = 20
RatedPowerDevice = 1



UserOutputTable = None

Machine = WP2_MachineData(Type,
                          lCS,
                          Clen,
                          YawAngle,
                          Float_flag,
                          InstalDepth,
                          MinDist,
                          OpThreshold,
                          UserArray,
                          RatedPowerArray,
                          RatedPowerDevice,
                          tidal_power_curve=Cp,
                          tidal_thrust_curve=Ct,
                          tidal_bidirectional=Bidirection,
                          tidal_cutinout=C_IO,
                          tidal_velocity_curve=X)

" Input assembly "
iWP2input = WP2input(Machine,Site)

if not iWP2input.stopWP2run:
    WPobj = WP2(iWP2input,debug=True)
    #Out = WPobj.optimisationLoop(FixedLayOut=FixedArrayLayout)
    Out = WPobj.optimisationLoop()
    
    if not Out == -1:
        Out.printRes()
    else:
        raise RuntimeError("Error code -1")

