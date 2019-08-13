# -*- coding: utf-8 -*-
"""
Example 1 tidal case
"""

import numpy as np

from dtocean_hydro import start_logging
from dtocean_hydro.input import WP2_SiteData, WP2_MachineData, WP2input
from dtocean_hydro.main import WP2
from dtocean_hydro.utils.Coordinates import *

# Start the logging system
start_logging()

# # Lease Area
leaseAreaVertexUTM = np.array([[0, 0],[1500., 0],[1500., 230.],[0, 230.]],dtype=float)

MinDist = 150
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


FixedArrayLayout = np.array(pos,dtype = float)+0.1



Nogoareas_wave = [np.array([[0, 0],[250., 0],[250., 50.],[0, 50.]],dtype=float),
                  np.array([[0, 0],[250., 20],[550., 450.],[100,600],[0, 50.]],dtype=float)]


# #x,y coordinate of the statistical analysis
# # Statistical analysis generation
# # --------------------------------           
x = np.linspace(-50.,1700.,75)
y = np.linspace(-50,300.,11)    
        
p = np.array([1.])
nx = len(x)
ny = len(y)
N = len(p)

TI = np.array([0.1])
V = 1.*np.ones((ny,nx,N))
U = 2.*np.ones((ny,nx,N))
SSH = 3.0*np.ones((N))
# END of Statistical analysis generation
# --------------------------------------- 
Meteocean = {'V':V,'U':U,'p':p,'TI':TI,'x':x,'y':y,'SSH':SSH} 
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
             ])
        
        
Ct = 0.4*np.ones((100))
y = np.empty((1))
z = np.empty((1))
Bidirection = False
C_IO = np.array([1.,4.])



Type = 'Tidal'
lCS = np.array([0,0,30])
Clen = (20,)
YawAngle = 0./180*np.pi#Need to be clarified
Float_flag = False
InstalDepth = [-np.inf,0]
MinDist = (150,)
OpThreshold = 0
UserArray = {'Option':1,'Value':'staggered'}
#UserArray = {'Option':1,'Value':'full'}
#UserArray = {'Option':3,'Value':np.random.rand(20,2)*200}
RatedPowerArray = 50
RatedPowerDevice = 1

#This will trigger the WAMIT or Nemoh calculation for the wave case
#and it is used in the tidal case for check on the user database format
#compability
InternalModel_flag = False

#The input folder defined the place were the software will search for 
#the user tidal database or the WAMIT folder or the WP2Input.wp2 file or the iWP2Pickler.pkl file
#InputFolder = 'C:\\Users\\fferri\\Desktop\\dtocean_tidal'
#DATAfolder = (InternalModel_flag,InputFolder)

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
                          UserOutputTable,
			  "",
			  tidal_power_curve=Cp,
                          tidal_thrust_curve=Ct,
                          tidal_bidirectional=Bidirection,
                          tidal_cutinout=C_IO,
                          tidal_velocity_curve=X)

" Input assembly "
iWP2input = WP2input(Machine,Site)

if not iWP2input.stopWP2run:
    WPobj = WP2(iWP2input,debug=False)
    #Out = WPobj.optimisationLoop(FixedLayOut=FixedArrayLayout)
    Out = WPobj.optimisationLoop()
    
    if not Out == -1:
        Out.printRes()
    else:
        raise RuntimeError("Error code -1")

