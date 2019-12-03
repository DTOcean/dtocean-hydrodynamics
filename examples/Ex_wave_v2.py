# -*- coding: utf-8 -*-
"""
Example 2 wave case
"""
import os

import numpy as np

from dtocean_hydro import start_logging
from dtocean_hydro.input import WP2_SiteData, WP2_MachineData, WP2input
from dtocean_hydro.main import WP2
from dtocean_hydro.utils.Timer import Timer

# Start the logging system
start_logging()

local = os.path.abspath(os.path.dirname(__file__))

" Lease Area "
leaseAreaVertexUTM = np.array([[0, 0],[1000., 0],[1000., 300.],[0, 300.]],dtype=float)-1.

MinDist = 60
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

FixedArrayLayout = np.array([[0,100.],[100,100.],[200,100.]],'f')
Nogoareas_wave =  []

"Statistical analysis generation"
"--------------------------------"
B=  np.array([0,60,120,180,240,300], "f")/180*np.pi
H=  np.array([1.,1.])
T=  np.array([4.,4.])

#this needs to be checked cos so far the definition was T,H,dirs while it was changed later
p= 0.0/len(B)/len(H)/len(T)* np.ones((len(T),len(H),len(B)))
p[0,0,0] = 1

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
Bathymetry[0,-1] = np.nan
Bathymetry[1,-1] = 10.

BR = np.empty(0)
electrical_connection_point = (-1000.0, -4000.0)

Site = WP2_SiteData(leaseAreaVertexUTM,
                    Nogoareas_wave,
                    Meteocean_wave,
                    None,
                    None,
                    MainDirection,
                    Bathymetry,
                    BR,
                    electrical_connection_point,
                    boundary_padding=10.)


" Machine data "
Te = np.array([3,4,5,6], "f")
Hm0 = np.array([0.5,1.5], "f")
dirs = np.array([0,60,120,180,240,300], "f")/180*np.pi
ScatterDiagram = {'p':np.zeros((4,2,3), "f"),'specType':specType}
PM =  1e4*np.ones((4,2,3), "f")
mooringM =  [np.zeros((4,2,3), "f")]
PTOM =  [1e4*np.ones((4,2,3), "f")]
SingleMachine = {'pf':PM,'Lf2':mooringM,'Lf1':PTOM,'x':Te,'y':Hm0,'z':dirs,'Add':ScatterDiagram}

Type = 'Wave'
lCS = np.array([0, 0., 0.])
Clen = (1,)
YawAngle = 0./180*np.pi

Float_flag = True
InstalDepth = [-np.inf,0]
MinDist = (80,)
OpThreshold = 0.8

#This will trigger the WAMIT or Nemoh calculation for the wave case
#and it is used in the tidal case for check on the user database format
#compability
InternalModel_flag = True

#The input folder defined the place were the software will search for 
#the user tidal database or the WAMIT folder or the WP2Input.wp2 file or the iWP2Pickler.pkl file
InputFolder = os.path.join(local,'inputs_wave')
#InputFolder = r'C:\Users\francesco\Documents\standalonewec\examples\inputs\cylinder\outputs'
#InputFolder = r"C:\Users\francesco\Desktop\dtocean_wave_test"
#DATAfolder = (InternalModel_flag,InputFolder)

#UserArray = {'Option':1,'Value':'rectangular'}
#UserArray = {'Option':1,'Value':'staggered'}
#UserArray = {'Option':1,'Value':'full'}
UserArray = {'Option':2,'Value':FixedArrayLayout}
#UserArray = {'Option':3,'Value':np.random.rand(20,2)*200}

RatedPowerArray = 12
RatedPowerDevice = 6

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
                          wave_data_folder=InputFolder)


" Input assembly "
#Site.printInput(indent=1)
#help(Site)

iWP2input = WP2input(Machine,Site)

if not iWP2input.stopWP2run:
    WPobj = WP2(iWP2input,
                Cfit=np.array([0.]),
                Kfit=np.array([0.]),
                pickup=True,
                debug=True)

    Out = WPobj.optimisationLoop()
    
    if not Out == -1:
        Out.printRes()
    else:
        raise RuntimeError("Error code -1")

