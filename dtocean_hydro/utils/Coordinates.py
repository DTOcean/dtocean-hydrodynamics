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
Created on Mon Feb 23 21:38:48 2015

Lat and Lon are inputed as numpy

.. moduleauthor:: Francesco Ferri <ff@civil.aau.dk>
"""
from math import pi
import numpy as np

def wgs2utm(Lat,Lon,*arg):

    
    if len(arg) > 0:
        if not len(arg) == 2:
            raise IOError('ERROR: Incorrect number of inputs!')
        else:
            utmzone = arg[0]
            utmhemi = arg[1]
    else:
        utmzone = np.empty(0)
        utmhemi = np.empty(0)
           
    
    #if ~sum(double(nargin==[2,4]))
    #    error('Wrong number of input arguments');return
    
    n1=Lat.shape
    n2=Lon.shape
    
    if not (n1==n2):
        print 'ERROR: Lat and Lon should have same size'
    
    
    n3=utmzone.shape
    n4=utmhemi.shape
    if not (n3==n4):
        print 'ERROR: utmzone and utmhemi should have same size'
            
    lat = Lat*pi/180
    lon = Lon*pi/180
    
    # WGS84 parameters
    a = 6378137           #semi-major axis
    b = 6356752.314245    #semi-minor axis
    e = np.sqrt(1-(b/a)**2) #eccentricity
    
    # UTM parameters
    
    if min(n3) > 0:
        Lon0 = 6*utmzone-183 # reference longitude in degrees
    else:
        Lon0 = np.floor(Lon/6)*6+3; # reference longitude in degrees
    
    lon0 = Lon0*pi/180
    k0 = 0.9996               # scale on central meridian
    
    FE = 500000              # false easting
    if min(n4) > 0:
        FN = float(utmhemi=='S')*10000000
    else:
        FN = (Lat < 0)*10000000 # false northing 
    
    # Equations parameters
    eps = e**2/(1-e**2);  # e prime square
    
    # N: radius of curvature of the earth perpendicular to meridian plane
    # Also, distance from point to polar axis
    N = a/np.sqrt(1-e**2*np.sin(lat)**2)
    T = np.tan(lat)**2                
    C = ((e**2)/(1-e**2))*(np.cos(lat))**2
    A = (lon-lon0)*np.cos(lat)                           
    # M: true distance along the central meridian from the equator to lat
    M = a*(  ( 1 - e**2./4 - 3*e**4./64 - 5*e**6./256 )  * lat         \
             -( 3*e**2./8 + 3*e**4./32 + 45*e**6./1024 ) * np.sin(2*lat) \
             +( 15*e**4./256 + 45*e**6./1024 )            * np.sin(4*lat) \
             -(35*e**6./3072 )                             * np.sin(6*lat) )
    
    # easting
    x = FE + k0*N*(                                  A       \
                     + (1-T+C)                      * A**3./6 \
                     + (5-18*T+T**2+72*C-58*eps) * A**5./120 )
                     
    # northing 
    # M(lat0) = 0 so not used in following formula
    y = FN + k0*M + k0*N*np.tan(lat)*( A**2./2  \
                                       + (5-T+9*C+4*C**2)              * A**4./24 \
                                       + (61-58*T+T**2+600*C-330*eps) * A**6./720 )
                                     
    # UTM zone
    if max(n3)>0:
        utmzone = utmzone
        utmhemi = utmhemi
    else:
       utmzone = np.floor(Lon0/6)+31
       utmASCII = ( 83* (Lat < 0) + 78* (Lat >= 0) )
       temp = []
       for ii in range(len(utmASCII)):
           temp.append(chr(utmASCII[ii]))
       utmhemi = np.array(temp)

    return x,y #,utmzone,utmhemi
    
# Pre-logic
def latDD(x):
    temp = []
    for ii in range(len(x)):
        D = int(x[ii][0])
        M = int(x[ii][1])
        S = float(x[ii][2])
        temp.append(D + float(M)/60 + float(S)/3600)
    return np.array(temp)

