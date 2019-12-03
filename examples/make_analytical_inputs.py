#!/usr/bin/python2.7
# encoding: utf-8
from __future__ import division

import os 
import cPickle as pkl

import numpy as np
import matplotlib.pyplot as plt

def main(input_dir):

    #Gridded data
    #in space
    #xmax = 2000.0
    #xmin = 0.0
    #Dx = 200
    #ymax = 1000
    #ymin = 0.0
    #Dy = 100
    
    #high density
    xmax = 1500.0
    xmin = 0.0
    Dx = 20
    ymax = 230
    ymin = 0.0
    Dy = 20
    
    x = np.linspace(xmin, xmax, (xmax-xmin)//Dx)
    y = np.linspace(ymin, ymax,(ymax-ymin)// Dy)
    X, Y = np.meshgrid(x, y)
    #and time
    tmax = 12.0*60*60
    tmin = 0.0
    Dt = 300
    t = np.linspace(tmin, tmax, Dt)
    period = tmax
    
    umax = 3.1
    vmax = 0.0
    sshmax = 1.0
    
    #time dependent
    #phase = np.pi/2.0
    #U = np.ones((t.shape[0],X.shape[0],X.shape[1]))
    #V = np.ones((t.shape[0],Y.shape[0],Y.shape[1]))
    #SSH  = np.ones((t.shape[0],Y.shape[0],Y.shape[1]))
    #I = 0
    #for i in t:
    #    U[I,:,:] = umax * np.cos(2.0*np.pi*(i/period))
    #    V[I,:,:] = vmax * np.cos(2.0*np.pi*(i/period))
    #    SSH[I,:,:] = sshmax * np.cos(2.0*np.pi*(i/period)+phase)
    #    I += 1 
    
    #time independent
    #V = np.ones(X.shape)
    #for i in range(X.shape[0]):
    #    V[i,:] = vmax * np.cos(2.0*np.pi*(i/(X.shape[0]/2.0)))
    V = vmax * np.ones(X.shape)
    
    #U = np.ones(Y.shape)
    #for i in range(Y.shape[0]):
    #    U[i,:] = umax * np.cos(2.0*np.pi*(i/Y.shape[0]))
    U = umax * np.ones(Y.shape)
    
    SSH  = sshmax * np.ones(Y.shape)
    
    #Plot quiver
#    plt.figure()
#    #Q = plt.quiver(X, Y, U[0,:,:], V[0,:,:], scale=100.0)
#    Q = plt.quiver(X, Y, U, V, scale=100.0)
#    #plt.figure() 
#    #P1 = plt.plot(t, U[:,0,0].T)
#    #P2 = plt.plot(t, SSH[:,0,0].T)
#    plt.show()
    
    ##Default values
    #Bathymetry
    bathy = -60.0 * np.ones(X.shape)
    #Turbulence intensity
    ti = 0.05 * np.ones(X.shape)
    
    #Save as pickle
    data = {}
    data['TI'] = ti
    data['X'] = x
    data['Y'] = y
    data['U'] = U
    data['V'] = V
    data['SSH'] = SSH
    data['bathy'] = bathy
    fpath = os.path.join(input_dir, 'simple_inputs.p')
    f = file(fpath, 'wb')
    pkl.dump(data, f, protocol=pkl.HIGHEST_PROTOCOL)
    f.close()

if __name__ == '__main__':
    
    input_dir = 'inputs_tidal'
    main(input_dir)