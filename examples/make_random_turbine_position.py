#!/usr/bin/python2.7
# encoding: utf-8

import os
import cPickle as pkl

import numpy as np

#Global variables
#from global_variables import *
#Global variables

def main(input_dir, InpName='simple_inputs.p'):

    #Load data from pickle file
    fpath = os.path.join(input_dir, InpName)
    f = file(fpath, 'rb')
    data = pkl.load(f)
    f.close()
    
    #Turbine array features
    nb = 50
    #randIndX = np.random.random_integers(0.0, data['X'].shape[1]-1, size=10)
    #randIndY = np.random.random_integers(0.0, data['Y'].shape[0]-1, size=10)
    #Y = data['Y'][randIndY,0]
    #X = data['X'][0,randIndX]
    #ind_turb = np.asarray((randIndX, randIndY)).T
    #Z = np.ones(X.shape) * -30.0
    #xyz_turb = np.asarray((X, Y, Z)).T
    #data = {}
    #data['position'] = xyz_turb
    #data['index'] = ind_turb
    
    Data = {}
    xmin = data['X'].min().min()
    xmax = data['X'].max().max()
    ymin = data['Y'].min().min()
    ymax = data['Y'].max().max()
    
    #random positions
    #
    #n_digits = len(str(nb))
    #
    #for i in range(nb):
    #    turb_name = 'turbine{:0{width}d}'.format(i, width=n_digits)
    #    randX = (np.random.rand()*(xmax-xmin)) + xmin
    #    randY = (np.random.random()*(ymax-ymin)) + ymin
    #    Z = -20.0
    #    xyz_turb = np.asarray((randX, randY, Z))
    #    xyz_turb = np.asarray((randX, randY, Z))
    #    Data[turb_name] = {}
    #    Data[turb_name]['position'] = xyz_turb
    
    #staggered layout1
    #Xa = [15.0, 115.0, 215.0]
    #Xb = [65.0, 165.0]
    #Y = [50.0, 300.0, 550.0, 800.0, 1050.0, 1300.0]
    
    #staggered layout2
    #Xa = np.arange(10.0, 230.0, 70.0)
    #Xb = np.arange(45.0, 230.0, 70.0)
    #Y = np.arange(50.0, 1500.0, 175.0)
    
    #staggered layout3
    Xa = np.arange(15.0, 230.0, 50.0)
    Xb = np.arange(40.0, 230.0, 50.0)
    Y = np.arange(50.0, 1500.0, 125.0)
    
    row = 1.0
    nb = 0
    
    n_digits = len(str(len(Y)))
    
    for i in range(len(Y)):
        
        turb_name = 'turbine{:0{width}d}'.format(nb, width=n_digits)
        
        if row > 0.0:
            for x in Xa:
                xyz = np.asarray((Y[i], x, -20.0))
                Data[turb_name] = {}
                Data[turb_name]['position'] = xyz
                print x, Y[i], nb
                print xyz
                nb += 1
        else:
            for x in Xb:
                xyz = np.asarray((Y[i], x, -20.0))
                Data[turb_name] = {}   
                Data[turb_name]['position'] = xyz
                print x, Y[i], nb
                print xyz
                nb += 1
        
        row = row * (-1.0)
    
    print "Number of turb. :", nb
    
    #Save as pickle
    fpath = os.path.join(input_dir, 'turb_pos.p')
    f = file(fpath, 'wb')
    pkl.dump(Data, f, protocol=pkl.HIGHEST_PROTOCOL)
    f.close()

if __name__ == '__main__':
    
    input_dir = 'inputs_tidal'
    main(input_dir)





