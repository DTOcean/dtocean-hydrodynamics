#!/usr/bin/python2.7
# encoding: utf-8
import os
import numpy as np
import cPickle as pkl
#Global variables
#from global_variables import *
#Global variables

def main(input_dir):

    #Turbine array features
    nb = 54
    
    #For now same features for all turbines
    cut_in=1.0
    cut_out=4.5
    
    #Ct curve
    x = np.linspace(0, 10.0, 100)
    y = np.zeros(x.shape)
    for i, speed in enumerate(x):
        if speed<cut_in:
            y[i] = 0.0
        else:
            y[i] = 0.2457*speed - 0.2457
            #y[i] = a4*(speed**4) + a3*(speed**3) + a2*(speed**2) + a1*speed + a0
        if speed>cut_out:
            #y[i] = a4*(cut_out**4) + a3*(cut_out**3) + a2*(cut_out**2) + a1*cut_out + a0
            y[i] = 0.86
    
    Ct = [x, y]
    
    #Cp curve
    X = np.linspace(0, 10.0, 100)
    Y = np.ones(x.shape)*0.4#constant Cp
    Cp = [X, Y]
    
    #Features dictionaries
    turbParams = {}
    turbParams['HAS'] = 0.0 #heading angle span (deg.)
    turbParams['Cp'] = Cp #0.86 #Ct = Power coefficient
    turbParams['Ct'] = Ct #0.86 #Ct = thrust coefficient
    turbParams['Diam'] = 30.0 #Diam = rotor diameter (m)
    turbParams['cutIO'] = np.array([cut_in, cut_out])
    turbParams['floating'] = False
    
    #Same features for all turbines here but could be different
    Data = {}
    n_digits = len(str(nb))
    for i in range(nb):
        turb_name = 'turbine{:0{width}d}'.format(i, width=n_digits)
        Data[turb_name] = turbParams
    
    #Save as pickle
    fpath = os.path.join(input_dir, 'turb_fea.p')
    f = file(fpath, 'wb')
    pkl.dump(Data, f, protocol=pkl.HIGHEST_PROTOCOL)
    f.close()

if __name__ == '__main__':
    
    input_dir = 'inputs_tidal'
    main(input_dir)