# -*- coding: utf-8 -*-
"""
Created on Tue Mar 17 12:38:08 2015

@author: 108630
"""

import make_analytical_inputs
import make_random_turbine_features
import make_random_turbine_position

def main(input_dir):
    
    make_analytical_inputs.main(input_dir)
    make_random_turbine_features.main(input_dir)
    make_random_turbine_position.main(input_dir)
    
if __name__ == '__main__':
    
    input_dir = 'inputs_tidal'
    main(input_dir)
    
