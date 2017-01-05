import numpy as np
import sys
sys.path.append( r'C:\Users\francesco\Documents\standalonewec')

from main import BemSolution, Visualiser


# usage of the external wave module
input_type = 4
project_folder = r"C:\Users\francesco\Desktop\Pelamis_out" # absolute path where the output files will be stored
precompiled_case = 1  # select between BEM solution stored in the DB
data_folder = r"C:\Users\francesco\Desktop\Pelamis_input"  # absolute path where the input data is stored
l_cs = np.array([0.,0.,0.])
water_depth = 50

perform_pm_fitting = True
additional_stiffness = None
additional_damping = None

site_spec = {'spec_shape':'Jonswap'}
site_spec['spec_gamma'] = 1.0
site_spec['spec_spreading'] = -1
site_spec['te'] = np.linspace(3,10,1)
site_spec['hm0'] = np.linspace(0.5, 3.5, 1)
site_spec['wave_angles'] = np.linspace(0,0,1)
site_spec['probability_of_occurence'] = np.ones((1,1,1))

machine_spec = {'c_pto': 1e5}
machine_spec['k_mooring'] = 1e3
machine_spec['power_matrix'] = np.zeros((1,1,1))+1e3
machine_spec['yaw'] = 0.

wec_obj = BemSolution(input_type,
                      project_folder,
                      precompiled_case=precompiled_case,
                      data_folder=data_folder,
                      debug=True)
wec_obj.load_data(l_cs=l_cs, water_depth=water_depth)

wec_obj.pm_fitting(machine_spec, site_spec,
                   ext_k=additional_stiffness,
                   ext_d=additional_damping)
wec_obj.generate_outputs()

# TODO: check the shape of the cpto and c_ext in input
# the visualiser will be hooked at the BEM object, which unifies the 4 possible cases
plotter = Visualiser(wec_obj)
plotter.show_rao(0,0,0,0,0)