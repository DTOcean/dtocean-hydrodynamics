# -*- coding: utf-8 -*-
"""
Created on Mon Mar 21 15:56:17 2016
"""
import os
import pickle

from dtocean_hydro import start_logging
from dtocean_hydro.main import WP2

# Start the logging system
start_logging()

mod_path = os.path.realpath(__file__)
mod_dir = os.path.dirname(mod_path)

# Pick up the pickled inputs
input_dict_file = os.path.join(mod_dir, "hydrodynamics_inputs.pkl")

with open(input_dict_file, "rb") as fstream:
    iWP2input = pickle.load(fstream)

main = WP2(iWP2input,
           pickup=True,
           debug=True)

result = main.optimisationLoop()

raw_input("Press Enter to continue...")
