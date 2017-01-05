### Test case using analytical inputs ###
####  Note:  to  run  into  ipython  ####

##### Change input_dir path to your local version ###

import os
import cPickle as pkl

from dtocean_tidal.main import Hydro
from dtocean_tidal.main import Array
from dtocean_tidal.submodel.WakeInteraction import WakeInteraction
from dtocean_tidal.modules import ArrayYield, HydroImpact


debug = True
debug_plot = True

input_dir='inputs_tidal'
#Upload inputs
##Velocity field
fpath = os.path.join(input_dir, 'simple_inputs.p')
f = file(fpath, 'rb')
data = pkl.load(f)
f.close()

##Turbines positions
fpath = os.path.join(input_dir, 'turb_pos.p')
f = file(fpath, 'rb')
turbines = pkl.load(f)
f.close()

##Turbines features
fpath = os.path.join(input_dir, 'turb_fea.p')
f = file(fpath, 'rb')
features = pkl.load(f)
f.close()

hydro = Hydro(data, debug=debug, debug_plot=debug_plot)
array = Array(hydro, turbines, features, debug=debug, debug_plot=debug_plot)
"""
interaction = WakeInteraction(hydro, array, debug=debug, debug_plot=debug_plot)
interaction.solv_induction(debug=debug, debug_plot=debug_plot)
arrayYield = ArrayYield(array, debug=debug, debug_plot=debug_plot)
arrayYield.performance()
impacts=HydroImpact(array, hydro, interaction, arrayYield, debug=debug, debug_plot=debug_plot)
impacts.dissipated_over_available_flux_ratio(debug=debug)
pow_perf_array = arrayYield.array_capacity / 1e6
ratio = impacts.diss_avai_mass_flow_rate
print "Array performance: " + str(pow_perf_array) + " MWatts"
print "Ratio of dissipated-over-available mass flow rate: " + str(ratio*100) + " %"
"""
#### or simply
from dtocean_tidal import main
#main.test(input_dir)

