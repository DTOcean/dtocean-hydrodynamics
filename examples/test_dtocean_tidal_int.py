### Test case using analytical inputs ###
####  Note:  to  run  into  ipython  ####

##### Change paths if needed in f =file('where_ever_that_is...*.p', 'rb')
 
import os
import cPickle as pkl
 
from dtocean_tidal.main import wp2_tidal

def main(input_dir):

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
    
    ##Check dosumentation
    
    ##Spin model in debug mode
    (pow_perf_dev_no_int,
     pow_perf_dev,
     pow_perf_array_no_int,
     pow_perf_array,
     ratio,
     ti) = wp2_tidal(data, turbines, features, debug=True, debug_plot=True)

if __name__ == '__main__':
    
    input_dir = 'inputs_tidal'
    main(input_dir)
    
