import numpy as np
import h5py
import os
['data_folder']

def create_empty_project(folder, prj_name):
#    general_inputs = {'angle_def':[1],
#					  'frequency_def':[1,1,1],
#					  'mooring_dof':[0],
#					  'pto_dof':[0],
#					  'ndof':[0],
#					  'water_depth':[0],
#					  'get_array_mat':[0],
#					  'cyl_nzeta':[0],
#					  'cyl_ntheta':[0],
#                               'data_folder': '',
#                               'input_type': 2}
#	
#    body = {'ID': [0],
#			'axis_angles':[0,0,0],
#			'cog': [0,0,0],
#			'mass': [1],
#			'inertia': [[0,0,0],[0,0,0],[0,0,0]],
#			'mesh': '',
#			'parent_body': [-1],
#			'child_dof_position': [[]],
#			'dof_with_parent': [[]]}
    general_inputs = {'angle_def':[],
					  'frequency_def':[],
					  'mooring_dof':[],
					  'pto_dof':[],
					  'ndof':[],
					  'water_depth':0,
					  'get_array_mat':0,
					  'cyl_nzeta':[],
					  'cyl_ntheta':[],
                               'data_folder': '',
                               'input_type': 1}
	
    body = -1
			
    bodies = body
    body_inputs = {'shared_dof':[0,0,0,0,0,0],
					'local_cs': [],
					'body': bodies}
    filename = os.path.join(folder, ''.join([prj_name, '_data_collection.hdf5']))

    dic = {'inputs_hydrodynamic': {'general_input':general_inputs,
									'body_inputs': body_inputs},
			'prj_filename': filename,
			'prj_folder': folder,
			'prj_name': prj_name}
    save_dict_to_hdf5(dic, filename)
     
    return load_dict_from_hdf5(filename)

def save_dict_to_hdf5(dic, filename):
    """
    ....
    """
    with h5py.File(filename, 'w') as h5file:
        recursively_save_dict_contents_to_group(h5file, '/', dic)

def recursively_save_dict_contents_to_group(h5file, path, dic):
    """
    ....modes error
    """
    for key, item in dic.items():
        if isinstance(item, (np.ndarray, np.int64, np.float64, float, int, list, str, bytes, unicode)):
			try: 
				h5file[path + key] = item
			except:
				print key
				print item
				print(['*']*100)
        elif isinstance(item, dict):
            recursively_save_dict_contents_to_group(h5file, path + key + '/', item)
        else:
            raise ValueError('Cannot save %s type'%type(item))

def load_dict_from_hdf5(filename):
    """
    ....
    """
    with h5py.File(filename, 'r') as h5file:
        return recursively_load_dict_contents_from_group(h5file, '/')

def recursively_load_dict_contents_from_group(h5file, path):
    """
    ....
    """
    ans = {}
    for key, item in h5file[path].items():
        if isinstance(item, h5py._hl.dataset.Dataset):
            ans[key] = item.value
        elif isinstance(item, h5py._hl.group.Group):
            ans[key] = recursively_load_dict_contents_from_group(h5file, path + key + '/')
    return ans

if __name__ == '__main__':

#    data = {'x': 'astring',
#            'y': np.arange(10),
#            'd': {'z': np.ones((2,3)),
#                  'b': b'bytestring'}}
#    print(data)
#    filename = 'test.h5'
#    save_dict_to_hdf5(data, filename)
#    dd = load_dict_from_hdf5(filename)
#    print(dd)
    
    
#    data = {'x':[[1,0,0,0,0,0],[1,0,0,0,0,0,0]]}
#    print(data)
#    filename = 'test.h5'
#    save_dict_to_hdf5(data, filename)
#    dd = load_dict_from_hdf5(filename)
#    print(dd)
    
    filename = r'C:\Users\francesco\Desktop\test_project\wec_solution.h5'
    dd = load_dict_from_hdf5(filename)
    
    
    # should test for bad type