# -*- coding: utf-8 -*-
"""
Created on Mon May 30 16:40:43 2016

@author: francesco
"""
import os
import numpy as np
import file_utilities as f_u

class DataStructure():
    def __init__(self, data_dic, force_read_flag=False):
        self.data_dic = data_dic
        self.general_inputs = data_dic['inputs_hydrodynamic']['general_input']
        self.input_type = self.general_inputs['input_type']
        self.project_folder = data_dic['prj_folder']
        self.precompiled_case = None
        self.get_array_mat = self.general_inputs['get_array_mat']
        
        if self.input_type == 2:
            self.data_folder = os.path.join(self.project_folder, 'hydrodynamic')
        else:
            self.data_folder = data_dic['inputs_hydrodynamic']['general_input']['data_folder']        
        
        if self.input_type == 4:
            self.general_inputs = data_dic['inputs_hydrodynamic']['general_input']
        elif self.input_type == 3:
            nbody = 0
            freq_def = []
            ang_def = []
            if os.path.isfile(os.path.join(self.data_folder,'nemoh.cal')):
                with open(os.path.join(self.data_folder,'nemoh.cal'),'r') as fid:
                    lines = fid.readlines()
                nbody = f_u.split_string(lines[6], int, sep=None)[0]
                water_depth = f_u.split_string(lines[3], sep=None)[0]
                line_ind = 7
                for bb in range(nbody):
                    ndof = int(f_u.split_string(lines[line_ind+3], sep=None)[0])
                    line_ind += 3+ndof*2+3
                freq_def = f_u.split_string(lines[line_ind+1], sep=None)
                ang_def = f_u.split_string(lines[line_ind+2], sep=None)
                line_ind += 2+6
                cyl_fea = [0,0,0]
                if len(lines) >= line_ind:
                    cyl_fea = f_u.split_string(lines[line_ind], sep=None)
                    
                    
                
            data_dic['inputs_hydrodynamic']['general_input']['frequency_def'] = np.array(freq_def)
            data_dic['inputs_hydrodynamic']['general_input']['angle_def'] = np.array(ang_def)
            data_dic['inputs_hydrodynamic']['general_input']['water_depth'] = np.array(water_depth)
            data_dic['inputs_hydrodynamic']['general_input']['cyl_nzeta'] = int(cyl_fea[1])
            data_dic['inputs_hydrodynamic']['general_input']['cyl_ntheta'] = int(cyl_fea[2])
            
            self.data_dic = data_dic
            self.general_inputs = data_dic['inputs_hydrodynamic']['general_input']
            self.body_inputs = data_dic['inputs_hydrodynamic']['body_inputs']
            for bi_k in self.body_inputs['body'].keys():
                ch_p = self.body_inputs['body'][bi_k]['child_dof_position']
                if isinstance(ch_p, int):
                    ch_p = self.body_inputs['body'][bi_k]['number_child'] = 0
                else:
                    ch_p = self.body_inputs['body'][bi_k]['number_child'] = len(ch_p)
            
            for bi_k in self.body_inputs['body'].keys():
                dof_p = self.body_inputs['body'][bi_k]['dof_with_parent']
                if isinstance(dof_p, int):
                    dof_p = self.body_inputs['body'][bi_k]['extended_dof'] = 0
                else:
                    dof_p = self.body_inputs['body'][bi_k]['extended_dof'] = len(dof_p)
            
            for bi_k in self.body_inputs['body'].keys():
                rel_path = self.body_inputs['body'][bi_k]['mesh']
                self.body_inputs['body'][bi_k]['mesh'] = os.path.join(self.project_folder, "raw_data", rel_path)
            
        else:
            self.body_inputs = data_dic['inputs_hydrodynamic']['body_inputs']
                      
            ang_v = np.linspace(0, 360.0, self.general_inputs['angle_def'][0], endpoint=False)
            self.general_inputs['angle_def'] = [int(self.general_inputs['angle_def'][0]), ang_v.min(), ang_v.max()] 
            
            for bi_k in self.body_inputs['body'].keys():
                ch_p = self.body_inputs['body'][bi_k]['child_dof_position']
                if isinstance(ch_p, int):
                    ch_p = self.body_inputs['body'][bi_k]['number_child'] = 0
                else:
                    ch_p = self.body_inputs['body'][bi_k]['number_child'] = len(ch_p)
            
            for bi_k in self.body_inputs['body'].keys():
                dof_p = self.body_inputs['body'][bi_k]['dof_with_parent']
                if isinstance(dof_p, int):
                    dof_p = self.body_inputs['body'][bi_k]['extended_dof'] = 0
                else:
                    dof_p = self.body_inputs['body'][bi_k]['extended_dof'] = len(dof_p)
                
            for bi_k in self.body_inputs['body'].keys():
                rel_path = self.body_inputs['body'][bi_k]['mesh']
                self.body_inputs['body'][bi_k]['mesh'] = os.path.join(self.project_folder, "raw_data", rel_path)
                
            if force_read_flag:
                self.input_type=3
        
    def check_inputs(self):
        # call the specific check for the different case
        errStr = []
        if self.input_type == 1:
            pass
        elif self.input_type == 2:
            
            if self.general_inputs['get_array_mat'] and self.general_inputs['water_depth'] == 0:
                errStr.append('the array interaction cannot be assessed with infinite water depth')
                return (False, errStr)

        elif self.input_type == 3:
            pass
        elif self.input_type == 4:
            pass
        else:
            return (False,)
            
        return (True,)            
        
    def set_inputs(self):
        folder = r"C:\Users\francesco\Desktop\test_gui\test_prj"
        gnr = os.path.join(folder, 'WEC_gnr.wp2')
        mds = os.path.join(folder, 'WEC_mds.wp2')
        mmx = os.path.join(folder, 'WEC_mmx.wp2')
        
        body_d = self.data_dic['inputs_hydrodynamic']['body_inputs']['body']
        bodies = [el for el in self.data_dic['inputs_hydrodynamic']['body_inputs']['body']]
        n_body = len(bodies)
        
        fid = open(gnr, 'w')
        fid.write("--- PTO and Mooring definition ---------------------------------------------------------------\n")
        fid.write("{}		! Number of degree of freedom\n".format(self.data_dic['inputs_hydrodynamic']['general_input']['ndof'][0]))
        for pto_i in range(len(self.data_dic['inputs_hydrodynamic']['general_input']['pto_dof'])):
            fid.write("{} ".format(self.data_dic['inputs_hydrodynamic']['general_input']['pto_dof'][pto_i]))
        fid.write("		! #DOF used to extract energy\n")
        for moo_i in range(len(self.data_dic['inputs_hydrodynamic']['general_input']['mooring_dof'])):
            fid.write("{} ".format(self.data_dic['inputs_hydrodynamic']['general_input']['mooring_dof'][moo_i]))
        fid.write("		! #DOF conneted to the mooring\n")
        fid.write("--- frequency and direction to be analysed ---------------------------------------------------\n")
        fid.write("{} {} {} ! number of frequency, min frequency, max frequency [rad/s]\n".format(*self.data_dic['inputs_hydrodynamic']['general_input']['frequency_def']))
        fid.write("{} ! Number of wave directions, for simple body 10 is ok for more complex the number should be set to 20 or bigger\n".format(self.data_dic['inputs_hydrodynamic']['general_input']['angle_def'][0]))        
        fid.close()
        
        fid = open(mds, 'w')
        fid.write("--- Shared degree of freedom -----------------------------------------------------------------------------------------------\n")
        fid.write("{} {} {} {} {} {}\n".format(*self.data_dic['inputs_hydrodynamic']['body_inputs']['shared_dof']))
        fid.write("--- Description of floating bodies -----------------------------------------------------------------------------------------------\n")
        fid.write("{}    ! Number of bodies\n".format(n_body))
        for iel, el_k in enumerate(bodies):
            el = body_d[el_k]
            dofs = el['dof_with_parent']
            child = el['child_dof_position']
            ndof_bi = 0
            nch_bi = 0
            if not isinstance(dofs, int): ndof_bi = len(dofs)
            if not isinstance(child, int): nch_bi = len(child)
            
            fid.write("------- Body {} - modes description --------\n".format(iel))
            fid.write("{}    ! body ID\n".format(el['ID'][0]))
            fid.write("{}    ! mesh file\n".format(el['mesh']))
            fid.write("{} {} {}    ! body coordinate system orientation\n".format(*el['axis_angles']))
            fid.write("{} {} {}    ! body CoG wer to the point of application\n".format(*el['cog']))
            fid.write("{}    ! parent body ID (-1 if the body is the root)\n".format(el['parent_body'][0]))
            fid.write("{}			 	! number of dof between the body and its parent\n".format(ndof_bi))
            for dof_bi in range(ndof_bi):            
                fid.write("{} {} {} {}\n".format(*dofs[dof_bi]))
            fid.write("{} 				! number of child of the body\n".format(nch_bi))
            for ch_bi in range(nch_bi):
                fid.write("{} {} {} {}\n".format(*child[ch_bi]))
        fid.write("\n")
        fid.close()
        
        fid = open(mmx, 'w')
        for iel, el_k in enumerate(bodies):
            el = body_d[el_k]
            fid.write("------- Body {} - modes description --------\n".format(iel))
            fid.write("{}    ! body ID\n".format(el['ID'][0]))
            fid.write('body mass\n')
            fid.write('{}\n'.format(el['mass'][0]))
            fid.write('inertia tensor\n')
            for rr in range(3):
                fid.write('{} {} {}\n'.format(*el['inertia'][rr,:]))
        fid.close()
        
    def get_outputs(self, res_obj):
        self.outputs = self.data_dic
        
    def __str__(self):
        return str(self.data_dic)
    

if __name__ == "__main__":
    import dtocean_wave.utils.hdf5_interface as h5i
    data = h5i.load_dict_from_hdf5(r"C:\Users\francesco\Desktop\test_gui\test_prj\test_prj_data_collection.hdf5")
    
    dataobj = DataStructure(data)
    dataobj.set_inputs()
    