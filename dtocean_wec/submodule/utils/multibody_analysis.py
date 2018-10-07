# -*- coding: utf-8 -*-

#    Copyright (C) 2016 Francesco Ferri
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.

"""
Created on Mon May 16 14:35:04 2016

.. moduleauthor:: Francesco Ferri <ff@civil.aau.dk>
"""
# TODO: include translation link dofs

import numpy as np

class MultiBodyAnalysis():
    def __init__(self, bodies_description, shared_dof_bin, local_cs):
        self.bodies = sorted(bodies_description, key=lambda k: k['ID'])
        self.shared_dof_binary = shared_dof_bin
        if not isinstance(local_cs, list):
            if isinstance(local_cs, np.ndarray):
                local_cs = local_cs.tolist()
                
        self.local_cs = local_cs
        
        self._n_bodies = len(bodies_description)
        self.body_dof = len(np.where(shared_dof_bin)[0])
        
        self.b_mat = None
        self._input_bodies = None
        self.path_matrix = None
        self.igs = None
        self._n_joints = None
        self.joints_pos = None
        self.joints_pos_abs = None 
        self.rr = None
             
    
    def set_multibody_data(self, second_stage=False):
        bodies = self.bodies
        lCS = self.local_cs
        
        n_joints = sum([el['number_child'] for el in bodies]) + 1  # the +1 is due to the floating joint at the root
        joints_pos = [[-1,[-1,]+lCS]]
        if second_stage:
            joints_pos[0][1][0] = 0
                
        for el in bodies:
            if el['number_child'] > 0:
                for eli in el['child_dof_position']:
                    el_id = el['ID']
                    if isinstance(el_id, list): 
                        joints_pos.append([el_id[0], eli])
                    else:
                        joints_pos.append([el_id, eli])
        
        self._n_joints = n_joints
        self.joints_pos = joints_pos
        
        
        
    def get_multibody_topology(self):
        
        # evaluate the preliminary path matrix
        # this is used to reorganise the bodies (nodes) tree from the root to the leafs
        preliminary_path_matrix = self.get_path_matrix( self.joints_pos, self._n_bodies, self._n_joints )

        # renumber the bodies in order to obtain an easier path_matrix                  
        o_bodies = self.ordering_multibody_tree(self.bodies, preliminary_path_matrix)
        
        # reorganise the new body list and the old one in the object
        self._input_bodies = [el for el in self.bodies]
        self.bodies = o_bodies
        
        # reset the initial datga, the ID of the root is now 0
        self.set_multibody_data(second_stage=True)
        
        # evaluate the final path_matrix
        path_matrix = self.get_path_matrix( self.joints_pos, self._n_bodies, self._n_joints )
        
        # evaluate the secondary path matrix
        Igs = self.get_secondary_path_matrix( path_matrix, self._n_bodies )
        
        self.path_matrix = path_matrix
        self.igs = Igs
        
        
    def set_joints_position_orientation(self):
        n_bodies = self._n_bodies
        n_joints = self._n_joints
        new_bodies = self.bodies
        Igs = self.igs
        
        joints_pos_abs = []
        rr = np.zeros((n_bodies,)*2).tolist()
        for j in range(n_joints):
            base_node_j = np.where(Igs[:,j]==-1)[0][0]
            b_j = new_bodies[base_node_j]
            R = self.rot_matrix(b_j['axis_angles'])
            for i in range(n_bodies):
                if Igs[i,j]==0:
                    continue
                rr[i][j] = R
        
        j_m = n_joints
        j_c = 0
        edge_j = 0
        
        p_partial = np.array([0,0,0.])
        while j_c < j_m:
            node_i = np.where(Igs[:,edge_j]==-1)[0][0]
            if node_i == 0:
                j_c += 1
                edge_j = j_c+0
                joints_pos_abs.append(p_partial+np.asarray(self.joints_pos[0][1][1:]))
                p_partial = np.array([0,0,0.])
            else:
                parent_index = new_bodies[node_i]['parent_body']
                if isinstance(parent_index, list): parent_index = parent_index[0]
                parent = new_bodies[parent_index]
                R = self.rot_matrix(parent['axis_angles'])
                J_pos = next(x for x in parent['child_dof_position'] if x[0]==node_i)
                J_pos_parent = R.dot(J_pos[1:])
                p_partial = R.dot(p_partial) + J_pos_parent
                edge_j = np.where(Igs[node_i,:].T==-2)[0][0]
                
        self.joints_pos_abs = joints_pos_abs
        self.rr = rr

    def eval_velocity_transformation_matrix(self):
        joints_pos = self.joints_pos
        Igs = self.igs
        new_bodies = self.bodies
        loose_dof_b = self.body_dof
        n_joints = self._n_joints
        joints_pos_abs = self.joints_pos_abs
        rr = self.rr
        n_bodies = self._n_bodies
        shared_dof_bin = self.shared_dof_binary
        
        loose_dof = loose_dof_b*self._n_bodies
        
        
        joint_dof = []
        for joint in joints_pos:
            dof_with_parent = []
            if new_bodies[int(joint[1][0])]['extended_dof'] > 0:
                dof_with_parent = new_bodies[joint[1][0]]['dof_with_parent']
            joint_dof.append(dof_with_parent)
            
        Dij = np.zeros(Igs.shape).tolist()
        I = np.eye(3)
        Z = np.eye(3)*0.
        Bm = np.zeros((loose_dof, loose_dof_b))
        for edge in range(n_joints):
            actual_joint = joints_pos_abs[edge]
            actual_joint_dof = joint_dof[edge]
            n_joint_dof = len(actual_joint_dof)
            b_col = np.zeros((loose_dof, n_joint_dof))
            
            for node in range(n_bodies):
                if Igs[node, edge] == 0: continue
                
                joint_ind = np.where(Igs[node]==-1)[0][0]
                joint = joints_pos_abs[joint_ind]
                bi = new_bodies[node]
                
                Rbi2I = np.eye(3)
                for j in range(joint_ind+1):
                    if not type(rr[node][j]) == np.ndarray:
                        continue
                    Rbi2I = Rbi2I.dot(rr[node][j])
                dij = Rbi2I.dot(bi['cog']) + joint - actual_joint
                Dij[node][edge] = dij
                
                if edge == 0:  # floating platform            
                    Bij = np.c_[np.r_[I,Z], np.r_[-self.tilde(dij), I]]
                    # reduce the order accordlying to the real number of dofs of the floating
                    # platform                
                    bij = Bij[np.ix_(np.where(shared_dof_bin)[0], np.where(shared_dof_bin)[0])]
                    Bm[loose_dof_b*node:loose_dof_b*(node+1), :loose_dof_b] = bij
                else:
                    for idof, dof in enumerate(actual_joint_dof):
                        ui = dof[1:]
                        ui = Rbi2I.dot(ui)
                        if dof[0]==1:  # prismatic joint
                            Bij = np.r_[ui, Z[:,0]]
                        else:  # revolute joint
                            Bij = np.r_[-self.tilde(dij).dot(ui), ui]
                        
                        bij = Bij[np.ix_(np.where(shared_dof_bin)[0])]
                        b_col[loose_dof_b*node:loose_dof_b*(node+1), idof] = bij
                            
            Bm = np.c_[Bm, b_col]
            
            self.b_mat = Bm
            self.dij = Dij
            
    def get_dofs_floating_bodies(self):
        base_dofs = self.set_shared_dof(self.shared_dof_binary)
        dofs_floating_bodies = []
        for ibi, bi in enumerate(self.bodies):
            cog = self.dij[ibi][0]
            body_dofs = [el for el in base_dofs if el[0] == 1]
            body_dofs += [el[:4]+cog.tolist() for el in base_dofs if el[0] == 2]
            
            dofs_floating_bodies.append(body_dofs)
            
        return dofs_floating_bodies
        
    @staticmethod
    def ordering_multibody_tree(bodies, path_matrix):
        
        if len(bodies) == 1:
            return [el for el in bodies]
            
        root_ind = [ind for ind, el in enumerate(path_matrix) if np.all(el>=0)][0]
        n_bodies = len(bodies)
        new_num = range(n_bodies)
        old_num = [root_ind]
        base_ind = [root_ind,]
        
        iteration = -1
        while len(old_num) < n_bodies:     
            new_b_i = []    
            for b_i in base_ind:
                child = np.where(path_matrix[:,path_matrix[b_i,:]==1]==-1)[0].tolist()
                old_num += child
                new_b_i += child
            base_ind = [el for el in new_b_i]
            iteration += 1
            if iteration > 1000:
                print("Ops... something went worng: could not find the multibody path matrix")
                break

        numb = [old_num, new_num]       
        new_bodies = []
        for body in bodies:
            bn = body.copy()
            body_id = body['ID']
            if isinstance(body_id, list): body_id = body_id[0]
            bn['ID'] = numb[1][[ii for ii, el in enumerate(numb[0]) if el==body_id][0]]
            if bn['number_child'] > 0:
                old_pos = bn['child_dof_position']
                new_pos = []
                for ch_p in old_pos:
                   ch_p[0] =  numb[1][[ii for ii, el in enumerate(numb[0]) if el==ch_p[0]][0]]
                   new_pos.append(ch_p)
            new_bodies.append(bn)
        
        new_bodies = sorted(new_bodies, key=lambda k: k['ID'])
        
        return new_bodies
    
    @staticmethod
    def get_secondary_path_matrix( path_matrix, n_bodies ):
        Igs = path_matrix.copy() 
        j = 2
        while j <= n_bodies:
            col = path_matrix[:,j-1]
            ms = np.where(col==1)
            ns = np.where(col==-1)
            for m in ms:
                for k in range(j-1):
                    if Igs[m,k] < 0:
                        for n in ns:
                            Igs[n,k] += Igs[m,k]-1
            j+=1
        Igs[Igs==1] = 0
        
        return Igs
    
    @staticmethod
    def get_path_matrix( joints_pos, n_bodies, n_joints ):
        
        path_matrix = np.zeros((n_bodies, n_joints))
        for ij, joint in enumerate(joints_pos):
            for b in range(n_bodies):
                if b==joint[0]:
                    path_matrix[b, ij] = 1
                
                if b==joint[1][0]:
                    path_matrix[b, ij] = -1
                    
        return path_matrix
        
    @staticmethod
    def set_shared_dof(shared_dof_binary):
        u = [0,0,0]
        shared_dof = []
        for iel, el in enumerate(shared_dof_binary):
            if el==0: continue
            u = [0,0,0]
            u[iel%3] = 1
            shared_dof.append([1+(iel//3)]+u+[0,0,0])
            
        return shared_dof
    
    @staticmethod
    def tilde(vec):
        tilde_mat = np.array([[0, -vec[2], vec[1]],[vec[2], 0, -vec[0]],[-vec[1], vec[0], 0.]])
        return tilde_mat
    
    @staticmethod
    def rot_matrix(ang):
        Rz = lambda angle: np.array([[np.cos(angle), -np.sin(angle), 0], [np.sin(angle), np.cos(angle), 0], [0, 0, 1]])
        Ry = lambda angle: np.array([[np.cos(angle), 0, np.sin(angle)], [0, 1, 0], [-np.sin(angle), 0, np.cos(angle)]])
        Rx = lambda angle: np.array([[1, 0, 0], [0, np.cos(angle), -np.sin(angle)], [0, np.sin(angle), np.cos(angle)]])
        
        return Rz(ang[2]).dot(Ry(ang[1]).dot(Rx(ang[0])))
    
if __name__ == "__main__":
    
    n_bodies = 2
    loose_dof_b = 3
    loose_dof = loose_dof_b*n_bodies

    lCS = [0.,0.,0.]
    shared_dof_bin = [1, 0, 1, 0, 1, 0]
    floating_flag = 1 
    
    # wavebob
    b0 = {'ID':[0], 'extended_dof':0, 'cog':[0,0,-30.], 'axis_angles':[0,0,0.], 
          'parent_body':[-1], 'dof_with_parent':[], 'number_child':1,
          'child_dof_position': [[1,0.,0,0]]}
          
    b1 = {'ID':[1], 'extended_dof':1, 'cog':[0.,0,0.], 'axis_angles':[0,0,0.], 
          'parent_body':[0], 'dof_with_parent':[[1, 0, 0, 1.]], 'number_child':0,
          'child_dof_position': []}
         
    
    bodies = [b0, b1]
    obj = MultiBodyAnalysis(bodies, shared_dof_bin, lCS)
    obj.set_multibody_data()
    obj.get_multibody_topology()
    obj.set_joints_position_orientation()
    obj.eval_velocity_transformation_matrix()
    
    shared_dofs = obj.set_shared_dof(shared_dof_bin)
    bodies_dofs = obj.get_dofs_floating_bodies()
