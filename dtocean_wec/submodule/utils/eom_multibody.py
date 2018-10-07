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

from sympy import symbols, Matrix, cos, sin
from numpy import ndarray

def linearisation(M, q, x):
    # linear analysis
    x0 = {key: value for (key, value) in zip(x,[0]*len(x))}
    Mlin = M.subs(x0)
    #Clin = - h.jacobian(q).subs(x0)
    
    return Mlin
    
def lagrangian_equation(T, q, dq):
    # non-potential energy is considered as external force applied to the generalised dofs
    u = Matrix([0]*len(q))
    intermediate = T.jacobian(dq).T
    
    M = intermediate.jacobian(dq)
    # removed the unuesed time derivative and the simplify which take up a lot of memory
    #M = simplify(intermediate.jacobian(dq))
    #h = simplify(-T.jacobian(dq).T.jacobian(q) * dq -T.jacobian(dq).T.diff(t) + T.jacobian(q).T - V.jacobian(q).T + u )
#    h1 = intermediate.jacobian(q) * dq
#    h2 = T.jacobian(q).T
#    h3 = V.jacobian(q).T
#    h =  - h1 + h2 - h3 + u 
    
    return M
    
def set_global_poa_and_local_cs(bodies, platform_tr, lCS):
    Ir0P1 = platform_tr + lCS
    # local Position from point of connection to child point of connection
    kn_r_PnPn1 = []
    for body in bodies:
        if body['number_child'] > 0:
            for child in range(body['number_child']):
                kn_r_PnPn1.append(body['child_dof_position'][child])
    
    return Ir0P1, kn_r_PnPn1

def set_data_for_mb_analysis(bodies, lCS, shared_dof):
    bodies = sorted(bodies, key=lambda k: k['ID'])
    lCS = Matrix(lCS)        
    extended_dof = 0
    for el in bodies:
        extended_dof += el['extended_dof']
    ndof = len(shared_dof)+extended_dof
    
    q = Matrix(symbols(['q'+str(el) for el in range(ndof)]))
    dq = Matrix(symbols(['dq'+str(el) for el in range(ndof)]))
    t = symbols('t')
    x = q.col_join(dq)
    if isinstance(shared_dof, ndarray):
        shared_dof = shared_dof.tolist()
    
    return bodies, lCS, ndof, q, dq, x, t, shared_dof


def evaluate_potential_energy(bodies, Ir0Sn, Iomega0Sn, k_hst):
    VV = Matrix([0.])
    for el in Ir0Sn:
        rotations = [ell[-1] for ell in Iomega0Sn if ell[0]==el[0]][0]
        state = el[-1].col_join(rotations)
        VV += state.T * k_hst[el[0]] * state
    
    VV = Matrix([.5*VV])
    return VV

def evaluate_kinetic_energy(IvSn, KnomegaIKn, m, inertia):    
    # assess the kinetic energy of the system
    TT = Matrix([0.])
    for el in IvSn:
        K_omega = [ell[-1] for ell in KnomegaIKn if ell[0]==el[0]][0]
        TT += m[el[0]] * el[1].T*el[1] + K_omega.T * inertia[el[0]] * K_omega
        
    TT *= 0.5
    
    return TT

def set_global_velocity_of_cog(Ir0Sn, q, dq, t):
    IvSn = []
    for el in Ir0Sn:
        IvSn.append([el[0], el[1].jacobian(q) * dq + el[1].diff(t)])
    return IvSn
    
def set_global_state_of_cog(bodies, rot_matrix, Ir0Pn, Iomega0Pn):
    Ir0Sn = []
    Iomega0Sn = []
    for body in bodies:
        R = [el[-1] for el in rot_matrix if el[0]==body['ID']][0]
        P_pos = [el[-1] for el in Ir0Pn if el[0]==body['ID']][0]
        P_rot = [el[-1] for el in Iomega0Pn if el[0]==body['ID']][0]
        cog_pos = Matrix(body['cog'])
        Ir0Sn.append([body['ID'], P_pos + R * cog_pos]) 
        Iomega0Sn.append([body['ID'], P_rot])
        
    return Ir0Sn, Iomega0Sn

def set_global_state_of_application_points(bodies, rot_matrix, dof_rot4mat, Ir0P1):
    Ir0Pn = [[0, Ir0P1]]
    Iomega0Pn = [next(x for x in dof_rot4mat if x[0]==0)]
    for body in bodies:
        if body['number_child'] > 0:
            R = [el[-1] for el in rot_matrix if el[0]==body['ID']][0]
            for el in body['child_dof_position']:
                chID = el[0]
                Ir0Pn_parent = next(x for x in Ir0Pn if x[0]==body['ID'])
                rot_child = next(x for x in dof_rot4mat if x[0]==chID)
                Iomega0Pn_parent = next(x for x in Iomega0Pn if x[0]==body['ID'])
                Ir0Pn.append([chID, Ir0Pn_parent[-1] + R* Matrix(el[1:])])
                Iomega0Pn.append([chID, Iomega0Pn_parent[-1] + rot_child[-1]])
    return Ir0Pn, Iomega0Pn


def set_local_angular_velocity(bodies, shared_dof, q, dq, ind):
    KnomegaIKn = []
    dof_rot4mat = []
    for body in bodies:
        if body['extended_dof'] == 0:
            platform_rotational_dof = [el for el in shared_dof if el[0]==2]
            platform_rot = Matrix([0,0,0])
            dof_rot = Matrix([0,0,0])
            if platform_rotational_dof:
                for r_dof in platform_rotational_dof:
                    if r_dof[1:4] == [1,0,0]:
                        platform_rot[0] = dq[ind]
                        dof_rot[0] = q[ind]
                        ind += 1
                    elif r_dof[1:4] == [0,1,0]:
                        platform_rot[1] = dq[ind]
                        dof_rot[1] = q[ind]
                        ind += 1
                    elif r_dof[1:4] == [0,0,1]:
                        platform_rot[2] = dq[ind]
                        dof_rot[2] = q[ind]
                        ind += 1
                    else:
                        raise IOError("Shared degree of freedom not understood")
                KnomegaIKn.append([body['ID'], platform_rot])
                dof_rot4mat.append([body['ID'], dof_rot])
            else:
                KnomegaIKn.append([body['ID'], Matrix([0,0,0])])
                dof_rot4mat.append([body['ID'], Matrix([0,0,0])])
        else:
            rotational_dof = [el for el in body['dof_with_parent'] if el[0]==2]
            bodies_rot = Matrix([0,0,0])
            dof_rot = Matrix([0,0,0])
            if rotational_dof:
                for r_dof in rotational_dof:
                    if r_dof[1:4] == [1,0,0]:
                        bodies_rot[0] = dq[ind]
                        dof_rot[0] = q[ind]
                        ind += 1
                    elif r_dof[1:4] == [0,1,0]:
                        bodies_rot[1] = dq[ind]
                        dof_rot[1] = q[ind]
                        ind += 1
                    elif r_dof[1:4] == [0,0,1]:
                        bodies_rot[2] = dq[ind]
                        dof_rot[2] = q[ind]
                        ind += 1
                    else:
                        raise IOError("Shared degree of freedom not understood")
                KnomegaIKn.append([body['ID'], bodies_rot])
                dof_rot4mat.append([body['ID'], dof_rot])
        
    rot_matrix = []         
    for el in dof_rot4mat:
        base_angle = Matrix(next(x['axis_angles'] for x in bodies if x['ID']==el[0]))
        rot_matrix.append([el[0], get_rot_matrix_from_euler(el[-1]+base_angle)])
    
    return KnomegaIKn, dof_rot4mat, rot_matrix 
    
    
def set_platform_translational_dof(shared_dof, q):
    # assuming the platform will always have only the canonical dofs
    platform_translational_dof = [el for el in shared_dof if el[0]==1]
    platform_tr = Matrix([0,0,0])
    ind = 0
    if platform_translational_dof:
        for tr_dof in platform_translational_dof:
            if tr_dof[1:4] == [1,0,0]:
                platform_tr[0] = q[ind]
                ind += 1
            elif tr_dof[1:4] == [0,1,0]:
                platform_tr[1] = q[ind]
                ind += 1
            elif tr_dof[1:4] == [0,0,1]:
                platform_tr[2] = q[ind]
                ind += 1
            else:
                raise IOError("Shared degree of freedom not understood")
    return platform_tr, ind

def get_rot_matrix_from_euler(ang):
    Az = lambda angle: Matrix([[cos(angle), -sin(angle), 0], [sin(angle), cos(angle), 0], [0, 0, 1]])
    Ay = lambda angle: Matrix([[cos(angle), 0, sin(angle)], [0, 1, 0], [-sin(angle), 0, cos(angle)]])
    Ax = lambda angle: Matrix([[1, 0, 0], [0, cos(angle), -sin(angle)], [0, sin(angle), cos(angle)]])
    
    return Ax(ang[0])*Ay(ang[1])*Az(ang[2])