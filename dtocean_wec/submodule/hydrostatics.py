# -*- coding: utf-8 -*-

#    Copyright (C) 2016 Francesco Ferri, Pau Mercadez Ruiz
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

# Imports
import numpy as np
from scipy.linalg import block_diag
import mesh as _msh
import os

def Hydrostatics_Nemohcal(directory):
    """
    """
    (geoms, modes) = ReadNemohcal_mshs_mds(directory)
    geoms = [_msh.readDAT(geom) for geom in geoms]
    for geom in geoms:
        geom.refine_mesh(option=2, triangles=False)
    cb = [[False,]]*1*len(geoms)
    cg = [[False,]]*1*len(geoms)
    Kb, Kg = Hydrostatics(geoms, modes, cb, cg, ro=1000.00, g=9.81, fixedA=False)[::2]
    return Kb+Kg

def ReadNemohcal_mshs_mds(directory) :
    """
    Read setting file Nemoh.cal and return mesh files and generalized modes.

    directory (str): path where Nemoh.cal or nemoh.cal
    """
    meshes = list()
    modes = list()
    with open(os.path.join(directory , 'Nemoh.cal') , 'r') as f_N :
        content = f_N.readlines()[6:]
        Nb = int(content[0].split('!')[0]) # Number of bodies
        strt = 1
        for body in range(Nb):
            strt += 1 #--- Body body -------------------------------------------------------------------------------------------------------------------
            mesh_file = content[strt].split('!')[0].split()[0]
            mesh_file = mesh_file.replace('"', '')
            mesh_file = mesh_file.replace("'", '')
            if mesh_file[0] == '.':
                meshes.append(directory+mesh_file[1:])
            elif meshes[0][0] == os.path.join('a','a')[1]:
                meshes.append(directory+mesh_file)
            else:
                meshes.append(os.path.join(directory, mesh_file))
            Ndof = int(content[strt+2].split('!')[0])
            DoF = [np.array(content[strt+3+ndof].split('!')[0].split(), float) for ndof in range(Ndof)]
            Nforce = int(content[strt+3+Ndof].split('!')[0])
            FoRcE = [np.array(content[strt+4+Ndof+nforce].split('!')[0].split(), float) for nforce in range(Nforce)]
            modes.append([DoF, FoRcE])
            Nadd = int(content[strt+4+Ndof+Nforce].split('!')[0])
            strt += 5+Ndof+Nforce+Nadd
    return (meshes, modes)

def Hydrostatics(meshes, modes, cb, cg, ro=1000.00, g=9.81, fixedA=False):
    """
    """
    ##
    Kbuoy, Fbuoy, Kgrav, Fgrav = list(), list(), list(), list()
    for nb, mesh in enumerate(meshes):
        ## generate mesh instance in case mesh is a GDF or DAT file
        if type(mesh) == str:
            if mesh[-3:] == 'GDF' or mesh[-3:] == 'gdf':
                mesh = _msh.readGDF(mesh)
            elif mesh[-3:] == 'DAT' or mesh[-3:] == 'dat':
                mesh = _msh.readDAT(mesh)
        ## center of gravity
        mesh.NCA()
        if len(cg[nb]) == 3:
            xg = cg[nb]
        else:
            xg, Volg = center_vol(mesh)
            Volg = Volg.mean()
            xg = 1./2./Volg*xg
        ## get rid of the mesh above z = 0
        for p, pan in enumerate(mesh.coord):
            if any(pan[:, -1] > 0.):
                mesh.coord[p] = np.nan
        mesh.coord = (mesh.coord[mesh.coord == mesh.coord]).reshape((-1, 4, 3))
        ## center of buoyancy (submerged volume)
        mesh.NCA()
        xb, Vol = center_vol(mesh)
        Vol = Vol.mean()
        xb = 1./2./Vol*xb
        if len(cb[nb]) == 3:
            xb = cb[nb]
        ## calculate integrals, coordinates of center of buoyancy and submerged volume
        cond = mesh.norm[:, 2] != 0
        S = integrals(mesh.norm[cond, 2], mesh.center[cond, 0], mesh.center[cond, 1], mesh.area[cond])
        ## calculate Buoyancy force
        DOF, FORCE = len(modes[nb][0]), len(modes[nb][1])
        kbuoy, kgrav = np.zeros((FORCE, DOF), dtype = float), np.zeros((FORCE, DOF), dtype = float)
        fbuoy, fgrav = np.zeros((FORCE, DOF), dtype = float), np.zeros((FORCE, DOF), dtype = float)
        for do, dof in enumerate(modes[nb][0]):
            for fo, force in enumerate(modes[nb][1]):
                kbuoy[fo, do], fbuoy[fo, do] =  ForceB(dof, force, Vol, xb, S, fixedA = fixedA)
                kgrav[fo, do], fgrav[fo, do] =  ForceG(dof, force, xg, fixedA = fixedA)
        ##
        Kbuoy.append(ro*g*kbuoy)
        Fbuoy.append(ro*g*fbuoy)
        Kgrav.append(ro*Vol*g*kgrav)
        Fgrav.append(ro*Vol*g*fgrav)
    return block_diag(*Kbuoy), Fbuoy, block_diag(*Kgrav), Fgrav

def ForceB(dof, force, Vol, xb, S, fixedA = False):
    """
    normalized as force/rho/gravity
    """
    ## Determine if the dof is a rotation or a displacement
    rot = True
    if dof[0] == 1:
        rot = False
    ##
    moment = True
    if force[0] == 1:
        moment = False
    ##
    axf = force[1:4]
    xA, yA, zA = force[4:]
    axd = dof[1:4]
    xB, yB, zB = dof[4:]
    xb, yb, zb = xb
    Sx, Sy, Sxx, Syy, Sxy, Sw = S
    ##
    K = np.zeros(3, dtype = float)
    f = np.zeros(3, dtype = float)
    if not moment:
        f[2] = -Vol
        if not rot:
            K[2] = axd[2]*Sw
        elif rot:
            K[2] = -axd[1]*(Sx-xB*Sw) + axd[0]*(Sy-yB*Sw)
    elif moment:
        if not rot:
            K[0] = axd[2]*(Sy-yA*Sw)
            K[1] = -axd[2]*(Sx-xA*Sw)
            if fixedA:
                K[0] += axd[1]*Vol
                K[1] += -axd[0]*Vol
            f[0] = -(yb-yA)*Vol
            f[1] = (xb-xA)*Vol
        elif rot:
            if fixedA:
                K[0] = axd[0]*(-(zb-zB)*Vol+(Syy-(yB+yA)*Sy+yB*yA*Sw)) - axd[1]*(Sxy-yA*Sx-xB*Sy+xB*yA*Sw) + axd[2]*(xb-xB)*Vol # axd[0]*(-(zb-zB)*Vol+np.dot((y-yB)*(y-yA), nz*Aw)) - axd[1]*np.dot((x-xB)*(y-yA), nz*Aw) + axd[2]*(xb-xB)*Vol
                K[1] = -axd[0]*(Sxy-xA*Sy-yB*Sx+yB*xA*Sw) + axd[1]*(-(zb-zB)*Vol+(Sxx-(xB+xA)*Sx+xB*xA*Sw)) + axd[2]*(yb-yB)*Vol # -axd[0]*np.dot((y-yB)*(x-xA), nz*Aw) + axd[1]*(-(zb-zB)*Vol+np.dot((x-xB)*(x-xA), nz*Aw)) + axd[2]*(yb-yB)*Vol
            else:
                K[0] = axd[0]*(-(zb-zA)*Vol+(Syy-(yB+yA)*Sy+yB*yA*Sw)) - axd[1]*(Sxy-yA*Sx-xB*Sy+xB*yA*Sw) + axd[2]*(xb-xA)*Vol # axd[0]*(-(zb-zA)*Vol+np.dot((y-yB)*(y-yA), nz*Aw)) - axd[1]*np.dot((x-xB)*(y-yA), nz*Aw) + axd[2]*(xb-xA)*Vol
                K[1] = -axd[0]*(Sxy-xA*Sy-yB*Sx+yB*xA*Sw) + axd[1]*(-(zb-zA)*Vol+(Sxx-(xB+xA)*Sx+xB*xA*Sw)) + axd[2]*(yb-yA)*Vol # -axd[0]*np.dot((y-yB)*(x-xA), nz*Aw) + axd[1]*(-(zb-zA)*Vol+np.dot((x-xB)*(x-xA), nz*Aw)) + axd[2]*(yb-yA)*Vol
            f[0] = -(yb-yA)*Vol
            f[1] = (xb-xA)*Vol
    return -np.dot(K, axf), -np.dot(f, axf)

def ForceG(dof, force, xg, fixedA = False):
    """
    norndimensionalized as force/mass/gravity
    """
    ## Determine if the dof is a rotation or a displacement
    rot = True
    if dof[0] == 1:
        rot = False
    ##
    moment = True
    if force[0] == 1:
        moment = False
    ##
    axf = force[1:4]
    xA, yA, zA = force[4:]
    axd = dof[1:4]
    xB, yB, zB = dof[4:]
    xg, yg, zg = xg
    ##
    K = np.zeros(3, dtype = float)
    f = np.zeros(3, dtype = float)
    if not moment:
        f[2] = -1.
    elif moment:
        if not rot:
            if fixedA:
                K[0] = axd[1]
                K[1] = -axd[0]
            f[0] = -(yg-yA)
            f[1] = xg-xA
        elif rot:
            if fixedA:
                K[0] = axd[2]*(xg-xB) - axd[0]*(zg-zB)
                K[1] = axd[2]*(yg-yB) - axd[1]*(zg-zB)
            else:
                K[0] = axd[2]*(xg-xA) - axd[0]*(zg-zA)
                K[1] = axd[2]*(yg-yA) - axd[1]*(zg-zA)
            f[0] = -(yg-yA)
            f[1] = xg-xA
    return np.dot(K, axf), np.dot(f, axf)

def integrals(nz, x, y, Aw):
    """
    """
    Sx = np.dot(x, nz*Aw)
    Sy = np.dot(y, nz*Aw)
    Sxx = np.dot(x*x, nz*Aw)
    Syy = np.dot(y*y, nz*Aw)
    Sxy = np.dot(x*y, nz*Aw)
    Sw = np.dot(nz, Aw)
    return Sx, Sy, Sxx, Syy, Sxy, Sw

def center_vol(mesh):
    """
    """
    ## Submerged volume
    vol = (mesh.center.T*mesh.norm.T*mesh.area).sum(axis = 1)
    ## Buoyancy coordinates
    xb = (mesh.center.T**2*mesh.norm.T*mesh.area).sum(axis = 1)
    return xb, vol