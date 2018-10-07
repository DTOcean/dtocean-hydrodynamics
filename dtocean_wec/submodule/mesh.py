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

import os
from math import pi

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


class SurF(object) :

    def __init__(self, coord) :
        """
        coord (3D numpy array) : (Number of panels,
                                  Number of nodes for each panel,
                                  Number of coordinates for each node == 3)
        """
        self.coord = coord

    def triangles(self) :
        """
        Identify triangular and quadrilateral panels as well as spurious panels
        """
        acl = np.allclose
        quad = list()
        tri = list()
        for n, pan in enumerate(self.coord) :
            if acl(pan[0],pan[1]) :
                tri.append(n)
                self.coord[n] = pan[(0,1,2,3),:]
            elif acl(pan[0],pan[3]) :
                tri.append(n)
                self.coord[n] = pan[(3,0,1,2),:]
            elif acl(pan[1],pan[2]) :
                tri.append(n)
                self.coord[n] = pan[(1,2,3,0),:]
            elif acl(pan[2],pan[3]) :
                tri.append(n)
                self.coord[n] = pan[(2,3,0,1),:]
            else :
                quad.append(n)
        self.quad = np.array(quad, dtype = int)
        self.tri = np.array(tri, dtype = int)

    def translation(self, disp) :
        """
        Translation is applied to coord according to
        a the displacement given by disp.

        disp (1D numpy array) : [translation in x,
                                 translation in y,
                                 translation in z]
        """
        self.coord += disp

    def rotation(self, rot) :
        """
        Rotation is applied to coord according to
        rot. The rotation is performend with respect to the (0,0,0)
        point in self.coord.

        rot (1D numpy array) : [rotation in x,
                                rotation in y,
                                rotation in z] (radians)
        """
        nodes = self.coord.reshape((-1,3))
        R = [np.array([[np.cos(angl),-np.sin(angl)],
                       [np.sin(angl),np.cos(angl)]],
                        dtype = float) for angl in rot]
        for ind0 in range(len(nodes)) :
            for ind1, ax in enumerate(((1,2),(0,2),(0,1))) :
                nodes[ind0, ax] = np.dot(R[ind1],nodes[ind0, ax])
        self.coord = nodes.reshape((-1,4,3))

    def Show_Panels(self) :
        """
        Plot the panels of the surface
        """
        fig = plt.figure()
        ax = fig.add_subplot(111, projection = '3d')
        for pan in self.coord:
            ax.plot(pan[:,0],pan[:,1],pan[:,2],'r-')
            # connect with a red line the last node with the first one
            ax.plot(pan[(0,-1),0],pan[(0,-1),1],pan[(0,-1),2],'r-')
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        set_axes_equal(ax)

    def Show_Norms(self, delta):
        """
        Plot the panels of the surface with norms

        delta (float) : scaling factor to plot the norm
        """
        ## plot
        self.NCA()
        fig = plt.figure()
        ax = fig.add_subplot(111, projection = '3d')
        for n, pan in enumerate(self.coord):
            ax.plot(pan[:,0],pan[:,1],pan[:,2],'r-')
            # connect with a red line the last node with the first one
            ax.plot(pan[(0,-1),0],pan[(0,-1),1],pan[(0,-1),2],'r-')
            # plot norms
            ax.plot([delta*self.norm[n,0]+self.center[n,0],self.center[n,0]],
                    [delta*self.norm[n,1]+self.center[n,1],self.center[n,1]],
                    [delta*self.norm[n,2]+self.center[n,2],self.center[n,2]],
                    'b-')
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        set_axes_equal(ax)

    def NCA(self) :
        """
        Compute normal vector (N), center point (C)
        and area (A) of each panel
        """
        ##
        self.triangles()
        S1 = self.coord[:,2]-self.coord[:,1]
        S2 = self.coord[:,3]-self.coord[:,1]
        self.norm = np.cross(S1,S2)
        self.norm = (self.norm.T/np.sqrt(np.sum(self.norm*self.norm, axis = 1))).T
        # compute centre point panel n
        self.center = np.zeros((len(self.coord), 3), float)
        self.center[self.quad] = self.coord[self.quad].mean(axis = 1)
        self.center[self.tri] = self.coord[self.tri,1:].mean(axis = 1)
        # compute area panel
        self.area = np.zeros(len(self.coord), dtype = float)
        for n, pan in zip(self.quad, self.coord[self.quad,:,:]) :
            self.area[n] = A4Pol(pan)
        for n, pan in zip(self.tri, self.coord[self.tri,:,:]) :
            S1, S2, S3 = pan[2]-pan[1], pan[3]-pan[1], pan[3]-pan[2]
            S1, S2, S3 = np.sqrt(np.dot(S1,S1)), np.sqrt(np.dot(S2,S2)), np.sqrt(np.dot(S3,S3))
            s = (S1+S2+S3)*.5
            self.area[n] = np.sqrt(s*(s-S1)*(s-S2)*(s-S3))

    def Chkud(self, refpoints) :
        """
        Check the outwardness of the normal vectors using reference points.
        If the distance between the reference point and the head of the
        normal vector is smaller than that of shifted 180 deg, this is
        interpreted as the normal vector is pointing inwards and the one
        pointing outwards is that shifted 180 deg.

        refpoints (2D numpy.array) : x-y-z coordinates of each reference point.
                                    Since they are used to decide upon the
                                    outwardness the user should select them
                                    so they are unequivocally inside the geometry.
        """
        self.NCA()
        for n , pan in enumerate(self.coord) :
            # look for closest refpoint to pan
            dist = np.sqrt(((refpoints-self.center[n])**2).sum(axis = 1))
            Pref = refpoints[dist == dist.min()].reshape(-1)
            Pref -= self.center[n] # Pref coordinates from self.center[n]
            # distance**2 between Pref and head of norm
            m1 = np.dot(Pref-self.norm[n],Pref-self.norm[n])
            # distance**2 between Pref and head of norm shifted 180 deg.
            m2 = np.dot(Pref--self.norm[n],Pref--self.norm[n])
            if m1 < m2 :
                print 'panel {:} points inwards'.format(n)
                # normal is pointing to Pref
                # user must know this is interpreted as normal is pointing inwards
                if n in self.tri: #if it is a triangle the first nodes must be equal to the second
                    self.coord[n] = pan[(3,3,2,1),:]
                elif n in self.quad:
                    self.coord[n] = pan[(3,2,1,0),:]

    def GDF(self, fn, ULEN = 1, GRAV = 9.81, ISX = 0, ISY = 0) :
        """
        Generate a .GDF file for the surface

        fn (string) : direction and name of the .GDF file that will be
                 generated. e.g. ".\\DesiredFolder\\DesiredName.GDF"
        ULEN, GRAV, ISX, ISY are parameters of the GDF file formatting.
        See "WAMIT Manual .GDF files" for their definition.
        """
        fn = fn.split('.GDF')[0]
        fn += '.GDF'
        Npanels = self.coord.shape[0]
        nodes = self.coord.reshape((-1,3))
        with open(fn , 'w') as fid:
            fid.write('.GDF file format\n{:} {:} 	ULEN GRAV\n'.format(ULEN, GRAV))
            fid.write('{:} {:} 	ISX ISY\n{:}\n'.format(ISX, ISY, Npanels))
            for node in nodes :
                fid.write('  {:.6f} {:.6f} {:.6f}\n'.format(*node))
        self.dirGDF = fn

    def dat(self, fn, refined = True, sym = 0) :
        """
        Generate a .dat file for the surface

        fn (string) : direction and name of the .dat file that will be
                 generated. e.g. ".\\DesiredFolder\\DesiredName.dat"
        refinement (bool) : if True (default), the format will follow that
                            of Nemoh mesh file (the one given in Nemoh.cal).
                            If False, the format will follow that prior
                            refinement for mesh generation using Mesh.exe.
        sym (int): if refined then sym is 1 if a symmetry about the (xOz)
                   plane is used. 0 otherwise
        """
        fn = fn.split('.dat')[0]
        fn += '.dat'
        Npanels = self.coord.shape[0]
        nodes = list()
        conectivity = np.zeros((Npanels, 4), dtype = int)
        for i0, nodeori in enumerate(self.coord.reshape((-1,3))) :
            nodes.append(nodeori)
            for i1, node in enumerate(nodes) :
                if all(nodeori == node) :
                    conectivity[i0/4, int((i0/4.-i0/4)*4)] = i1+1
                    if i1+1 < len(nodes) :
                        nodes.pop()
                        break
        if refined :
            with open(fn , 'w') as fid:
                fid.write('2 {:}\n'.format(sym))
                for i0, node in enumerate(nodes) :
                    fid.write('{:} '.format(i0+1))
                    fid.write('{:} {:} {:}\n'.format(*node))
                fid.write('0 0. 0. 0.\n')
                for panel in conectivity :
                    fid.write('{:} {:} {:} {:}\n'.format(*panel))
                fid.write('0 0. 0. 0.\n')
        else :
            with open(fn , 'w') as fid:
                fid.write('{:}\n{:}\n'.format(len(nodes),Npanels))
                for node in nodes :
                    fid.write('{:} {:} {:}\n'.format(*node))
                for panel in conectivity :
                    fid.write('{:} {:} {:} {:}\n'.format(*panel))
        self.dirDAT = fn
        self.Npanels = Npanels
        self.Nnodes = len(nodes)

    def refine_mesh(self, option=2, triangles=False):
        """
        """
        self.NCA()
        crd_q, cen_q = self.coord[self.quad], self.center[self.quad]
        crd_t, cen_t = self.coord[self.tri], self.center[self.tri]
        if option == 1:
            stp = len(crd_q)*2
        else:
            stp = len(crd_q)*4
        if triangles:
            crd_ref = np.zeros((stp+len(crd_t)*3, 4, 3), dtype=float)
        else:
            crd_ref = np.zeros((stp+len(crd_t), 4, 3), dtype=float)
            crd_ref[stp:] = crd_t
        # quad
        if option == 1:
            crd_ref[:stp:2] = np.transpose(np.array([crd_q[:,0], crd_q[:,0], crd_q[:,1], crd_q[:,2]], float), axes=(1,0,2))
            crd_ref[1:stp:2] = np.transpose(np.array([crd_q[:,0], crd_q[:,0], crd_q[:,2], crd_q[:,3]], float), axes=(1,0,2))
        else:
            crd_ref[:stp:4] = np.transpose(np.array([crd_q[:,0], crd_q[:,0], crd_q[:,1], cen_q], float), axes=(1,0,2))
            crd_ref[1:stp:4] = np.transpose(np.array([cen_q, cen_q, crd_q[:,1], crd_q[:,2]], float), axes=(1,0,2))
            crd_ref[2:stp:4] = np.transpose(np.array([crd_q[:,3], crd_q[:,3], cen_q, crd_q[:,2]], float), axes=(1,0,2))
            crd_ref[3:stp:4] = np.transpose(np.array([crd_q[:,0], crd_q[:,0], cen_q, crd_q[:,3]], float), axes=(1,0,2))
        # tri
        if triangles:
            crd_ref[stp::3] = np.transpose(np.array([crd_t[:,1], crd_t[:,1], crd_t[:,2], cen_t], float), axes=(1,0,2))
            crd_ref[stp+1::3] = np.transpose(np.array([cen_t, cen_t, crd_t[:,2], crd_t[:,3]], float), axes=(1,0,2))
            crd_ref[stp+2::3] = np.transpose(np.array([crd_t[:,1], crd_t[:,1], cen_t, crd_t[:,3]], float), axes=(1,0,2))
        self.coord = crd_ref

def mesher(Nth, Nz, radius, draught, geometry, thf = 2*pi, zcut = 0., axes = (0,1,2)) :
    """
    Generates a mesh for few geometries (see geometry input) and
    exports it into .GDF format file

    Nth (integer) : angular discretization
    Nz (integer) : vertical discretization
    radius (float) : radius
    draught (float) : total draught
    fn (string) : direction and name of the .GDF file that will be
                 generated. e.g. ".\\DesiredFolder\\DesiredName.GDF"
    geometry (string) : geometry of the body that is wanted to be meshed.
                    options are: 'Cylinder',
                                 'Hemisphere',
                                 'Sphere',
                                 'Cone'
                                 'Paraboloid'
                                 'PerAWat'
             (2D numpy.array) : x-y coordinates of the nodes (rows) of the section
                                which will be further extruded along -z direction.
                                The last column is the degree of the polynomial
                                for vertical extrusion (0, 1 or 2). This must
                                be the same for all nodes (rows). The user should
                                notice that for a closed geometry, the first node
                                must be repeated the last node. Moreover, the user
                                should notice that Nth must be larger than the
                                total number of nodes (rows).
    thf (float) : end angle (in radians) as of which no more mesh
                  will be generated, leading to a sharp edge. Note
                  that for a fully mesh generation thf must
                  be equal to 2*pi
    zcut (float): zcut is substracted to draught.
    axes (tuple) : order of the coordinates x = 0 (theta = 0 deg),
                                            y = 1 (theta = 90 deg)
                                            z = 2 (axis of revolution).
                   By default is set to (0,1,2), i.e. (x,y,z)
    """
    ## Connectivity
    Id = np.array(range(Nz*Nth), dtype = int).reshape((Nz,Nth))
    con = np.array([Id[:-1,1:].reshape(-1),
                    Id[:-1,:-1].reshape(-1),
                    Id[1:,:-1].reshape(-1),
                    Id[1:,1:].reshape(-1)], dtype = int).T
    ## Coordinates
    z = np.linspace(0, - (draught - zcut), Nz)
    if type(geometry) == str :
        th = np.linspace(0, thf, Nth)
        x = radius*np.cos(th)
        y = radius*np.sin(th)
    ## Non circular sectional area
    else :
        x = np.zeros(Nth, dtype = float)
        y = np.zeros(Nth, dtype = float)
        Nn = geometry.shape[0]
        x = np.zeros(Nth, dtype = float)
        y = np.zeros(Nth, dtype = float)
        Nn = geometry.shape[0]
        x[:Nn] = geometry[:,0].copy()
        y[:Nn] = geometry[:,1].copy()
        inds = np.array(range(Nth-1), dtype = int) # number of intervals
        for i in range(Nth-Nn) :
            L = (x[1:Nn+i]-x[:Nn+i-1])**2+(y[1:Nn+i]-y[:Nn+i-1])**2
            n0 = inds[:Nn+i-1][L == L.max()][0] # initial node of the largest interval
            buff = (x[n0+1:-1].copy(),y[n0+1:-1].copy())
            if abs(x[n0+1]-x[n0]) < 1E-6:
                x[n0+1] = x[n0+1]/2-x[n0]/2.+x[n0]
                y[n0+1] = y[n0+1]/2+y[n0]/2.+y[n0]
            else :
                m = (y[n0+1]-y[n0])/(x[n0+1]-x[n0])
                x[n0+1] = x[n0+1]/2-x[n0]/2.+x[n0]
                y[n0+1] = m*x[n0+1] + y[n0]-m*x[n0]
            x[n0+2:] = buff[0]
            y[n0+2:] = buff[1]
        if geometry[0,2] == 0 :
            geometry = 'Cylinder'
        elif geometry[0,2] == 1 :
            geometry = 'Cone'
        else :
            geometry = 'Paraboloid'
    ## Non circular sectional area
    Z, X = np.meshgrid(z, x, indexing = 'ij')
    Z, Y = np.meshgrid(z, y, indexing = 'ij')
    del(x,y,z) # free some space
    X = X.reshape(-1)
    Y = Y.reshape(-1)
    Z = Z.reshape(-1)
    ## Redifine end point connectivity
    if (geometry == 'Hemisphere' or geometry == 'Cone' or\
        geometry == 'PerAWat' or geometry == 'Sphere' or\
        geometry == 'Paraboloid') and zcut == 0. :
        # in such cases the last Nth points are actually a unique point,
        # X = 0 , Y = 0, Z = -draught, since radius goes to 0 at z = -draught.
        # Thus, we are only going to take one of these Nth points,
        # the rest are discarded.
        X = X[:-(Nth-1)]; Y = Y[:-(Nth-1)]; Z = Z[:-(Nth-1)]
        # the connectivity of the last Nth-1 panels will be different.
        Npp = (Nth-1)/2 # number of panels that will be generated instead of the Nth-1 panels
        conp = np.array([Id[-2,2::2],
                         Id[-2,1:-1:2],
                         Id[-2,:-2:2],
                         Nth*(Nz-1)*np.ones(Npp)], dtype = int).T
        con = con[:Npp-(Nth-1)]
        con[-Npp:] = conp
        if geometry == 'Sphere' :
            # Do the same also for the first Nth points
            X = X[Nth-1:]; Y = Y[Nth-1:]; Z = Z[Nth-1:]
            conp = np.array([(Nth-1)*np.ones(Npp),
                             Id[1,:-2:2],
                             Id[1,1:-1:2],
                             Id[1,2::2]], dtype = int).T
            con = con[(Nth-1)-Npp:]
            con[:Npp] = conp
            con -= Nth-1
    ## List of nodes for the GDF file format
    conGDF = con.reshape(-1)
    ##
    if geometry == 'Cylinder' :
        pass
    if geometry == 'Hemisphere' :
        if abs(draught) != abs(radius) :
            print 'draught must be equal to the radius, i.e. draught == radius = True'
        X *= np.sqrt(1-Z**2/radius**2)
        Y *= np.sqrt(1-Z**2/radius**2)
    if geometry == 'Sphere' :
        if abs(draught) != abs(2*radius) :
            print 'draught must be equal to the diameter, i.e. draught == 2*radius = True'
        tr = draught/2
        X *= np.sqrt(1-(Z--tr)**2/radius**2)
        Y *= np.sqrt(1-(Z--tr)**2/radius**2)
    if geometry == 'Cone' :
        X *= (Z/draught+1)
        Y *= (Z/draught+1)
    if geometry == 'Paraboloid' :
        X *= np.sqrt(Z/draught+1)
        Y *= np.sqrt(Z/draught+1)
    if geometry == 'PerAWat' :
        H = draught - radius
        sep = Z < -H
        X[sep] *= np.sqrt(1-(Z[sep]--H)**2/radius**2)
        Y[sep] *= np.sqrt(1-(Z[sep]--H)**2/radius**2)
    ## Store X, Y, Z in aux and order it
    aux = np.array([X,Y,Z] , dtype = float).T[:,axes]
    del(X,Y,Z)
    ## Generate the mesh, with panels info and so on using the
    ## .GDF file format convention
    coord = np.zeros(conGDF.shape + (3,), dtype = float)
    for ind, node in enumerate(conGDF) :
        coord[ind,] = aux[node,]
    ## instantiate the surface using coord as (Npanels, Nnodes, Ncoord)
    coord = coord.reshape((-1,4,3))
    surface = SurF(coord)
    return surface

def circumference(Nl, radius, z, upwards = True) :
    """
    Generates a mesh for a 2D circumference

    Nl (int) : number of stripes, discretization of the circumscribed square.
               Striped pattern, with same number of vertical and horizontal lines
               and equal to Nl.
    radius (float) : radius of the circumference.
    z (float) : vertical position of the circumference.
    upwards (bool) : determines whether or not the normal vectors point upwards,
                     i.e. towards the positive z-axis.
    """
    a = 1./np.sqrt(2.)*radius
    b = a
    xn = np.linspace(-a, a, Nl)
    yn = np.linspace(-b, b, Nl)[::-1]
    # central rectangle
    (X, Y) = np.meshgrid(xn, yn, indexing = 'xy')
    X = X.reshape(-1)
    Y = Y.reshape(-1)
    aux = np.array(range(Nl*Nl), dtype = int).reshape((Nl, Nl))
    con = np.array([aux[:-1, 1:].reshape(-1),
                    aux[:-1, :-1].reshape(-1),
                    aux[1:, :-1].reshape(-1),
                    aux[1:, 1:].reshape(-1)], dtype = int).T
    coord0 = np.ones(((Nl-1)*(Nl-1), 4, 3), dtype = float)
    coord0[:, :, 0] = X[con]
    coord0[:, :, 1] = Y[con]
    # borders
    def borders(X, Y, Xe, Ye, N) :
        """
        """
        aux = np.array(range(N*8), dtype = int).reshape((8, N))
        con =  np.array([aux[:-1:2, 1:].reshape(-1),
                         aux[:-1:2, :-1].reshape(-1),
                         aux[1::2, :-1].reshape(-1),
                         aux[1::2, 1:].reshape(-1)], dtype = int).T
        conxe = np.array([[0, 2, 2, 0, 4, 6, 6, 4],
                          [1, 3, 3, 1, 5, 7, 7, 5],
                          [1, 2, 2, 1, 5, 6, 6, 5],
                          [0, 2, 2, 0, 4, 6, 6, 4]], dtype = int).T
        conye = np.array([[0, 2, 4, 6, 6, 5, 3, 0],
                          [0, 3, 4, 7, 6, 4, 3, 0],
                          [1, 3, 5, 6, 7, 4, 2, 1],
                          [0, 2, 4, 6, 6, 5, 3, 0]], dtype = int).T
        Np = con.shape[0] + conxe.shape[0]
        coord = np.ones((Np, 4, 3), dtype = float)
        coord[:con.shape[0], :, 0] = X[con]
        coord[:con.shape[0], :, 1] = Y[con]
        coord[con.shape[0]:, :, 0] = Xe[conxe]
        coord[con.shape[0]:, :, 1] = Ye[conye]
        return coord
    #
    bordx = np.sqrt(radius**2-yn[1:-1]**2)
    bordy = np.sqrt(radius**2-xn[1:-1]**2)
    Xe1 = np.array([xn[0], xn[1], xn[0], -bordx[0], xn[-1], xn[-2], xn[-1], bordx[0]], dtype = float)
    Ye1 = np.array([yn[0], bordy[0], yn[0], yn[1], yn[-2], yn[-1], yn[-1], -bordy[0]], dtype = float)
    if (radius-a) / (xn[1]-xn[0]) >= 1.5 :
        ones = np.ones(Nl-2, dtype = int)
        X21 = np.array([xn[1:-1], xn[1:-1], bordx[0]*ones, xn[-1]*ones, xn[-2:0:-1], xn[-2:0:-1], -bordx[0]*ones, xn[0]*ones], dtype = float).reshape(-1)
        Y21 = np.array([bordy[0]*ones, yn[0]*ones, yn[1:-1], yn[1:-1], -bordy[0]*ones, yn[-1]*ones, yn[-2:0:-1], yn[-2:0:-1]], dtype = float).reshape(-1)
        coord1 = borders(X21, Y21, Xe1, Ye1, Nl-2)
        ones = np.ones(Nl-4, dtype = int)
        X22 = np.array([xn[2:-2], xn[2:-2], bordx[1:-1], bordx[0]*ones, xn[-3:1:-1], xn[-3:1:-1], -bordx[-2:0:-1], -bordx[0]*ones], dtype = float).reshape(-1)
        Y22 = np.array([bordy[1:-1], bordy[0]*ones, yn[2:-2], yn[2:-2], -bordy[-2:0:-1], -bordy[0]*ones, yn[-3:1:-1], yn[-3:1:-1]], dtype = float).reshape(-1)
        Xe2 = np.array([xn[1], xn[2], -bordx[0], -bordx[1], xn[-2], xn[-3], bordx[0], bordx[1]], dtype = float)
        Ye2 = np.array([bordy[0], bordy[1], yn[1], yn[2], yn[-3], yn[-2], -bordy[0], -bordy[1]], dtype = float)
        coord2 = borders(X22, Y22, Xe2, Ye2, Nl-4)
    else :
        ones = np.ones(Nl-2, dtype = int)
        X2 = np.array([xn[1:-1], xn[1:-1], bordx, xn[-1]*ones, xn[-2:0:-1], xn[-2:0:-1], -bordx[::-1], xn[0]*ones], dtype = float).reshape(-1)
        Y2 = np.array([bordy, yn[0]*ones, yn[1:-1], yn[1:-1], -bordy[::-1], yn[-1]*ones, yn[-2:0:-1], yn[-2:0:-1]], dtype = float).reshape(-1)
        coord1 = borders(X2, Y2, Xe1, Ye1, Nl-2)
        coord2 = np.zeros((0, 4, 3), dtype = float) # dummy
    # assemble
    coord = z*np.ones((coord0.shape[0] + coord1.shape[0] + coord2.shape[0], 4, 3))
    coord[:coord0.shape[0], :, :2] = coord0[:,:,:2]
    coord[coord0.shape[0]:coord0.shape[0] + coord1.shape[0], :, :2] = coord1[:,:,:2]
    coord[coord0.shape[0] + coord1.shape[0]:, :, :2] = coord2[:,:,:2]
    if not upwards :
        coord = coord[:, ::-1, :]
    surface = SurF(coord)
    return surface

def readGDF(fn, translation = np.zeros(3), rotation = np.zeros(3)) :
    """
    Reads .GDF files and returns the x,y,z coordinates of the nodes
    of each panel through an instance of SurF. If translation or rotation
    are non zero arrays, a translated and/or rotated SurF
    instance will be generated instead.

    fn (string) : direction and name of the .GDF file that will be
                 generated. e.g. ".\\DesiredFolder\\DesiredName.GDF"
    translation (1D numpy array) : [translation in x,
                                    translation in y,
                                    translation in z]
    rotation (1D numpy array) : [rotation in x,
                                 rotation in y,
                                 rotation in z]
                                 with respect to the (0,0,0)
                                 point given in the mesh file. In this
                                 regard, it is worth noticing that
                                 the translation is performed before
                                 the rotation.
    """
    ## read GDF
    with open(fn,'r') as fIlE :
        lines = fIlE.readlines()
        lines = lines[4:]
        Npoints = len(lines)
        Npanels = Npoints/4
        coord = np.zeros((Npanels,4,3),dtype=float)
        for panel in range(Npanels):
            for node in range(4):
                coord[panel,node] = np.array(lines[panel*4+node].split(), dtype=float)
                if any(translation != 0) :
                    coord[panel,node] += translation
                if any(rotation != 0) :
                    for d,ax in enumerate(((1,2),(0,2),(0,1))) :
                        R = np.array([[np.cos(rotation[d]),-np.sin(rotation[d])],
                                      [np.sin(rotation[d]),np.cos(rotation[d])]], dtype = float)
                        coord[panel,node,ax] = np.dot(R,coord[panel,node,ax])
    ## instantiate the surface using coord as (Npanels, Nnodes, Ncoord)
    surface = SurF(coord)
    return surface

def readDAT(fn, translation = np.zeros(3), rotation = np.zeros(3)) :
    """
    Reads .dat mesh files and returns the x,y,z coordinates of the nodes
    of each panel through an instance of SurF. If translation or rotation
    are non zero arrays, a translated and/or rotated SurF
    instance will be generated instead.

    fn (string) : direction and name of the .dat file that will be
                 generated. e.g. ".\\DesiredFolder\\DesiredName.dat"
    translation (1D numpy array) : [translation in x,
                                    translation in y,
                                    translation in z]
    rotation (1D numpy array) : [rotation in x,
                                 rotation in y,
                                 rotation in z]
                                 with respect to the (0,0,0)
                                 point given in the mesh file. In this
                                 regard, it is worth noticing that
                                 the translation is performed before
                                 the rotation.
    """
    ## read dat
    with open(fn,'r') as fIlE :
        nodes = list()
        conectivity = list()
        line = fIlE.readline().split('!')[0].split()
        if int(line[0]) == 2 :
            sym = int(line[1])
        else :
            Nn = int(line[0])
            sym = -1
        if sym > -1 :
            trigger = False
            for fl in fIlE.readlines() :
                line = np.array(fl.split('!')[0].split(), dtype = float)
                if not trigger :
                    nodes.append(line)
                elif trigger :
                    conectivity.append(line)
                if all(line == 0) :
                    if not trigger:
                        nodes.pop()
                    elif trigger :
                        conectivity.pop()
                    trigger = True
            nodes = np.array(nodes, dtype = float).reshape((-1,4))[:,1:]
        else :
            Np = int(fIlE.readline().split('!')[0].split()[0])
            for i , fl in enumerate(fIlE.readlines()) :
                line = np.array(fl.split('!')[0].split(), dtype = float)
                if i < Nn :
                    nodes.append(line)
                elif i >= Nn and i < Np + Nn :
                    conectivity.append(line)
            nodes = np.array(nodes, dtype = float).reshape((-1,3))
        conectivity = np.array(conectivity, dtype = int).reshape((-1,4))-1
        coord = np.zeros((len(conectivity),4,3), dtype = float)
        for panel , i0 in enumerate(conectivity):
            for node in range(4) :
                coord[panel,node] = nodes[i0[node]]
                if any(translation != 0) :
                    coord[panel,node] += translation
                if any(rotation != 0) :
                    for d,ax in enumerate(((1,2),(0,2),(0,1))) :
                        R = np.array([[np.cos(rotation[d]),-np.sin(rotation[d])],
                                      [np.sin(rotation[d]),np.cos(rotation[d])]], dtype = float)
                        coord[panel,node,ax] = np.dot(R,[panel,node,ax])
    ## instantiate the surface using coord as (Npanels, Nnodes, Ncoord)
    surface = SurF(coord)
    return surface

def mergeGDFs(fns, mfn) :
    """
    fns (tuple) : each element of the tuple is the direction and
                  name of the .GDF file to be merged with the others
                  (string).
    mfn (string) : direction and name of the .GDF file that will be
                   generated. e.g. ".\\DesiredFolder\\DesiredName.GDF"
    """

    mfn = mfn.split('.GDF')[0]
    mfn += '.GDF'

    with open(mfn, 'w') as F :
        for fn in fns :
            ## read GDF
            with open(fn,'r') as fIlE :
                for line in fIlE.readlines()[4:] :
                    F.write(line)

    return mfn

def mergesurfs(Surfs) :
    """
    Surfs (tuple) : each element of the tuple is an instance of SurF
                    to be merged with the others.
    """
    nnodes = np.array([Surf.coord.shape[0]*4 for Surf in Surfs], dtype = int)
    coord = np.zeros((nnodes.sum(),3), dtype = float)
    for ind, Surf in enumerate(Surfs) :
        ini = nnodes[:ind].sum()
        fini = ini + nnodes[ind]
        coord[ini:fini,] = Surf.coord.reshape((-1,3))
    ## instantiate the surface using coord as (Npanels, Nnodes, Ncoord)
    coord = coord.reshape((-1,4,3))
    surface = SurF(coord)
    return surface

def Bounds(point, coord, norm, center, area, betaf, alphai, alphaf, Nb, Na) :
    """
    Compute the bounds of point according to a body
    defined through coord.

    point (1D numpy array) : x-y-z coordinates of the point.
    coord (3D numpy array) : x-y-z coordinates of nodes for each panel.
    norm (2D numpy array) : x-y-z coordinates normal vector for each panel.
    center (2D numpy array) : x-y-z coordinates of center point for each panel.
    area (1D numpy array) : area for each panel.
    betaf (float) : final azimuthal angle of search path.
    alphai (float) : initial zenital angle of search path.
    alphaf (float) : final zenital angle of search path.
    Nb (int) : total number of azimuths.
    Na (int) : total number of zenits.
    """
    #
    if betaf == 2*pi :
        beta = np.linspace(0,betaf,Nb, endpoint = False)
    else : beta = np.linspace(0,betaf,Nb)
    alpha = np.linspace(alphai,alphaf,Na)
    alpha, beta = np.meshgrid(alpha, beta, indexing = 'ij')
    # initialize
    S = np.inf*np.ones(alpha.shape, dtype = float)
    B = np.inf*np.ones(alpha.shape, dtype = int)
    spherical = np.array([np.cos(alpha)*np.cos(beta),
                          np.cos(alpha)*np.sin(beta),
                          np.sin(alpha)], dtype = float)
    for ip , pan in enumerate(coord) :
        # compute constant D for the equation of the plane panj
        Dj = -np.dot(norm[ip],center[ip])
        # intersection through parameter s with plane defined by panj
        denominator = norm[ip,0]*spherical[0]+norm[ip,1]*spherical[1]+norm[ip,2]*spherical[2]
        filt = abs(denominator) > 1E-6
        sij = (-Dj-np.dot(norm[ip],point))/denominator[filt]
        cond = (sij > 0) * (sij < S[filt])
        sij = sij[cond]
        filt[filt] *= cond
        # spherical-cartesian centred at Pi
        cij = (sij*spherical[:,filt]).T + point
        # intersects with panj ?
        dmax = (np.sqrt(((pan-center[ip])**2).sum(axis = 1))).max()
        dsij = np.sqrt(((cij-center[ip])**2).sum(axis = 1))
        cij = cij[dsij <= dmax]
        sij = sij[dsij <= dmax]
        filt[filt] *= dsij <= dmax
        # Area panel j with extra vertex cij called "e"
        cond = np.zeros(sij.shape , dtype = bool)
        for ie, e in  enumerate(cij):
            pane = np.zeros(pan.shape, dtype = float)
            pane[0] = e # e
            pane[1:] = pan[:-1] # a, b and c so eabc
            Areap = A4Pol(pane)
            pane[2] = pan[-1] # change b by d instead so eadc
            Areap += A4Pol(pane)
            # condition of inclusion
            if abs(area[ip]-Areap) < 1E-6 :
                cond[ie] = True
        filt[filt] *= cond
        S[filt] = sij[cond]
        B[filt] = ip
    return S, B, spherical

def A4Pol(coord) :
    """
    Calculates the area of an irregular polygon
    consisting of 4 vertices

    coord (2D numpy.array) : each row is a vertex and
                             the columns are the coordinates.
                             All the coordinates must be
                             given according to the same
                             coordinate system. Vertices must
                             be ordered either clockwise or
                             counter-clockwise.
    """
    # sides length
    (ab, bc, cd) = np.sqrt(((coord[1:]-coord[:-1])**2).sum(axis = 1))
    ad = np.sqrt(((coord[-1]-coord[0])**2).sum())
    ac = np.sqrt(((coord[2]-coord[0])**2).sum())
    Area = 0.0
    abc = (ab,bc) # triangle consisting of vertices a, b and c
    cda = (ad,cd) # triangle consisting of vertices c, d and a
    for pivot in (abc, cda) :
        aca = .5/ac*(pivot[0]**2+ac**2-pivot[1]**2) # ac is divided into aca and acc so that aca + acc == ac
        h = np.sqrt(pivot[0]**2-aca**2)
        Area += .5*ac*h
    return Area

def set_axes_equal(ax):
    '''Make axes of 3D plot have equal scale so that spheres appear as spheres,
    cubes as cubes, etc..  This is one possible solution to Matplotlib's
    ax.set_aspect('equal') and ax.axis('equal') not working for 3D.

    Input
      ax: a matplotlib axis, e.g., as output from plt.gca().
    '''

    x_limits = ax.get_xlim3d()
    y_limits = ax.get_ylim3d()
    z_limits = ax.get_zlim3d()

    x_range = x_limits[1] - x_limits[0]; x_mean = np.mean(x_limits)
    y_range = y_limits[1] - y_limits[0]; y_mean = np.mean(y_limits)
    z_range = z_limits[1] - z_limits[0]; z_mean = np.mean(z_limits)

    # The plot bounding box is a sphere in the sense of the infinity
    # norm, hence I call half the max range the plot radius.
    plot_radius = 0.5*max([x_range, y_range, z_range])

    ax.set_xlim3d([x_mean - plot_radius, x_mean + plot_radius])
    ax.set_ylim3d([y_mean - plot_radius, y_mean + plot_radius])
    ax.set_zlim3d([z_mean - plot_radius, z_mean + plot_radius])



# Example
if __name__ == '__main__' :

    # Example: Meshing a single hinged pelamis

    ## Nth, Nz, radius and axes for all the bodies
    Nth = 10
    Nz = 10
    radius = 2.5
    ax = (2,0,1)

    ## Mesh beam
    beam = mesher(Nth, Nz, radius, draught = -25,
                  geometry = 'Cylinder', thf = -pi, zcut = 0., axes = ax)

    ## hinged end
    hinge = mesher(Nth, Nz, radius, draught = -5.7,
                   geometry = 'Cone', thf = -pi, zcut = -.2, axes = ax)

    ## head
    head = mesher(Nth, Nz, radius, draught = 2.5,
                  geometry = 'Hemisphere', thf = -pi, zcut = 0., axes = ax)

    ## translate bodies to its global position
    beam.translation(np.array((-30.,0.,0.)))
    hinge.translation(np.array((-5.,0.,0.)))
    head.translation(np.array((-30.,0.,0.)))

    ## merge bodies
    pelamisL = mergesurfs((beam, hinge, head))

    ## create the other part of the pelamis
    pelamisR = SurF(pelamisL.coord.copy())
    pelamisR.rotation(np.array((pi,-pi,0.)))

    ## merge pelamis
    pelamis = mergesurfs((pelamisL,pelamisR))

    ## compute norms, centers and areas for each panel
    pelamis.Show_Norms(1)

    # shift inward normal vectors using refpoints
    S = pelamis.Chkud(refpoints = np.array([[-25,0,-1],
                                            [-12.5,0,-1],
                                            [0,0,0],
                                            [12.5,0,-1],
                                            [25,0,-1]]))
    pelamis.Show_Norms(1)

    ## Make a directory to store the .GDF files
    depot = 'ExMesh'
    try :
        os.mkdir(depot)
    except :
        pass

    ## generate a GDF file
    pelamis.GDF(os.path.join(depot,'pelamis'))

    # Example: Meshing a weird shape
    deg = 2
    geometry = np.array([[5,0.5,deg],[3.5,5,deg],[0,1.125,deg],
                         [-3.5,5,deg],[-5,0.5,deg],[-2,0,deg],
                         [-5,-0.5,deg],[-3.5,-5,deg],[0,-1.125,deg],
                         [3.5,-5,deg],[5,-0.5,deg],[2,0,deg],
                         [5,0.5,deg]], dtype = float)

    weird = mesher(Nth = 25, Nz = 10, radius = 5, draught = 5,
                  geometry = geometry, thf = 2*pi)

    weird.Show_Norms(0.5)

    ## generate .dat file
    weird.dat(os.path.join(depot,'weird'))

    ## read .dat file
    weird2 = readDAT(weird.dirDAT)

    weird2.Show_Norms(0.5)