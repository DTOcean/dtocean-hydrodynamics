# -*- coding: utf-8 -*-

#    Copyright (C) 2016 Francesco Ferri
#    Copyright (C) 2022 Mathew Topper
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
import re

import matplotlib.pyplot as plt
from numpy import (array,
                   cos,
                   dot,
                   cross,
                   int_,
                   mean,
                   pi,
                   reshape,
                   sin,
                   vstack,
                   zeros)
from numpy import linalg as LA
from mpl_toolkits.mplot3d import axes3d


class MeshBem():
    """
    Mesh_BEM: class used to open, visualise, transform and save structured
        meshes
    
    Args:
        file_name (str): file name of the mesh to be read
    
    Optional args:
        path (str): location of the mesh file
    
    Attributes:
        file_name (str): file name of the mesh to be read
        path (str): location of the mesh file
        mesh_fn (str): full path name to the mesh file
        xsim (int): index identifying the symmetry around the x-axis
        Connectivity (numpy.ndarray): vertex ID of each panel, listed in a 
            counter-clockwise direction from the fluid perspective
        Vertex (numpy.ndarray): x,y,z coordinates of the mesh vertex
        panels (list): list of panel objects
        nP (int): number of panel of the mesh
        nV (int): number of vertex of the mesh
    """
    
    def __init__(self, file_name, path=""):
        
        self.file_name = file_name
        self.path = path
        self.mesh_fn = os.path.join(path, file_name)
        
        extension = self.mesh_fn[-3:].lower()
        
        if extension == "gdf":
            xsim, vertices, connectivity = read_WAMIT(self.mesh_fn)
        elif extension == "dat":
            xsim, vertices, connectivity = read_NEMOH(self.mesh_fn)
        else:
            raise IOError("Mesh file type not supported. Use GDF or dat.")
        
        self.xsim = xsim
        self.Vertex = vertices
        self.Connectivity = int_(connectivity)
        self.nP = len(connectivity)
        self.nV = len(vertices)
        self.panels = [Panel(self.Vertex[panel, :])
                                                for panel in self.Connectivity]
    
    def translate(self, x, y, z):
        """
        translate: translates the mesh to the given point
        
        Args:
            x (float) [m]: x-axis translation
            y (float) [m]: y-axis translation
            z (float) [m]: z-axis translation
        """
        self.Vertex += array([x,y,z])
        self.panels = [Panel(self.Vertex[panel, :])
                                                for panel in self.Connectivity]
    
    def rotate(self, rotZ, pivot):
        """
        rotate: rotates the mesh around the given pivoting point
        
        Args:
            rotZ (float)[rad]: rotation angle
            pivot (numpy.ndarray) [m]: pivoting point coordinates
        """
        v = self.Vertex.copy()
        v -= pivot 
        
        R = array([[cos(rotZ), -sin(rotZ), 0],
                   [sin(rotZ), cos(rotZ), 0],
                   [0, 0, 1]])
        
        self.Vertex = dot(R,v.T).T+pivot
        self.panels = [Panel(self.Vertex[panel, :])
                                                for panel in self.Connectivity]
    
    def invertNorm(self,index=-1):
        """
        invertNorm: used to invert the direction of the specified panel norm
        
        Args:
            index: index of the panel subject to the transformation. If -1 all
                   the panels will be transformed
        """
        if index==-1:
            index = range(self.nP)
            
        for pn in index:
            self.panels[pn].invert_direction()
        
        Vertex = zeros((self.nP*4,3))
        ind_v = -1
        
        for _p in self.panels:
            for _v in range(4):
                ind_v +=1
                Vertex[ind_v,:] = array([_p.x[_v], _p.y[_v], _p.z[_v]],
                                        dtype=float)
       
        Connectivity = zeros((self.nP,4))
        vertex = 1
        
        for panel in range(self.nP):
            Connectivity[panel,:] = array([vertex,vertex+1,vertex+2,vertex+3])
            vertex +=4
        
        self.Connectivity = int_(Connectivity)
        self.Vertex = array(Vertex)
    
    def visualise_mesh(self, a=None, f=None):
        """
        visualise_mesh: plot the mesh grid
        
        Optional args:
            a (matplotlib axes): parent figure axes pointer
            f (matplotlib figure): parent figure pointer
        """
        if a is None:
            f, a = plt.subplots(2, 2)
            a[1,1] = f.add_subplot(224, projection="3d")
        
        a[0,0].margins(0.05)
        for elm in self.panels:
            elm.show(a[0,0],dimension=1)
        
        a[0,0].set_aspect('equal')
        plt.grid()
        
        a[0,1].margins(0.05)
        for elm in self.panels:
            elm.show(a[0,1],dimension=2)
            
        a[0,1].set_aspect('equal')
        plt.grid()
        
        a[1,0].margins(0.05)
        for elm in self.panels:
            elm.show(a[1,0],dimension=3)
        
        a[1,0].set_aspect('equal')
        plt.grid()
        
        for elm in self.panels:
            elm.show(a[1,1])
        
        a[1,1].set_aspect('equal')
        plt.show()
    
    def visualise_norm(self, scale=1):
        """
        visualise_mesh: plot the mesh grid
        
        Optional args:
            scale (float): scaling factor for the norm. Only for visualisation
        """
        f = plt.figure()
        a = f.add_subplot(111, projection='3d')
        for elm in self.panels:
            #elm.show(a)
            elm.norm(scale)
            elm.show_norm(a)
        a.set_aspect('equal')
        plt.show()
    
    def mesh_generation(self, output_format, output_path=None):
        """
        mesh_generation: generate the mesh file specified in the output_format
        
        Args:
            output_format (str): defines the output format of the saved files.
                            gdf: wamit standard
                            nemoh: nemoh dat file
                            mesh: nenoh dat file for hydrostatic calculation
        
        Optional args:
            output_path (str): name of the location where to save the file
        """
        
        if output_path is None: output_path = self.path
        
        if output_format == "gdf":
            
            file_n = '{}gdf'.format(self.file_name[0:-3])
            file_p = os.path.join(output_path, file_n)
            
            with open(file_p, 'w') as f:
                
                f.write('WAMIT mesh file.\n')
                f.write('1 9.82    ULEN GRAV\n')
                f.write('{} 0 ISX ISY\n'.format(self.xsim))
                f.write('{}\n'.format(self.nP))
                
                for vertex in range(self.nV):
                    f.write('{} {} {}\n'.format(self.Vertex[vertex, 0],
                                                self.Vertex[vertex, 1],
                                                self.Vertex[vertex, 2]))
        
        elif output_format == "nemoh":
            
            file_n = '{}dat'.format(self.file_name[0:-3])
            file_p = os.path.join(output_path, file_n)
            
            with open(file_p, 'w') as f:
                
                f.write('2 0\n')
                
                for vertex in range(self.nV):
                    f.write('{} {} {} {}\n'.format(vertex+1,
                                                   self.Vertex[vertex, 0],
                                                   self.Vertex[vertex, 1],
                                                   self.Vertex[vertex, 2]))
                
                f.write('0 0.00 0.00 0.00\n')
                
                for panel in range(self.nP):
                    f.write('{} {} {} {}\n'.format(*(
                                            self.Connectivity[panel,:] + 1)))
                
                f.write('0 0 0 0\n')
        
        elif output_format == "mesh":
            
            file_n = '{}_mesh.dat'.format(self.file_name[0:-4])
            file_p = os.path.join(output_path, file_n)
            
            with open(file_p, 'w') as f:
            
                f.write('{}\n'.format(self.nV))
                f.write('{}\n'.format(self.nP))
                
                for vertex in range(self.nV):
                    f.write('{} {} {}\n'.format(self.Vertex[vertex, 0],
                                                self.Vertex[vertex, 1],
                                                self.Vertex[vertex, 2]))
                
                for panel in range(self.nP):
                    msg = '{} {} {} {}\n'.format(self.Connectivity[panel, 0],
                                                 self.Connectivity[panel, 1],
                                                 self.Connectivity[panel, 2],
                                                 self.Connectivity[panel, 3])
                    f.write(msg)
        
        else:
            
             raise ValueError("Unsupported output format. Use gdf, nemoh or "
                              "mesh")


class Panel():
    """
    Panel: the class is the base object of the mesh class. Each panel is defined by the four vertexs

    Args:
        vertex (numpy.ndarray): coordinates of the vertexs

    Attributes:
        x (numpy.ndarray): x coordinates of the vertexs
        y (numpy.ndarray): y coordinates of the vertexs
        z (numpy.ndarray): z coordinates of the vertexs
        centroid (numpy.ndarray): panel centroid
        area (float): panel area
        dr (numpy.ndarray): normal direction
        d (numpy.ndarray): normal direction applied at the panel centroid
    """
    def __init__(self,vertex):
        self.x = vertex[[0,1,2,3],0]
        self.y = vertex[[0,1,2,3],1]
        self.z = vertex[[0,1,2,3],2]
        self.centroid = mean(vertex,axis=0)
        self.norm()
    
    def invert_direction(self):
        """
        invert_direction: invert the direction of the panel vertex, inverting the normal direction

        """
        self.x = self.x[-1::-1]
        self.y = self.y[-1::-1]
        self.z = self.z[-1::-1]
        
        self.norm()

    def norm(self, scale=1):
        """
        norm: calculated the panel norm

        Optional args:
            scale: scaling factor of the normal direction
        """
        x = self.x
        y = self.y
        z = self.z
        
        v12 = array([x[1]-x[0],y[1]-y[0],z[1]-z[0]])
        v14 = array([x[3]-x[0],y[3]-y[0],z[3]-z[0]])
        v34 = array([x[3]-x[2],y[3]-y[2],z[3]-z[2]])
        v32 = array([x[1]-x[2],y[1]-y[2],z[1]-z[2]])
        
        
        n1 = cross(v12,v14)
        n2 = cross(v34,v32)
        d1 = n1/LA.norm(n1)
        d2 = n2/LA.norm(n2)
        
        self.area = (LA.norm(n1)+dot(d1,d2)*LA.norm(n2))/2.
        
        self.dr = (((d1+d2)/2.*self.area))*scale
        self.d = self.dr+self.centroid
        
    def show_norm(self,a):
        """
        show_norm: plot the norm or the panel into the given axes

        Args:
            a (matplotlib axes): axes of the parent figure pointer
        """
        a.plot_wireframe(self.x,self.y,self.z,color="#000000")
        a.plot_wireframe([self.centroid[0],self.d[0]],[self.centroid[1],self.d[1]],[self.centroid[2],self.d[2]],color="red")
        a.scatter(self.centroid[0],self.centroid[1],self.centroid[2],c="r", marker="o")
        #a.scatter(self.d[0],self.d[1],self.d[2],c="g", marker=">")
    
    def show(self,a,dimension=None):
        """
        show: plot the panel into the given axes

        Args:
            a (matplotlib axes): axes of the parent figure pointer

        Optional args:
            dimension (int):
        """
        if not dimension is None:
            V = [self.x,self.y,self.z]
            Vred = [el for ind,el in enumerate(V) if not ind==dimension-1]
            a.plot(Vred[0],Vred[1],'k')
        else:                
            a.plot_wireframe(self.x[[0,1,2,3,0]],self.y[[0,1,2,3,0]],self.z[[0,1,2,3,0]])
    
    def translate(self, delt):
        """
        translate: translates the panel by the amount specified in delt

        Args:
            delt (numpy.ndarray): transformation vector
        """
        self.x = self.x+delt[0]
        self.y = self.y+delt[1]
        self.z = self.z+delt[2]
        self.centroid = array([mean(self.x),
                               mean(self.y),
                               mean(self.z)])
        self.norm(scale=1)

    def rotate(self,rotZ, pivot):
        """
        rotate: rotates the panel by the amount specified in rotZ, arount the pivoting point

        Args:
            rotZ (float) [rad]: angle defining the transformation
            pivot (numpy.ndarray): pivoting point coordinates
        """
        vert = vstack((self.x, self.y, self.z)).T
        R = array([[cos(rotZ), -sin(rotZ), 0],
                   [sin(rotZ), cos(rotZ), 0],
                   [0, 0, 1]])
        vert_r = dot(R,(vert-pivot).T).T+pivot

        self.x = vert_r[:,0]
        self.y = vert_r[:,1]
        self.centroid = array([mean(self.x),
                               mean(self.y),
                               mean(self.z)])
        self.norm(scale=1)

    def update(self,delt):  # not used now
        """
        unused
        :param delt:
        :return:
        """
        self.x = self.x+delt[0]
        self.y = self.y+delt[1]
        self.z = self.z+delt[2]


def read_NEMOH(f_n):
    
    with open(f_n,'r') as mesh_f:
        msg = mesh_f.read()
    
    msg = strip_comments(msg)
    lines = msg.split('\n')
    
    first_line = array(lines.pop(0).split(), dtype=int)
    
    if len(first_line) == 1:
        
        xsim = 0
        nV = first_line[0]
        nP = int(lines.pop(0))
        vertices = zeros((nV, 3))
        
        for vertex in range(nV):
            vertices[vertex, :] = array(lines.pop(0).split(), dtype=float)
        
        connectivity = zeros((nP, 4))
        
        for panel in range(nP):
            connectivity[panel, :] = array(lines.pop(0).split(),
                                           dtype=float) - 1
    
    elif len(first_line) == 2:
        
        xsim = first_line[1]
        vertices = []
        connectivity = []
        pass_to_connectivity = False
        
        for line in lines:
            
            if not line: continue
            
            temp = array(line.split(), dtype=float)
            
            if not int(temp[0]) == 0 and not pass_to_connectivity:
                vertices.append(temp[1:].tolist())
            else:
                pass_to_connectivity = True
                connectivity.append([v - 1 for v in temp.tolist()])
        
        vertices = array(vertices)
        connectivity.pop(0)
        connectivity.pop(-1)
    
    connectivity = array(connectivity, dtype=int)
    
    return xsim, vertices, connectivity


def read_WAMIT(f_n):
    
    with open(f_n,'r') as mesh_f:
        msg = mesh_f.read()
    
    msg = strip_comments(msg)
    lines = msg.split('\n')
    lines = lines[2:]
    
    xsim = array(lines.pop(0).split()[0], dtype=int)
    nP = array(lines.pop(0), dtype=int)
    nV = nP * 4
    
    points = array([float(val) for entries in lines
                                               for val in entries.split()])
    vertices = reshape(points, (-1, 3))
    
    assert vertices.shape[0] == nV
    
    connectivity = zeros((nP, 4))
    vertex = 0
    
    for panel in range(nP):
        connectivity[panel, :] = array([vertex,
                                        vertex + 1,
                                        vertex + 2,
                                        vertex + 3])
        vertex +=4
    
    connectivity = array(connectivity, dtype=int)
    
    return xsim, vertices, connectivity


def strip_comments(code):
    msg = re.sub('\s*!.*', '', code, re.MULTILINE)
    msg = re.sub(r'^\n', '', msg)
    return re.sub(r'\n\s*\n', '\n', msg, re.MULTILINE)


if __name__== "__main__":

    path = r"C:\Users\francesco\Desktop\Nemoh_run"
    name = "Cylinder.GDF"

    m = MeshBem(name, path=path)
    m.visualise_mesh()
    m.translate(5,0,0)
    m.rotate(90./180*pi, array([0,0,0]))
    m.visualise_mesh()
    m.rotate(-90./180*pi, array([0,0,0]))
    m.rotate(90./180*pi, m.Vertex.mean(0))
    m.visualise_mesh()

