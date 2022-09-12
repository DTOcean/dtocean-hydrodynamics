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

import numpy as np
import matplotlib.pyplot as plt
from numpy import linalg as LA
from mpl_toolkits.mplot3d import axes3d # pylint: disable=unused-import


class MeshBem(object):
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
        xsim (int): index identifying the symmetry around the x-axis
        connectivity (numpy.ndarray): vertex ID of each panel, listed in a 
            counter-clockwise direction from the fluid perspective
        vertices (numpy.ndarray): x,y,z coordinates of the mesh vertex
        nP (int): number of panel of the mesh
        nV (int): number of vertex of the mesh
    """
    
    def __init__(self, file_name, path=""):
        
        self.file_name = file_name
        self.path = path
        
        mesh_fn = os.path.join(path, file_name)
        _, extension = os.path.splitext(mesh_fn)
        
        if extension.lower() == ".gdf":
            xsim, vertices, connectivity = read_WAMIT(mesh_fn)
        elif extension.lower() == ".dat":
            xsim, vertices, connectivity = read_NEMOH(mesh_fn)
        else:
            raise IOError("Mesh file type not supported. Use GDF or dat.")
        
        self.xsim = xsim
        self.vertices = vertices
        self.connectivity = np.int_(connectivity)
        self.nP = len(connectivity)
        self.nV = len(vertices)
    
    def translate(self, x, y, z):
        """
        translate: translates the mesh to the given point
        
        Args:
            x (float) [m]: x-axis translation
            y (float) [m]: y-axis translation
            z (float) [m]: z-axis translation
        """
        self.vertices += np.array([x,y,z])
    
    def rotate(self, rotZ, pivot):
        """
        rotate: rotates the mesh around the given pivoting point
        
        Args:
            rotZ (float)[rad]: rotation angle
            pivot (numpy.ndarray) [m]: pivoting point coordinates
        """
        
        v = self.vertices.copy()
        v -= pivot 
        
        R = np.array([[np.cos(rotZ), -np.sin(rotZ), 0],
                      [np.sin(rotZ), np.cos(rotZ), 0],
                      [0, 0, 1]])
        
        self.vertices = np.dot(R,v.T).T + pivot
    
    def invert_norm(self, index=None):
        """
        invert_norm: used to invert the direction of the specified panel norm
        
        Args:
            index: index of the panel subject to the transformation. If None
                   the panels will be transformed
        """
        
        if index is None:
            self.connectivity = np.fliplr(self.connectivity)
            return
        
        self.connectivity[index, :] = self.connectivity[index, ::-1]
    
    def visualise_mesh(self, ax=None, fig=None):
        """
        visualise_mesh: plot the mesh grid
        
        Optional args:
            ax (matplotlib axes): parent figure axes pointer
            fig (matplotlib figure): parent figure pointer
        """
        
        panels = [Panel(self.vertices[panel, :])
                                                for panel in self.connectivity]
        
        if ax is None:
            fig, ax = plt.subplots(2, 2)
            ax[1, 1] = fig.add_subplot(224, projection="3d")
        
        ax[0, 0].margins(0.05)
        
        for elm in panels:
            elm.show(ax[0, 0], dimension=1)
        
        ax[0, 0].set_aspect('equal')
        plt.grid()
        
        ax[0, 1].margins(0.05)
        
        for elm in panels:
            elm.show(ax[0, 1], dimension=2)
            
        ax[0, 1].set_aspect('equal')
        plt.grid()
        
        ax[1, 0].margins(0.05)
        
        for elm in panels:
            elm.show(ax[1, 0], dimension=3)
        
        ax[1, 0].set_aspect('equal')
        plt.grid()
        
        for elm in panels:
            elm.show(ax[1, 1])
        
        ax[1, 1].set_aspect('equal')
        plt.show()
    
    def visualise_norm(self, scale=1):
        """
        visualise_mesh: plot the mesh grid
        
        Optional args:
            scale (float): scaling factor for the norm. Only for visualisation
        """
        
        panels = [Panel(self.vertices[panel, :])
                                                for panel in self.connectivity]
        
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        
        for elm in panels:
            elm.show_norm(ax)
        
        ax.set_aspect('equal')
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
            output_path (str): name of the folder where to save the file
        """
        
        if output_path is None: output_path = self.path
        
        if output_format == "gdf":
            
            file_n = '{}gdf'.format(self.file_name[0:-3])
            file_p = os.path.join(output_path, file_n)
            
            with open(file_p, 'w') as f:
                
                f.write('WAMIT mesh file.\n')
                f.write('1 9.80665 ULEN GRAV\n')
                f.write('{} 0 ISX ISY\n'.format(self.xsim))
                f.write('{}\n'.format(self.nP))
                
                for panel in self.connectivity:
                    
                    for vertex in panel[:-1]:
                        f.write('{:15.6f} {:15.6f} {:15.6f}'.format(
                                                self.vertices[vertex, 0],
                                                self.vertices[vertex, 1],
                                                self.vertices[vertex, 2]))
                        f.write('  ')
                    
                    f.write('{:15.6f} {:15.6f} {:15.6f}'.format(
                                                self.vertices[panel[-1], 0],
                                                self.vertices[panel[-1], 1],
                                                self.vertices[panel[-1], 2]))
                    f.write('\n')
        
        elif output_format == "nemoh":
            
            file_n = '{}dat'.format(self.file_name[0:-3])
            file_p = os.path.join(output_path, file_n)
            float_format = '{:>10d} {:15.6f} {:15.6f} {:15.6f}\n'
            int_format = '{:>10d} {:>10d} {:>10d} {:>10d}\n'
            
            with open(file_p, 'w') as f:
                
                f.write('2 0\n')
                
                for vertex in range(self.nV):
                    f.write(float_format.format(vertex + 1,
                                                self.vertices[vertex, 0],
                                                self.vertices[vertex, 1],
                                                self.vertices[vertex, 2]))
                
                f.write(float_format.format(0, 0, 0, 0))
                
                for panel in range(self.nP):
                    f.write(int_format.format(
                                            *(self.connectivity[panel,:] + 1)))
                
                f.write(int_format.format(0, 0, 0, 0))
        
        elif output_format == "mesh":
            
            file_n = '{}dat'.format(self.file_name[0:-3])
            file_p = os.path.join(output_path, file_n)
            float_format = '{:15.6f} {:15.6f} {:15.6f}\n'
            int_format = '{:>10d} {:>10d} {:>10d} {:>10d}\n'
            
            with open(file_p, 'w') as f:
            
                f.write('{}\n'.format(self.nV))
                f.write('{}\n'.format(self.nP))
                
                for vertex in range(self.nV):
                    f.write(float_format.format(self.vertices[vertex, 0],
                                                self.vertices[vertex, 1],
                                                self.vertices[vertex, 2]))
                
                for panel in range(self.nP):
                    f.write(int_format.format(
                                            *(self.connectivity[panel,:] + 1)))
        
        else:
            
             raise ValueError("Unsupported output format. Use gdf, nemoh or "
                              "mesh")


class Panel(object):
    """
    Panel: a class for visualizing a 4-node BEM panel
    
    Args:
        vertex (numpy.ndarray): coordinates of the vertexs

    Attributes:
        x (numpy.ndarray): x coordinates of the vertexs
        y (numpy.ndarray): y coordinates of the vertexs
        z (numpy.ndarray): z coordinates of the vertexs
        centroid (numpy.ndarray): panel centroid
        n (numpy.ndarray): normal vector at the panel centroid
    """
    def __init__(self, vertex):
        
        self.x = vertex[[0, 1, 2, 3], 0]
        self.y = vertex[[0 ,1 ,2, 3], 1]
        self.z = vertex[[0, 1, 2, 3], 2]
        self.centroid = np.mean(vertex, axis=0)
        self.n = _get_panel_norm(self.x, self.y, self.z)
    
    def show_norm(self, ax):
        """
        show_norm: plot the norm or the panel into the given axes
        
        Args:
            ax (matplotlib axes): axes of the parent figure pointer
        """
        
        x = self.x
        y = self.y
        z = self.z
        
        v24 = np.array([x[3] - x[1],
                        y[3] - y[1],
                        z[3] - z[1]])
        v31 = np.array([x[0] - x[2],
                        y[0] - y[2],
                        z[0] - z[2]])
        
        scale = (LA.norm(v24) + LA.norm(v31)) / 8
        d = self.n * scale + self.centroid
        
        ax.plot_wireframe(x, y, z, color="#000000")
        ax.plot_wireframe([self.centroid[0], d[0]],
                          [self.centroid[1], d[1]],
                          [self.centroid[2], d[2]],
                          color="red")
        ax.scatter(self.centroid[0],
                   self.centroid[1],
                   self.centroid[2],
                   c="r",
                   marker="o")
    
    def show(self, ax, dimension=None):
        """
        show: plot the panel into the given axes
        
        Args:
            ax (matplotlib axes): axes of the parent figure pointer
        
        Optional args:
            dimension (int):
        """
        
        if not dimension is None:
            
            V = [self.x,self.y,self.z]
            Vred = [el for ind, el in enumerate(V) if ind != dimension - 1]
            ax.plot(Vred[0], Vred[1], 'k')
        
        else:
            
            ax.plot_wireframe(self.x[[0, 1, 2, 3, 0]],
                              self.y[[0, 1, 2, 3, 0]],
                              self.z[[0, 1, 2, 3, 0]])


def _get_panel_norm(x, y, z):
        
    v24 = np.array([x[3] - x[1],
                    y[3] - y[1],
                    z[3] - z[1]])
    v31 = np.array([x[0] - x[2],
                    y[0] - y[2],
                    z[0] - z[2]])
    
    dr = np.cross(v24, v31)
    
    return dr / LA.norm(dr)


def read_NEMOH(f_n):
    
    with open(f_n,'r') as mesh_f:
        msg = mesh_f.read()
    
    msg = strip_comments(msg)
    lines = msg.split('\n')
    
    first_line = np.array(lines.pop(0).split(), dtype=int)
    
    if len(first_line) == 1:
        
        xsim = 0
        nV = first_line[0]
        nP = int(lines.pop(0))
        vertices = np.zeros((nV, 3))
        
        for vertex in range(nV):
            vertices[vertex, :] = np.array(lines.pop(0).split(), dtype=float)
        
        connectivity = np.zeros((nP, 4))
        
        for panel in range(nP):
            connectivity[panel, :] = np.array(lines.pop(0).split(),
                                              dtype=float) - 1
    
    elif len(first_line) == 2:
        
        xsim = first_line[1]
        vertices = []
        connectivity = []
        pass_to_connectivity = False
        
        for line in lines:
            
            if not line: continue
            
            temp = np.array(line.split(), dtype=float)
            
            if not int(temp[0]) == 0 and not pass_to_connectivity:
                vertices.append(temp[1:].tolist())
            else:
                pass_to_connectivity = True
                connectivity.append([v - 1 for v in temp.tolist()])
        
        vertices = np.array(vertices)
        connectivity.pop(0)
        connectivity.pop(-1)
    
    connectivity = np.array(connectivity, dtype=int)
    
    return xsim, vertices, connectivity


def read_WAMIT(f_n):
    
    with open(f_n,'r') as mesh_f:
        msg = mesh_f.read()
    
    msg = strip_comments(msg)
    lines = msg.split('\n')
    lines = lines[2:]
    
    xsim = np.array(lines.pop(0).split()[0], dtype=int)
    nP = np.array(lines.pop(0), dtype=int)
    nV = nP * 4
    
    points = np.array([float(val) for entries in lines
                                               for val in entries.split()])
    vertices = np.reshape(points, (-1, 3))
    
    assert vertices.shape[0] == nV
    
    connectivity = np.zeros((nP, 4))
    vertex = 0
    
    for panel in range(nP):
        connectivity[panel, :] = np.array([vertex,
                                           vertex + 1,
                                           vertex + 2,
                                           vertex + 3])
        vertex +=4
    
    connectivity = np.array(connectivity, dtype=int)
    
    return xsim, vertices, connectivity


def strip_comments(code):
    msg = re.sub('\s*!.*', '', code, re.MULTILINE)
    msg = re.sub(r'^\n', '', msg)
    return re.sub(r'\n\s*\n', '\n', msg, re.MULTILINE)
