#!/usr/bin/python2.7
# encoding: utf-8
"""
This module contains the class used to define the array layouts.

.. module:: array
   :platform: Windows
   :synopsis: Array module for DTOcean WP2

.. moduleauthor:: Francesco Ferri <ff@civil.aau.dk>
"""

# Start logging
import logging
module_logger = logging.getLogger(__name__)

from numpy import array, zeros, inf, average, cross, meshgrid
from math import tan, cos, pi
from dtocean_wave.utils.WatWaves import len2
import matplotlib.pyplot as plt
from matplotlib.path import Path
import matplotlib.patches as patches
import numpy as np
from shapely.geometry import MultiPoint, MultiPolygon, Polygon
from dtocean_hydro.utils.Visualise_polygons import *


class Array_pkg(object):
    """
    Array_pkg: The class is used to generate nodes for the array grid and identify the feasible poisitons

    Args:
        Lease (numpy.ndarray)[m]: lease area (polygon) defined in northing-easting coordinates with mooring contraction strategy upplied
        Lease_unbuffered (numpy.ndarray)[m]: lease area (polygon) defined in northing-easting coordinates
        Dmin (tuple)[m]: minimun distance between devices in the x and y directions
        mainAngle (float)[rad]: orientation angle of the array
        NoGo_bathymetry (shapely Polygon/Multipolygon , optional)[m]: (multi)polygons object describing unfeasible areas for the device positioning
        debug (boolean, optional): if set to True, plots and additional comand line outputs are issued.

    Attributes:
        same as Args, plus:
        _lease_P (shapely.polygon): shapely polygon of the lease area verted specified in the Lease attribute
        minDist_constraint (bool): flag identifying whether the min disance constriant is violated or not
        centroidLease (numpy.ndarray)[m]: x and y coordiantes of the Lease area centroid
        nogo (list)[m]: list of polygon vertex
        _nogo_P (list): list of shapely polygons
        coord (numpy.ndarray)[m]: x,y coordiantes of the array nodes

    """
    def __init__(self, Lease, Lease_unbuffered, Dmin, mainAngle, NoGo_bathymetry=None, debug=False):
        dirwise = sum(cross(Lease[1:-1]-Lease[0], Lease[2:]-Lease[0]))
        if dirwise > 0:
            Lease[1:] = Lease[range(-1, -len(Lease), -1)]
        self.Lease = Lease
        self._lease_unbuffered = Lease_unbuffered
        self._lease_P = Polygon(Lease)
        self.Dmin = Dmin
        self.minDist_constraint = False
        self.mainAngle = mainAngle
        self.centroidLease = average(self.Lease, axis=0)
        self.Nogo_bathymetry = NoGo_bathymetry
        self.nogo = None
        self._nogo_P = None
        self.coord = None
        self._actual_mindist = None
        self._debug = debug

    def generator(self, NR, NC, IR, IC, beta, psi):
        """
        generator: Generates coordinates of the grid nodes, in function of the input arguments

        Args:
            NR (int): number of grid rows
            NC (int): number of grid columns
            IR (float)[m]: rows interdistance
            IC (float)[m]: columns interdistance
            beta (float)[rad]: angle between the rows and the main direction
            psi (float)[rad]: angle between the columns and the main direction

        Attributes:
            same as the parent class

        Note:
            The function update the following class attributes:
                coord (numpy.ndarray): x,y coordinates of the array nodes
                minDist_constraint (bool): flag identifying whether the min distance constraint is violated or not

        """
        if not np.all(self.check_grid_distance(IR, IC, beta, psi)):
            self.minDist_constraint = True
        else:
            self.minDist_constraint = False
            
        if not self.minDist_constraint:
        
            Rz = np.array([[np.cos(self.mainAngle),-np.sin(self.mainAngle)],
                                [np.sin(self.mainAngle),np.cos(self.mainAngle)]])  # 2D rotation matrix to apply main angle rotation
            i, j = np.meshgrid(np.arange(NR), np.arange(NC))
            i -= NR/2
            j -= NC/2
            
            x = IR*np.cos(beta)*i+IC*np.cos(psi)*j
            y = IR*np.sin(beta)*i+IC*np.sin(psi)*j
            
            coord_raw = np.zeros((2,NR*NC))
            coord_raw[0,:] = x.ravel()
            coord_raw[1,:] = y.ravel()
            coord = np.dot(Rz,coord_raw).T
    
            self.coord = coord+self.centroidLease  # devices translation up to Lease's centroid
        
        else:
            self.coord = np.zeros((1,2))

    
    def check_grid_distance(self, x1, x2, a1, a2):
        """
        check_grid_distance: check if the actual layout brake the constraint imposed by Dmin, using an ellipsoid
                    as base shape.
        Args:
            x1 (float)[m]: inter-column distance
            x2 (float)[m]: inter-row distance
            a1 (float)[rad]: column angle wrt the main direction
            a2 (float)[rad]: row angle wrt the main direction

        Returns:
            test (bool): if true the actual array layout is not valid

        """
        # alpha = a1-a2
        # d1 = np.sqrt(x1**2+x2**2+2*x1*x2*np.cos(alpha))
        Dmin = self.Dmin
        p1 = [x1*np.cos(a1), x1*np.sin(a1)]
        p2 = [x2*np.cos(a2), x2*np.sin(a2)]
        p3 = p2[0]-p1[0], p2[1]-p1[1]
        p4 =  p1[0]+p2[0], p1[1]+p2[1]
        
        test = ((p1[0]/(Dmin[0]/2.0))**2+(p1[1]/(Dmin[1]/2.0))**2 > 4,
                (p2[0]/(Dmin[0]/2.0))**2+(p2[1]/(Dmin[1]/2.0))**2 > 4,
                (p3[0]/(Dmin[0]/2.0))**2+(p3[1]/(Dmin[1]/2.0))**2 > 4,
                (p4[0]/(Dmin[0]/2.0))**2+(p4[1]/(Dmin[1]/2.0))**2 > 4)
                
        self._actual_mindist = np.min([x1, x2, np.sqrt(p4[0]**2+p4[1]**2)])
        return test

    def checkMinDist(self):
        """
        checkMinDist: warn the class whether the given grid nodes fulfil the minDist constraints

        Args:
            self: parent class

        Attributes:
            same as the parent class

        Note:
            the function update the following class attributes:
                minDist_constraint (bool): flag identifying whether the min distance constraint is violated or not
        """
        if self.coord is None:
            raise IOError("No coordinates provided")

        def distances(xy1, xy2):
            n = len(xy1)
            d0 = np.subtract.outer(xy1[:,0], xy2[:,0])
            d1 = np.subtract.outer(xy1[:,1], xy2[:,1])
            distance = np.hypot(d0, d1)
            return np.min(distance[distance > 0].reshape((n, n-1)), axis=1)

        dist = distances(self.coord, self.coord)
        self._actual_mindist = np.min(dist)
        if np.min(dist) < min(self.Dmin):
            self.minDist_constraint = True
        else:
            self.minDist_constraint = False

    def checkout(self, nogo_list=None):
        """
        checkout: return a boolean mask representing the feasible nodes in the self.coord attribute.
        The method is implemented base on the shapely.Polygon and shapely.MultiPoint classes.

        Args:
            nogo_list (list)[m]: list of nogo areas vertex given by the user or by other WPs

        Attributes:

        """
        if self.minDist_constraint:
            if self._debug:
                module_logger.warning('Warning[minDist]: violation of the '
                                      'minimum distance constraint in the '
                                      'actual array layout')
            machine_mask = np.zeros(self.coord.shape[0],dtype=bool)
        else:
            original_el = MultiPoint(self.coord)
            lease_mask = np.array([self._lease_P.intersects(el)
                                                for el in original_el],'bool')
            
            if np.any(lease_mask):  # identify the points inside the lease
                if not self.Nogo_bathymetry is None:  # bathymetry related nogo zones (WP2 generated)
                    # in this case the NoGo zones can have feasible zones inside,
                    # this is checked using MultiPolygons features of shapely.
                    reduced_el_array = self.coord[lease_mask]
                    reduced_mask = np.zeros(reduced_el_array.shape[0],'bool')
                    reduced_el = MultiPoint(reduced_el_array)
                    reduced_mask = np.array([
                            self.Nogo_bathymetry.intersects(el)
                                                for el in reduced_el],'bool')
                    nogo_bath_mask = lease_mask.copy()
                    nogo_bath_mask[nogo_bath_mask] = np.logical_not(
                                                                reduced_mask)
                else:
                    nogo_bath_mask = True
                                    
                if not nogo_list is None:  # other nogo zones (external)
                    self.nogo = nogo_list
                    nogo_polygons = []
                    for el in nogo_list:  # instantiates a list of Polygons for each nogo zone
                        nogo_polygons.append(Polygon(el))
                    nogo = MultiPolygon(nogo_polygons)
                    self._nogo_P = nogo_polygons
                    reduced_el_array = self.coord[lease_mask]
                    reduced_mask = np.zeros(reduced_el_array.shape[0],'bool')
                    reduced_el = MultiPoint(reduced_el_array)
                    for ng in nogo:
                        reduced_mask = np.array([ng.intersects(el)
                                for el in reduced_el],'bool') + reduced_mask
                    nogo_mask = lease_mask.copy()
                    nogo_mask[nogo_mask] = np.logical_not(reduced_mask)
                else:
                    nogo_mask = True
                    
                machine_mask = lease_mask*nogo_mask*nogo_bath_mask
            else:
                machine_mask = np.zeros(self.coord.shape[0],dtype=bool)

        return machine_mask

    def inner_region(self, points):
        """
        inner_region: used to verify whether the grid nodes adiacent to the first node passed
        are inside the ellipse generated using the Dmin attribute

        Args:
            points (numpy.ndarray): array of nodes. The first element is used as reference point

        Returns:
            out (list): list of boolean used to identify which point does not fulfil the Dmin attributes
        """
        a = self.Dmin[0]
        b = self.Dmin[1]
        nodedistance = np.zeros(4)
        for inp, point in enumerate(points[1:]):
            nodedistance[inp] = np.sqrt(((point[0]-points[0][0])**2+(point[1]-points[0][1])**2))
        nodedistance[-1] = np.sqrt(((points[1][0]-points[-1][0])**2+(points[1][1]-points[-1][1])**2))
            
        self._actual_mindist = nodedistance.min()
        return (nodedistance > np.sqrt((a*b)))

    def show(self, inside=None, ax=None):
        """
        show: visualise the feasible and unfeasible grid nodes along with the Lease area and Nogo areas

        Args:
            inside (numpy.ndarrray): boolean mask to differentiate feasible and unfeasible points
            ax (matplotlib.pyplot.axes): plot on the specified ax if given

        Returns:
            ax (matplotlib.pyplot.axes): anchor to the plot in which the data has been stacked.
        """
        if inside is None:
            inside = np.ones(len2(self.coord),dtype=bool)
        
        
        Leasep_unbuffered = zeros((len2(self._lease_unbuffered)+1,2),dtype=float)
        Leasep_unbuffered[:-1]= self._lease_unbuffered
        Leasep_unbuffered[-1]= self._lease_unbuffered[0]
        
        Leasep = zeros((len2(self.Lease)+1,2),dtype=float)
        Leasep[:-1]= self.Lease
        Leasep[-1]= self.Lease[0]
   
        if ax is None:        
            ax = plt.figure().add_subplot(111)

        if not self.Nogo_bathymetry is None:
            plotCompositePolygon(self.Nogo_bathymetry, ax)
        if not self._nogo_P is None:
            for NG in self._nogo_P:
                plotCompositePolygon(NG, ax)

        patterns = ['\\','/','-', '+', 'x', 'o', 'O', '.', '*']  # more patterns
        scale = np.mean(np.max(self.Lease,0)-np.min(self.Lease,0))
        ax.add_patch(
                        patches.Arrow(
                            self.centroidLease[0]-scale*np.cos(self.mainAngle-2*np.pi),            # x
                            self.centroidLease[1]-scale*np.sin(self.mainAngle-2*np.pi),            # y
                            0.5*scale*np.cos(self.mainAngle),            # dx
                            0.5*scale*np.sin(self.mainAngle),            # dy
                            width=0.3*scale,       # optional - defaults to 1.0
                            alpha=0.5,
                            #facecolor="#00ffff",#'none' no background
                            #linewidth=3,
                            #linestyle='dashdot',
                            #hatch=patterns[-2],
                            fill='#0000dd',
                            edgecolor='k'
                        )
                    )

        Nbody = np.shape(self.coord[inside,0])[0]
    
        for mac, Str in enumerate(range(Nbody)):
            ax.annotate(Str,(self.coord[inside,0][mac],
                             self.coord[inside,1][mac]))
        
        ax.plot(Leasep[:,0],Leasep[:,1], color='#6699cc', alpha=0.7,
                        linewidth=2, solid_capstyle='round', zorder=10000)
        
        ax.plot(self.coord[:,0],self.coord[:,1],'k+',
                    self.coord[inside,0],self.coord[inside,1],'bo',
                    Leasep_unbuffered[:,0],Leasep_unbuffered[:,1],'r',linewidth=2)
        
        

        ax.axis('equal')
        plt.show()
        return ax
        
if __name__ == "__main__":
    lease = np.array([[-100,-100],[100,-100],[100,100],[-100,100]],'f')
    arr = Array_pkg(lease, lease*1.05, (10, 10), 45.0/180*np.pi, Polygon(lease*0.2), True)
    arr.generator(10,10,20,20,0.0/180*np.pi,-45.0/180*np.pi)
    arr.show(arr.checkout())
    
    