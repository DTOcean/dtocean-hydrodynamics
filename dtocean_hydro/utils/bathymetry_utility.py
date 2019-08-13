# -*- coding: utf-8 -*-

#    Copyright (C) 2016 Francesco Ferri
#    Copyright (C) 2017-2018 Mathew Topper
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
Created on Wed May 04 16:11:29 2016

.. moduleauthor:: Francesco Ferri <ff@civil.aau.dk>
.. moduleauthor:: Mathew Topper <mathew.topper@dataonlygreater.com>
"""

import pickle
import logging

import numpy as np
from scipy import ndimage
import matplotlib.pyplot as plt
from matplotlib import _cntr as cntr
from shapely.ops import cascaded_union
from shapely.geometry import Polygon

# Start logging
module_logger = logging.getLogger(__name__)


def get_unfeasible_regions(xyz,
                           z_bounds,
                           area_thr=100,
                           g_fil=0.1,
                           debug=False,
                           debug_plot=False):

    xyz, z_bounds, stat = check_bathymetry_format(xyz, z_bounds)
    
    if stat == -1:
        errStr = ("Error[InstalDepth]: No feasible installation area "
                  "has been found for the device for the given bathymetry.")
        raise ValueError(errStr)
    elif stat == 1:
        module_logger.debug("No NOGO areas related to the machine depth "
                            "installation constraints have been found.")
        return None, False
    
    dx = np.max(xyz[1:, 0] - xyz[:-1, 0])
    dy = np.max(xyz[1:, 1] - xyz[:-1, 1])

    X, Y, Z = get_bathymetry_grid(xyz)
    safe_Z = np.nan_to_num(Z)
        
    unfeasible_mask = np.logical_not((safe_Z >= z_bounds[0]) *
                                                     (safe_Z <= z_bounds[1]))
                         
    # restructure the data in a 2d matrix and binarise the result
    data = unfeasible_mask.astype(np.int)
    label_im, labels = clustering(data,
                                  dx * dy,
                                  area_thr=area_thr,
                                  g_fil=g_fil,
                                  debug=debug)
    
    if debug:
        plt.figure(200)
        plt.imshow(label_im, cmap=plt.cm.spectral)
        plt.show()

    result_polys = get_shapely_polygons(X, Y, label_im, labels, debug=debug)
    
    return result_polys, unfeasible_mask
    

def get_shapely_polygons(X, Y, label_im, labels, debug=False):
    multi_polygon = []
    false_unfeasible = []
    
    for ind_label in labels:
        data_masked = label_im*(label_im == ind_label)
        if data_masked.max() < 1:
            continue
        c = cntr.Cntr(X, Y, data_masked)
        # trace a contour at z == 0.5
        res = c.trace(0.)
        nseg = len(res) // 2
        segments, codes = res[:nseg], res[nseg:]
        del codes
        multi_polygon.append(Polygon(segments[0]))
        if not len(segments) == 1:
            for iel in range(1,len(segments)):
                false_unfeasible.append(Polygon(segments[iel]))

    mp = None
    
    if multi_polygon:
        mp = cascaded_union(multi_polygon)
        
    if false_unfeasible:
        mp -= cascaded_union(false_unfeasible)
    
    return mp
    
    
def clustering(im, pixel_area, area_thr=10, g_fil=0.5, debug=False):
    
    # smooth the edges of the data edges
    im = ndimage.gaussian_filter(im, g_fil)
    
    # remove small areas with zeros and small areas with ones
    open_img = ndimage.binary_opening(im)
    close_img = ndimage.binary_closing(open_img)
    
    # cluster the all the independent areas
    label_im, nb_labels = ndimage.label(close_img)
    
    # calculate the areas of each cluster
    areas = ndimage.sum(im, label_im, range(nb_labels + 1))*pixel_area  # cos each pixels has an area of dx*dy

    # remove areas based on their size and the given threshold    
    mask_size = areas < area_thr
    remove_pixel = mask_size[label_im]
    label_im[remove_pixel] = 0
    
    # cleas up the IDs and short the areas
    labels = np.unique(label_im)
    label_im = np.searchsorted(labels, label_im)
    labels = np.unique(label_im) 
    
    return label_im, labels


def get_bathymetry_grid(xyz):
    
    # Get coordinates
    xi = np.unique(xyz[:, 0])
    yj = np.unique(xyz[:, 1])
    
    # Build default grids
    X, Y = np.meshgrid(xi, yj, indexing="ij")
    Z = np.zeros([len(xi), len(yj)]) * np.nan
    
    # Fill in the depth values
    for record in xyz:
        
        xidxs = np.where(xi == record[0])[0]
        assert len(xidxs) == 1
        xidx = xidxs[0]

        yidxs = np.where(yj == record[1])[0]
        assert len(yidxs) == 1
        yidx = yidxs[0]

        Z[xidx, yidx] = record[2]
    
    return X, Y, Z


def check_bathymetry_format(xyz, bound):
    status = 0  # bathymetry constraints active
    if not isinstance(xyz, np.ndarray):
        raise IOError("The data type of the bathymetry needs to be a numpy.ndarray")
    if len(xyz.shape) == 1:
        nx, ny = 1, len(xyz)
        xyz = xyz.reshape(nx,ny)
    else:
        nx, ny = xyz.shape
        
    if not (nx==3 or ny==3):
        raise IOError("The size of the bathymetry needs to be [nx3], where n is the number of datapoints")
        
    elif not ny==3:
        xyz = xyz.T
        
    if np.all(xyz[:,2] < 0):
        xyz = xyz*np.array([1,1,-1])
    
    bound = np.abs(bound).tolist()
    if bound[0] > bound[1]:
        bound = bound[::-1]
    
    # check trivial solutions
    if min(bound) == 0 and max(bound) == np.inf:  # the all area is feasible
        status = 1
    elif min(bound) == max(bound):  # none of the area is feasible
        status = -1
    else:
        unfeasible_mask = np.logical_not((xyz[:,2]>=bound[0])*(xyz[:,2]<=bound[1]))
        if np.all(unfeasible_mask):  # again none of the area is feasible
            status = -1
        elif not np.any(unfeasible_mask):  # again the all area is feasible
            status = 1
                
    return xyz, bound, status
    

if __name__ == "__main__":
    from Visualise_polygons import *
    from shapely.geometry import Point
    
    z_bounds = [0, np.inf]

    # load a test bathymetry
    xyz = pickle.load( open('bathy_test.p', 'rb') ).T
    xyz = xyz*np.array([1,1,-1],'f')
        
    # generate a random array
    test_array = np.random.rand(1000,2)*100.
    
    # test the bathymetry function
    mp, unf = get_unfeasible_regions(xyz, z_bounds, debug=True)
    mask = np.array([mp.contains(Point(xx, yy)) for xx, yy in test_array])
    
    fig, ax = plt.subplots(1, 1)
    plotCompositePolygon(mp, ax=ax)
    
    ax.plot(test_array[mask==False, 0], test_array[mask==False, 1], "rx")
    ax.plot(test_array[mask, 0], test_array[mask, 1], "ro")
    
    
    
    

