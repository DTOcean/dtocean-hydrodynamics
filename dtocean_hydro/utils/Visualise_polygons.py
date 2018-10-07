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
This module contains the functions used to plot the nogo zones polygons and bathymetry points

.. module:: Visualise_polygons
   :platform: Windows
   :synopsis: Polygons and points visualisation

.. moduleauthor:: Francesco Ferri <ff@civil.aau.dk>
"""
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from descartes import PolygonPatch
from shapely.geometry import Polygon

def plotPoints(ax,points,marker,alpha=1):
    """
    plotPoints: plot the Points (shapely) in the given axes

    Args:
        ax (matplotlib axes): pointer to the parent figure
        points (shapely Points): list of shapely Point to be visualise
        marker (str): marker type
        alpha (float): the alpha value defines the opacity of the points
    """
    ax.plot([p.x for p in points], [p.y for p in points],marker, alpha=alpha)


def plotConvexhull(hull, ax=None, col='#6699cc'):
    """
    plotConvexhull: plot the shape of the convex hull

    Args:
        ax (matplotlib axes): pointer to the parent figure
        hull (shapely Polygon): list of shapely Polygons to be visualise
        col (str): line color in hexadecimal format
    """
    close_flg = False
    if ax is None:
        fig, ax = plt.subplots(1, 1)
        close_flg = True
        
    if not np.shape(hull) == ():
        for el in hull:
            ax.add_patch(PolygonPatch(el, fc=col, ec='k', alpha=0.2))
    else:
        ax.add_patch(PolygonPatch(hull, fc=col, ec='k', alpha=0.2))
    
    if close_flg: plt.show()
        


def plotCompositePolygon(poly, ax=None):
    """
    plotCompositePolygon: plot the polygon describing the nogo areas

    Args:
        ax (matplotlib axes): pointer to the parent figure
        poly (shapely CompositePolygon): list of shapely CompositePolygons to be visualise
    """
    close_flg = False
    if ax is None:
        fig, ax = plt.subplots(1, 1)
        close_flg = True
    cmap = matplotlib.cm.get_cmap('Accent')
    n = 1
    if not isinstance(poly, Polygon):
        n = len(poly)
    else:
        poly=[poly,]
        
    for ic, tp in enumerate(poly):
        x,y = tp.exterior.xy
        ax.plot(x, y, color=cmap(ic/n)[:-1], alpha=0.7,
                    linewidth=3, solid_capstyle='round', zorder=2)
        ax.add_patch(PolygonPatch(tp, fc=cmap(ic/n)[:-1], ec='k', alpha=0.2))
    if close_flg: plt.show()