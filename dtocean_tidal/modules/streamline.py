# -*- coding: utf-8 -*-

# Copyright (c) 2011 Raymond Speth.
# Copyright (C) 2016 Thomas Roc
# Copyright (C) 2017-2018 Mathew Topper

# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from shapely.geometry import Point



class Streamlines:
    """
    Compute a set of streamlines covering the given velocity field.
    X and Y - 1D or 2D (e.g. generated by np.meshgrid) arrays of the
    grid points. The mesh spacing is assumed to be uniform
    in each dimension.
    U and V - 2D arrays of the velocity field.
    maxlen - The maximum length of an individual streamline segment.

    Plots are generated with the 'plot' or 'plotArrows' methods.
    """
    def __init__(self, data, turbPos, NbTurb, maxlen=np.inf, debug=False):
        #Define attributs 
        self._debug = debug
        res = 0.25  # Sets the distance between successive points in each
                    # streamline
        detectLoops = True  # Determines whether an attempt is made to stop 
                            # extending a given streamline before reaching 
                            # maxLen points if it forms a closed loop or 
                            # reaches a velocity node.
        self.detectLoops = detectLoops
        self.maxlen = maxlen
        self.res = res
        #allocate data
        self.lease = data['lease']
        self.x = data['X'][:]
        self.y = data['Y'][:]
        self.dx = (self.x[-1]-self.x[0])/(self.x.size-1)
        self.dy = (self.y[-1]-self.y[0])/(self.y.size-1)
        self.dr = self.res * np.sqrt(self.dx * self.dy)
        #initial drifter positions and velocities
        self.driftPos = turbPos
        self.u = data['U']
        self.v = data['V']
        self.interpU = data['interpU']
        self.interpV = data['interpV']
        
        # marker for which regions have contours
        self.used = np.zeros(self.u.shape, dtype=bool)
        self.used[0] = True
        self.used[-1] = True
        self.used[:,0] = True
        self.used[:,-1] = True
        
        # Don't try to compute streamlines in regions where there is no
        # velocity data
        for i in range(self.y.size):
            for j in range(self.x.size):
                if self.u[i,j] == 0.0 and self.v[i,j] == 0.0:
                    self.used[i,j] = True
        
        # Make the streamlines
        self.streamlines = []
        for i in range(NbTurb):
            # Make a streamline starting at the first unrepresented grid point
            self.streamlines.append(self._makeStreamline(self.driftPos[i,0],
                                                         self.driftPos[i,1]))
    
    def plot(self, NbTurb, lw=1, ax=None, size=16):
        """
        Draw the computed streamline segments with arrows indicating flow
        direction.
        
        Optional keyword arguments:
            lw - line width
            ax - axes to use for plotting
            size - size of the arrow head
        
        """
        if ax is None:
            ax = plt.gca()
        
        speed = np.sqrt(self.u ** 2 + self.v ** 2)
        cs = ax.contourf(self.x, self.y, speed,  extend='both')
        ax.streamplot(self.x, self.y, self.u, self.v, color='0.8', density=2)
        
        for streamline in self.streamlines:
            if (len(streamline[0]) > 1) or (len(streamline[1]) > 1):
                ax.plot(streamline[0], streamline[1], 'k', lw=lw)
        
        #Plot points and numbers
        ax.scatter(self.driftPos[:,0], self.driftPos[:,1])
        n = range(NbTurb)
        
        for i, txt in enumerate(n):
            ax.annotate(txt, (self.driftPos[i,0], self.driftPos[i,1]))
        
        for streamline in self.streamlines:
            
            if (len(streamline[0]) > 1) or (len(streamline[1]) > 1):
                
                path = mpl.path.Path(np.asarray((streamline[0],
                                                 streamline[1])).T)
                patch = mpl.patches.FancyArrowPatch(path=path,
                                                    arrowstyle='->',
                                                    mutation_scale=size,
                                                    lw=lw,
                                                    color='k')
                ax.add_patch(patch)
        
        fig = plt.gcf()
        cbar = fig.colorbar(cs)
        
        ax.axis('tight')
        ax.set_aspect('equal')
        plt.show()
        
        return
    
    def _makeStreamline(self, x0, y0):
        """
        Compute a streamline extending in both directions from the given point.
        """

        sx, sy = self._makeHalfStreamline(x0, y0, 1) # forwards
        #rx, ry = self._makeHalfStreamline(x0, y0, -1) # backwards

        #rx.reverse()
        #ry.reverse()

        return [x0]+sx, [y0]+sy
        #return rx+[x0]+sx, ry+[y0]+sy

    def _makeHalfStreamline(self, x0, y0, sign):
        """
        Compute a streamline extending in one direction from the given point.
        """
        
        xmin = self.x[0]
        xmax = self.x[-1]
        ymin = self.y[0]
        ymax = self.y[-1]
        
        sx = []
        sy = []
        slen = 0.
        
        x = x0
        y = y0
        
        while xmin < x < xmax and ymin < y < ymax:
            
            if not Point(x, y).within(self.lease):
                break
            
            if self.detectLoops and self._detectLoop(sx, sy):
                break
            
            u = self.interpU(x, y)[0][0]
            v = self.interpV(x, y)[0][0]
            
            if u == 0. or v == 0.:
                break
            
            theta = np.arctan2(v,u)
            
            xold = x
            yold = y
            
            x += sign * self.dr * np.cos(theta)
            y += sign * self.dr * np.sin(theta)
            
            slen += np.sqrt((x - xold) ** 2 + (y - yold) ** 2)
            
            if slen > self.maxlen:
                break
            
            sx.append(x)
            sy.append(y)
        
        return sx, sy
    
    def _detectLoop(self, x, y):
        """ Detect closed loops and nodes in a streamline. """
        
        if not x or not y: return False
        
        x0 = x[-1]
        y0 = y[-1]
        
        xd = x0 - x[:-1]
        yd = y0 - y[:-1]
        
        D = np.hypot(xd, yd)
        result = (D < 0.9 * self.dr).any()
        
        return result
