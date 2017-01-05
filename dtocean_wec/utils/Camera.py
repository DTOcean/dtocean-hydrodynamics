# -*- coding: utf-8 -*-

"""
Created on Mon May 30 15:11:56 2016
@author: at

"""


from MyGeom import Point3D, Vector3D, Matrix4x4
from OpenGL.GL import *
import math
import numpy as np

class Camera:

    def __init__(self):

        self.FIELD_OF_VIEW_IN_DEGREES = 10.0
        self.ORBITING_SPEED_IN_DEGREES_PER_RADIUS_OF_VIEWPORT = 300.0

        # These are in world-space units.
        self.nearPlane = 1.0
        self.farPlane = 10000.0

        # Perspective parameters
        self.aspectratio = 1.0

        # Those are the direction fo the three axis of the camera in
        # world coordinates, used to compute the rotations necessary
        self.a =  Vector3D(1.0, 0.0, 0.0)
        self.b =  Vector3D(0.0, 1.0, 0.0)
        self.c =  Vector3D(0.0, 0.0, 1.0)
        
        # During dollying (i.e. when the camera is translating into
        # the scene), if the camera gets too close to the target
        # point, we push the target point away.
        # The threshold distance at which such "pushing" of the
        # target point begins is this fraction of nearPlane.
        # To prevent the target point from ever being clipped,
        # this fraction should be chosen to be greater than 1.0.

        self.PUSH_THRESHOLD = 1.3

        # We give these some initial values just as a safeguard
        # against division by zero when computing their ratio.

        self.viewportWidthInPixels = 10
        self.viewportHeightInPixels = 10
        self.viewportRadiusInPixels = 5

        self.sceneRadius = 10

        # point of view, or center of camera; the ego-center; the eye-point
        self.position = Point3D()

        # point of interest; what the camera is looking at; the exo-center
        self.target = Point3D()

        # This is the up vector for the (local) camera space
        self.up = Vector3D()
        self.mouse_zoom(0)

        # This is the up vector for the (global) world space;
        # it is perpendicular to the horizontal (x,z)-plane
        self.ground = Vector3D(0,1,0)
        self.reset()
        

    def reset(self):

        tangent = math.tan( self.FIELD_OF_VIEW_IN_DEGREES/2.0 / 180.0 * math.pi )
        distanceFromTarget = self.sceneRadius / tangent
        self.position = Point3D(0,0,distanceFromTarget)
#        print(self.position)
        self.target = Point3D(0,0,0)
#        print(self.target)
        self.up = self.ground.returnCopy()
        t2p = self.position - self.target
        M1 = Matrix4x4.rotationAroundOrigin( np.pi/3, Vector3D(1,0,0) )
#        print('M1',M1)
#        M2 = Matrix4x4.rotationAroundOrigin( np.pi/2, Vector3D(1,0,0) )
#        print(M2)
        M3 = Matrix4x4.rotationAroundOrigin( -np.pi/6, Vector3D(0,1,0) )
#        print('M3',M3)
#        t2p = M2 * t2p
#        t2p = M1 * t2p

#        M_temp =  Matrix4x4()
#        M_temp.m = [ 0.0, 0.0, -1.0, 0.0,
#                   -1.0, 0.0, 0.0, 0.0,
#                   0.0, 1.0, 0.0, 0.0,
#                   0.0, 0.0, 0.0, 1.0 ]

        t2p = M1 * t2p
        t2p = M3 * t2p
        self.position = self.target + t2p
        self.ground = Vector3D(0,0,1)
#        print(self.position)

    def setViewportDimensions( self, widthInPixels, heightInPixels ):

        self.viewportWidthInPixels = widthInPixels
        self.viewportHeightInPixels = heightInPixels
        self.viewportRadiusInPixels = 0.5*widthInPixels if (widthInPixels < heightInPixels) else 0.5*heightInPixels



    def getViewportWidth(self):
        return self.viewportWidthInPixels

    def getViewportHeight(self):
        return self.viewportHeightInPixels

    def setSceneRadius(self,radius):
        self.sceneRadius = radius



    def transform(self):
        tangent = math.tan( self.FIELD_OF_VIEW_IN_DEGREES/2.0 / 180.0 * math.pi )
        viewportRadius = self.nearPlane * tangent
        if self.viewportWidthInPixels < self.viewportHeightInPixels:
            viewportWidth = 2.0*viewportRadius
            viewportHeight = viewportWidth * self.viewportHeightInPixels / float(self.viewportWidthInPixels)

        else:
            viewportHeight = 2.0*viewportRadius
            viewportWidth = viewportHeight * self.viewportWidthInPixels / float(self.viewportHeightInPixels)
        glFrustum(
            - 0.5 * viewportWidth,  0.5 * viewportWidth,    # left, right
            - 0.5 * viewportHeight, 0.5 * viewportHeight,   # bottom, top
            self.nearPlane, self.farPlane
            )
        M = Matrix4x4.lookAt(self.position, self.target, self.up, False)
#        print(self.position, self.target, self.up)
#        print('M', M)
        glMultMatrixf(M.get())



    # Causes the camera to "orbit" around the target point.
    # This is also called "tumbling" in some software packages.
    def orbit(self,old_x_pixels,old_y_pixels,new_x_pixels,new_y_pixels):
        pixelsPerDegree = self.viewportRadiusInPixels / float(self.ORBITING_SPEED_IN_DEGREES_PER_RADIUS_OF_VIEWPORT)
        radiansPerPixel = 1.0 / pixelsPerDegree * math.pi / 180.0
        t2p = self.position - self.target
        M = Matrix4x4.rotationAroundOrigin( (old_x_pixels - new_x_pixels) * radiansPerPixel, self.ground )
        t2p = M * t2p
        self.up = M * self.up
        right = (self.up ^ t2p).normalized()
        M = Matrix4x4.rotationAroundOrigin( (old_y_pixels - new_y_pixels) * radiansPerPixel, right )
        t2p = M * t2p
        self.up = M * self.up
        self.position = self.target + t2p



    # This causes the scene to appear to translate right and up
    # (i.e., what really happens is the camera is translated left and down).
    # This is also called "panning" in some software packages.
    # Passing in negative delta values causes the opposite motion.
    def translateSceneRightAndUp( self, delta_x_pixels, delta_y_pixels ):
        direction = self.target - self.position
        distanceFromTarget = direction.length()
        direction = direction.normalized()

        translationSpeedInUnitsPerRadius = distanceFromTarget * math.tan( self.FIELD_OF_VIEW_IN_DEGREES/2.0 / 180.0 * math.pi )
        pixelsPerUnit = self.viewportRadiusInPixels / translationSpeedInUnitsPerRadius

        right = direction ^ self.up

        translation = right*(- delta_x_pixels / pixelsPerUnit) + self.up*(- delta_y_pixels / pixelsPerUnit)

        self.position = self.position + translation
        self.target = self.target + translation



    # This causes the camera to translate forward into the scene.
    # This is also called "dollying" or "tracking" in some software packages.
    # Passing in a negative delta causes the opposite motion.
    # If ``pushTarget'' is True, the point of interest translates forward (or backward)
    # *with* the camera, i.e. it's "pushed" along with the camera; otherwise it remains stationary.
    def dollyCameraForward( self, delta_pixels, pushTarget ):
        direction = self.target - self.position
        distanceFromTarget = direction.length()
        direction = direction.normalized()

        translationSpeedInUnitsPerRadius = distanceFromTarget * math.tan( self.FIELD_OF_VIEW_IN_DEGREES/2.0 / 180.0 * math.pi )
        pixelsPerUnit = self.viewportRadiusInPixels / translationSpeedInUnitsPerRadius

        dollyDistance = delta_pixels / pixelsPerUnit
        if not pushTarget:
            distanceFromTarget -= dollyDistance
            if distanceFromTarget < self.PUSH_THRESHOLD * self.nearPlane:
                distanceFromTarget = self.PUSH_THRESHOLD * self.nearPlane

        self.position += direction * dollyDistance
        self.target = self.position + direction * distanceFromTarget
  
    def mouse_zoom(self, inc):
        '''Convenience function to implement a zoom function.
        This is achieved by moving ``Camera.position`` in the
        direction of the ``Camera.c`` vector.
        '''
#        print(self.target)
#        print(self.position)
#        print(self.target.x())

        # Square Distance from pivot
        dsq_prime = np.sqrt(self.position.distance(self.target))
        minsq = 1.0**2  # How near can we be to the pivot
        maxsq = self.sceneRadius*7.0**2/10 # How far can we go 

        scalefac = 0.25

        if dsq_prime > maxsq and inc < 0: 
            # We're going too far
            pass

        elif dsq_prime < minsq and inc > 0:
            # We're going too close
            pass
        else:
            # We're golden
            self.position += self.position.asVector3D()*inc*scalefac


#    def autozoom(self, points):

#        '''Fit the current view to the correct zoom level to display

#        all *points*.

#  

#        The camera viewing direction and rotation pivot match the

#        geometric center of the points and the distance from that

#        point is calculated in order for all points to be in the field

#        of view. This is currently used to provide optimal

#        visualization for molecules and systems

#  

#        **Parameters**

#  

#        points: np.ndarray((N, 3))

#             Array of points.

#  

#        '''

#        extraoff = 0.01

#  

#        # Project points on the plane defined by camera up and right

#        # vector. This is achieved by using dot product on camera a

#        # and b vectors

#        abc = np.array([self.a, self.b, self.c])

#  

#        old_geom_center = points.sum(axis=0)/len(points)

#        # Translate points

#        points = points.copy() + self.position

#  

#        # Translate position to geometric_center along directions

#        # a and b

#        geom_center = points.sum(axis=0)/len(points)

#        self.position += self.a * np.dot(geom_center, self.a)

#        self.position += self.b * np.dot(geom_center, self.b)

#  

#        # Translate pivot to the geometric center

#        self.pivot = old_geom_center

#  

#        # Get the bounding sphere radius by searching for the most

#        # distant point

#        bound_radius = np.sqrt(((points-geom_center) * (points-geom_center)).sum(axis=1).max())

#  

#        # Calculate the distance in order to have the most distant

#        # point in our field of view (top/bottom)

#        fov_topbottom = self.FIELD_OF_VIEW_IN_DEGREES*np.pi/180.0

#  

#        dist = (bound_radius + self.z_near)/np.tan(fov_topbottom * 0.5)

#  

#        # Set the c-component of the position at the calculated distance

#        # 1) translate the position on the pivot

#        self.position = self.pivot.copy() 

#        # 2) add the distance plus a little extra room

#        self.position -= self.c * (dist*(1 + extraoff))

#  