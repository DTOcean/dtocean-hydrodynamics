# -*- coding: utf-8 -*-
"""
Created on Mon May 30 15:08:09 2016
@author: at

"""

from Camera import Camera
from MyGeom import Point3D, Vector3D, Matrix4x4
from PyQt4 import QtGui, QtCore
from PyQt4.QtOpenGL import *
from OpenGL.GL import *
from OpenGL.GLU import *
#import math
import sys
#import numpy as np
from mesh_class import Mesh
#import colors
#import time

class Viewer3DWidget(QGLWidget):

    def __init__(self, parent, path, file_name):

        QGLWidget.__init__(self, parent)
        self.setMouseTracking(True)
        self.setMinimumSize(800, 800)
        self.mesh_obj = Mesh(path, file_name)
        self.scaling_factor = self.get_scale()
        self.camera = Camera()
        self.camera.setSceneRadius( 2*self.scaling_factor )
        self.camera.reset()
        self.isPressed = False
        self.oldx = self.oldy = 0
    
        self.drawMesh_cond = True
        self.drawNorms_cond = False
        self.drawWire_cond = False
#        self.init_view()

    def get_scale(self):
        return max(abs(self.mesh_obj.v.max(0)).max(), 
                   abs(self.mesh_obj.v.min(0)).max())
         
#    def init_view(self):
#        new = 0
#        self.camera.orbit(self.oldx,self.oldy,1,-1)
#        self.update()
#        self.oldx = -new
#        self.oldy = 0

    def paintGL(self):
        glMatrixMode( GL_PROJECTION )
        glLoadIdentity()
        self.camera.transform()
        glMatrixMode( GL_MODELVIEW )
        glLoadIdentity()

        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)

        glDepthFunc( GL_LEQUAL )
        glEnable( GL_DEPTH_TEST )
        glEnable( GL_CULL_FACE )
        glFrontFace( GL_CCW )
        glDisable( GL_LIGHTING )
        glShadeModel( GL_FLAT )


        glColor(1.0, 1.0, 1.0)
        glColor(1.0, 0.0, 0.0)
        glBegin(GL_LINES)
        glVertex( 0, 0, 0)
        glVertex( 1*self.scaling_factor/4, 0, 0)
        glEnd()

        glColor(0.0, 1.0, 0.0)
        glBegin(GL_LINES)
        glVertex( 0, 0, 0)
        glVertex( 0, 1*self.scaling_factor/4, 0)
        glEnd()

        glColor(0.0, 0.0, 1.0)
        glBegin(GL_LINES)
        glVertex( 0, 0, 0)
        glVertex( 0, 0, 1*self.scaling_factor/4)
        glEnd()
        
        glColor(1.0, 1.0, 1.0)
        self.renderText(0, 0, 0,"(0,0,0)")
        self.renderText(self.scaling_factor, self.scaling_factor, 
                        self.scaling_factor,"(%s,%s,%s)"%(self.scaling_factor, 
                                                          self.scaling_factor, 
                                                          self.scaling_factor))
        self.renderText(-self.scaling_factor, self.scaling_factor, 
                        self.scaling_factor,"(%s,%s,%s)"%(-self.scaling_factor, 
                                                          self.scaling_factor, 
                                                          self.scaling_factor))
        self.renderText(self.scaling_factor, -self.scaling_factor, 
                        self.scaling_factor,"(%s,%s,%s)"%(self.scaling_factor, 
                                                          -self.scaling_factor, 
                                                          self.scaling_factor))
        glColor(0.0, 1.0, 0.0)
        self.bounding_box()

        if self.drawMesh_cond:
            self.drawMesh()

        if self.drawNorms_cond:
            self.drawNorms()

        if self.drawWire_cond:
            self.drawMesh(wireframe=True)



    def drawNorms(self):
        self.makeCurrent()
        ind = -1

        while self.mesh_obj.p.shape[0]-1 > ind:
            ind += 1
            centr = self.mesh_obj.c[ind]
            cpn = self.mesh_obj.n[ind]+centr
            glColor3f(1,0,0)
            glBegin(GL_LINES)
            glVertex3f(centr[0], centr[1],  centr[2])
            glVertex3f(cpn[0], cpn[1],  cpn[2])
            glEnd();

    

    def drawMesh(self, wireframe=False):

        self.makeCurrent()

        ind = -1

        while self.mesh_obj.p.shape[0]-1 > ind:

            ind += 1
            mesh = self.mesh_obj.v[self.mesh_obj.p[ind,:],:]

            if not wireframe:

                glColor3f(0.2,0.8,0.5)

                glBegin(GL_QUADS)
                glVertex3f(mesh[0][0], mesh[0][1],  mesh[0][2])
                glVertex3f(mesh[1][0], mesh[1][1],  mesh[1][2])
                glVertex3f(mesh[2][0], mesh[2][1],  mesh[2][2])
                glVertex3f(mesh[3][0], mesh[3][1],  mesh[3][2])
                glEnd();

                glColor3f(0.05,0.55,0.3)
                glBegin(GL_QUADS)
                glVertex3f(mesh[3][0], mesh[3][1],  mesh[3][2])
                glVertex3f(mesh[2][0], mesh[2][1],  mesh[2][2])
                glVertex3f(mesh[1][0], mesh[1][1],  mesh[1][2])
                glVertex3f(mesh[0][0], mesh[0][1],  mesh[0][2])
                glEnd();

            

            glColor3f(1.0,1.0,1.0)
            glLineWidth(2)
            glBegin(GL_LINES)
            glVertex3f(mesh[0][0], mesh[0][1],  mesh[0][2])
            glVertex3f(mesh[1][0], mesh[1][1],  mesh[1][2])
            glEnd();

            glBegin(GL_LINES)
            glVertex3f(mesh[1][0], mesh[1][1],  mesh[1][2])
            glVertex3f(mesh[2][0], mesh[2][1],  mesh[2][2])
            glEnd();

            glBegin(GL_LINES)
            glVertex3f(mesh[2][0], mesh[2][1],  mesh[2][2])
            glVertex3f(mesh[3][0], mesh[3][1],  mesh[3][2])
            glEnd();

            glBegin(GL_LINES)
            glVertex3f(mesh[3][0], mesh[3][1],  mesh[3][2])
            glVertex3f(mesh[0][0], mesh[0][1],  mesh[0][2])
            glEnd();
            
            
    def bounding_box(self):
        self.makeCurrent()

        sf = 2*self.scaling_factor

        glLineWidth(0.5)
        glColor3f(1.0,1.0,1.0)

        glBegin(GL_LINES)
        glVertex3f(-0.5*sf, -0.5*sf, -0.5*sf)
        glVertex3f(0.5*sf, -0.5*sf, -0.5*sf)
        glEnd();
        
        glBegin(GL_LINES)
        glVertex3f(0.5*sf, -0.5*sf, -0.5*sf)
        glVertex3f(0.5*sf, 0.5*sf, -0.5*sf)
        glEnd();

        glBegin(GL_LINES)
        glVertex3f(0.5*sf, 0.5*sf, -0.5*sf)
        glVertex3f(-0.5*sf, 0.5*sf, -0.5*sf)
        glEnd();

        glBegin(GL_LINES)        
        glVertex3f(-0.5*sf, 0.5*sf, -0.5*sf)
        glVertex3f(-0.5*sf, -0.5*sf, -0.5*sf)
        glEnd();

        glBegin(GL_LINES)
        glVertex3f(-0.5*sf, -0.5*sf, -0.5*sf)
        glVertex3f(-0.5*sf, -0.5*sf, 0.5*sf)
        glEnd();
        
        glBegin(GL_LINES)
        glVertex3f(0.5*sf, -0.5*sf, -0.5*sf)
        glVertex3f(0.5*sf, -0.5*sf, 0.5*sf)
        glEnd();

        glBegin(GL_LINES)
        glVertex3f(0.5*sf, 0.5*sf, -0.5*sf)
        glVertex3f(0.5*sf, 0.5*sf, 0.5*sf)
        glEnd();

        glBegin(GL_LINES)        
        glVertex3f(-0.5*sf, 0.5*sf, -0.5*sf)
        glVertex3f(-0.5*sf, 0.5*sf, 0.5*sf)
        glEnd();
        
        glBegin(GL_LINES)
        glVertex3f(-0.5*sf, -0.5*sf, 0.5*sf)
        glVertex3f(0.5*sf, -0.5*sf, 0.5*sf)
        glEnd();
        
        glBegin(GL_LINES)
        glVertex3f(0.5*sf, -0.5*sf, 0.5*sf)
        glVertex3f(0.5*sf, 0.5*sf, 0.5*sf)
        glEnd();

        glBegin(GL_LINES)
        glVertex3f(0.5*sf, 0.5*sf, 0.5*sf)
        glVertex3f(-0.5*sf, 0.5*sf, 0.5*sf)
        glEnd();

        glBegin(GL_LINES)        
        glVertex3f(-0.5*sf, 0.5*sf, 0.5*sf)
        glVertex3f(-0.5*sf, -0.5*sf, 0.5*sf)
        glEnd();


    def resizeGL(self, widthInPixels, heightInPixels):
        self.camera.setViewportDimensions(widthInPixels, heightInPixels)
        glViewport(0, 0, widthInPixels, heightInPixels)

    def initializeGL(self):
        self.qglClearColor(QtGui.QColor(0, 0,  0))
        glClearDepth(1.0)

    def mouseMoveEvent(self, mouseEvent):
        if int(mouseEvent.buttons()) != QtCore.Qt.NoButton :
            # user is dragging
            delta_x = mouseEvent.x() - self.oldx
            delta_y = self.oldy - mouseEvent.y()
            if int(mouseEvent.buttons()) & QtCore.Qt.RightButton :
                if int(mouseEvent.buttons()) & QtCore.Qt.MidButton :
                    self.camera.dollyCameraForward( 3*(delta_x+delta_y), False )
                else:
                    self.camera.orbit(self.oldx,self.oldy,mouseEvent.x(),mouseEvent.y())
            elif int(mouseEvent.buttons()) & QtCore.Qt.LeftButton :
                self.camera.translateSceneRightAndUp( delta_x, delta_y )
            self.update()
        self.oldx = mouseEvent.x()
        self.oldy = mouseEvent.y()

    def mouseDoubleClickEvent(self, mouseEvent):
        print("double click")
        self.showFullScreen()

    def mousePressEvent(self, e):
        print("mouse press")
        self.isPressed = True

    def mouseReleaseEvent(self, e):
        print("mouse release")
        self.isPressed = False
       
    def wheelEvent(self, e):
        z = e.delta()
        self.camera.mouse_zoom(z*0.001)
        self.update()


class PythonQtOpenGLMeshViewer(QtGui.QMainWindow):

    def __init__(self, path, file_name):

        QtGui.QMainWindow.__init__(self)
        self.setWindowTitle('Visualisation of the mesh: %s'%('/' + path + '/' + file_name))
        self.statusBar().showMessage("Mesh of the structure")

        exit = QtGui.QAction("Exit", self)
        exit.setShortcut("Ctrl+Q")
        exit.setStatusTip('Exit application')
        self.connect(exit, QtCore.SIGNAL('triggered()'), QtCore.SLOT('close()'))

        help = QtGui.QAction("Help", self)        
        help.setShortcut("F1")
        help.setStatusTip('Help')
        self.connect(help, QtCore.SIGNAL('triggered()'), QtCore.SLOT(''))

        menubar = self.menuBar()
        fileMenu = menubar.addMenu('&File')
        fileMenu.addAction(help)
        fileMenu.addAction(exit)

        fileMenu1 = menubar.addMenu('&About')

        about = QtGui.QAction("About", self)        
        about.setShortcut("Ctrl+A")
        about.setStatusTip('About')
        self.connect(about, QtCore.SIGNAL('triggered()'), QtCore.SLOT(''))
        fileMenu1.addAction(about)

        self.setToolTip('Mesh plotter')

        self.viewer3D = Viewer3DWidget(self, path, file_name)
        createButtons = True
        if createButtons:
            parentWidget = QtGui.QWidget()

            self.button1 = QtGui.QPushButton("Structure", self)
            self.button1.setStatusTip('Visualise the structure')
            self.viewer3D.makeCurrent()
            self.connect(self.button1, QtCore.SIGNAL('clicked()'), self.button1Action)

            self.button2 = QtGui.QPushButton("Norms", self)
            self.button2.setStatusTip('Visualise the norms')
            self.connect(self.button2, QtCore.SIGNAL('clicked()'), self.button2Action)

            self.button3 = QtGui.QPushButton("Wireframe", self)
            self.button3.setStatusTip('Visualise the the structure wireframe')
            self.connect(self.button3, QtCore.SIGNAL('clicked()'), self.button3Action)

            vbox = QtGui.QVBoxLayout()
            vbox.addWidget(self.button1)
            vbox.addWidget(self.button3)
            vbox.addWidget(self.button2)
            vbox.addStretch(1)
            self.viewer3D.setSizePolicy( QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Expanding )
            hbox = QtGui.QHBoxLayout()
            hbox.addLayout(vbox)
            hbox.addWidget(self.viewer3D)

            parentWidget.setLayout(hbox)
            self.setCentralWidget(parentWidget)
        
        else:
            self.setCentralWidget(self.viewer3D)
#        self.resize(800,800)

    def closeEvent(self, event):
        reply = QtGui.QMessageBox.question(self, "Confirmation",
            "Are you sure you want to quit?", QtGui.QMessageBox.Yes | QtGui.QMessageBox.No, QtGui.QMessageBox.No)
        if reply == QtGui.QMessageBox.Yes:
            event.accept()
        else:
            event.ignore()

    def button1Action(self):
        if self.viewer3D.drawMesh_cond:
            self.viewer3D.makeCurrent()
            self.viewer3D.drawMesh_cond = False
        else:
            self.viewer3D.makeCurrent()
            self.viewer3D.drawMesh_cond = True
            self.viewer3D.drawWire_cond = False
        self.viewer3D.paintGL()
        self.viewer3D.updateGL()
        
    def button2Action(self):
        if self.viewer3D.drawNorms_cond:
            self.viewer3D.makeCurrent()
            self.viewer3D.drawNorms_cond = False
        else:
            self.viewer3D.makeCurrent()
            self.viewer3D.drawNorms_cond = True
        self.viewer3D.paintGL()
        self.viewer3D.updateGL()

    def button3Action(self):
        if self.viewer3D.drawWire_cond:
            self.viewer3D.makeCurrent()
            self.viewer3D.drawWire_cond = False
        else:
            self.viewer3D.makeCurrent()
            self.viewer3D.drawWire_cond = True
            self.viewer3D.drawMesh_cond = False
        self.viewer3D.paintGL()
        self.viewer3D.updateGL()


if __name__ == '__main__':
    file_name = "mesh_test.GDF"
    path = ""
    app = QtGui.QApplication(sys.argv)
    window = PythonQtOpenGLMeshViewer(path, file_name)
    window.show()
    window.viewer3D.bounding_box()
    sys.exit(app.exec_())


