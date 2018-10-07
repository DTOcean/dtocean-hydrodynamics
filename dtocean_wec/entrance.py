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
Created on Tue May 31 19:17:18 2016

.. moduleauthor:: Francesco Ferri <ff@civil.aau.dk>
"""

import sys
from PyQt4 import QtGui, QtCore

class Window(QtGui.QMainWindow):

    def __init__(self):
        super(Window, self).__init__()
        self.setGeometry(50, 50, 500, 300)
        self.setWindowTitle("PyQT tuts!")
        self.setWindowIcon(QtGui.QIcon('pythonlogo.png'))
        self.home()

    def home(self):
        btn1 = QtGui.QPushButton("1", self)
        btn1.clicked.connect(QtCore.QCoreApplication.instance().quit)
        btn1.resize(100,100)
        btn1.move(0,0)
        
        btn2 = QtGui.QPushButton("2", self)
        btn2.clicked.connect(QtCore.QCoreApplication.instance().quit)
        btn2.resize(100,100)
        btn2.move(0,100)
        
        btn3 = QtGui.QPushButton("3", self)
        btn3.clicked.connect(QtCore.QCoreApplication.instance().quit)
        btn3.resize(100,100)
        btn3.move(0,200)
        
        btn4 = QtGui.QPushButton("4", self)
        btn4.clicked.connect(QtCore.QCoreApplication.instance().quit)
        btn4.resize(100,100)
        btn4.move(100,0)
        
        
        btn5 = QtGui.QPushButton("5", self)
        btn5.clicked.connect(QtCore.QCoreApplication.instance().quit)
        btn5.resize(100,100)
        btn5.move(100,100)
        
        btn6 = QtGui.QPushButton("6", self)
        btn6.clicked.connect(QtCore.QCoreApplication.instance().quit)
        btn6.resize(100,100)
        btn6.move(100,200)
                
        
        self.show()

        
def run():
    app = QtGui.QApplication(sys.argv)
    GUI = Window()
    sys.exit(app.exec_())

run()