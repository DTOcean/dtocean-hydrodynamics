# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file '../designer/mdi_layout.ui'
#
# Created by: PyQt4 UI code generator 4.11.4
#
# WARNING! All changes made in this file will be lost!

from PyQt4 import QtCore, QtGui

try:
    _fromUtf8 = QtCore.QString.fromUtf8
except AttributeError:
    def _fromUtf8(s):
        return s

try:
    _encoding = QtGui.QApplication.UnicodeUTF8
    def _translate(context, text, disambig):
        return QtGui.QApplication.translate(context, text, disambig, _encoding)
except AttributeError:
    def _translate(context, text, disambig):
        return QtGui.QApplication.translate(context, text, disambig)

class Ui_MainWindow(object):
    def setupUi(self, MainWindow):
        MainWindow.setObjectName(_fromUtf8("MainWindow"))
        MainWindow.resize(1000, 800)
        self.centralwidget = QtGui.QWidget(MainWindow)
        self.centralwidget.setObjectName(_fromUtf8("centralwidget"))
        self.gridLayout = QtGui.QGridLayout(self.centralwidget)
        self.gridLayout.setObjectName(_fromUtf8("gridLayout"))
        self.mdi_area = QtGui.QMdiArea(self.centralwidget)
        self.mdi_area.setObjectName(_fromUtf8("mdi_area"))
        self.gridLayout.addWidget(self.mdi_area, 0, 1, 1, 1)
        self.lw_area = QtGui.QListView(self.centralwidget)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Expanding)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.lw_area.sizePolicy().hasHeightForWidth())
        self.lw_area.setSizePolicy(sizePolicy)
        self.lw_area.setMinimumSize(QtCore.QSize(200, 0))
        self.lw_area.setMaximumSize(QtCore.QSize(200, 16777215))
        self.lw_area.setObjectName(_fromUtf8("lw_area"))
        self.gridLayout.addWidget(self.lw_area, 0, 0, 1, 1)
        self.console = QtGui.QTextEdit(self.centralwidget)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.console.sizePolicy().hasHeightForWidth())
        self.console.setSizePolicy(sizePolicy)
        self.console.setObjectName(_fromUtf8("console"))
        self.gridLayout.addWidget(self.console, 1, 0, 1, 2)
        MainWindow.setCentralWidget(self.centralwidget)
        self.menubar = QtGui.QMenuBar(MainWindow)
        self.menubar.setGeometry(QtCore.QRect(0, 0, 1000, 21))
        self.menubar.setObjectName(_fromUtf8("menubar"))
        self.menuFile = QtGui.QMenu(self.menubar)
        self.menuFile.setObjectName(_fromUtf8("menuFile"))
        self.menuDTOcean = QtGui.QMenu(self.menubar)
        self.menuDTOcean.setObjectName(_fromUtf8("menuDTOcean"))
        self.menuWindows = QtGui.QMenu(self.menubar)
        self.menuWindows.setObjectName(_fromUtf8("menuWindows"))
        self.menuOpen_Window = QtGui.QMenu(self.menuWindows)
        self.menuOpen_Window.setObjectName(_fromUtf8("menuOpen_Window"))
        MainWindow.setMenuBar(self.menubar)
        self.statusbar = QtGui.QStatusBar(MainWindow)
        self.statusbar.setObjectName(_fromUtf8("statusbar"))
        MainWindow.setStatusBar(self.statusbar)
        self.actionNew_Project = QtGui.QAction(MainWindow)
        self.actionNew_Project.setObjectName(_fromUtf8("actionNew_Project"))
        self.actionLoad_Project = QtGui.QAction(MainWindow)
        self.actionLoad_Project.setObjectName(_fromUtf8("actionLoad_Project"))
        self.actionAbout = QtGui.QAction(MainWindow)
        self.actionAbout.setObjectName(_fromUtf8("actionAbout"))
        self.actionQuit = QtGui.QAction(MainWindow)
        self.actionQuit.setObjectName(_fromUtf8("actionQuit"))
        self.actionAbout_2 = QtGui.QAction(MainWindow)
        self.actionAbout_2.setObjectName(_fromUtf8("actionAbout_2"))
        self.actionQuit_2 = QtGui.QAction(MainWindow)
        self.actionQuit_2.setObjectName(_fromUtf8("actionQuit_2"))
        self.actionGenerate_array_hydrodynamic = QtGui.QAction(MainWindow)
        self.actionGenerate_array_hydrodynamic.setEnabled(False)
        self.actionGenerate_array_hydrodynamic.setObjectName(_fromUtf8("actionGenerate_array_hydrodynamic"))
        self.actionCascade = QtGui.QAction(MainWindow)
        self.actionCascade.setObjectName(_fromUtf8("actionCascade"))
        self.actionTiled = QtGui.QAction(MainWindow)
        self.actionTiled.setObjectName(_fromUtf8("actionTiled"))
        self.actionMinimise_All = QtGui.QAction(MainWindow)
        self.actionMinimise_All.setObjectName(_fromUtf8("actionMinimise_All"))
        self.actionSave_Project = QtGui.QAction(MainWindow)
        self.actionSave_Project.setEnabled(False)
        self.actionSave_Project.setObjectName(_fromUtf8("actionSave_Project"))
        self.actionHydrodynamic = QtGui.QAction(MainWindow)
        self.actionHydrodynamic.setObjectName(_fromUtf8("actionHydrodynamic"))
        self.actionPerformance_Fit = QtGui.QAction(MainWindow)
        self.actionPerformance_Fit.setObjectName(_fromUtf8("actionPerformance_Fit"))
        self.actionData_Visualisation = QtGui.QAction(MainWindow)
        self.actionData_Visualisation.setObjectName(_fromUtf8("actionData_Visualisation"))
        self.menuFile.addAction(self.actionNew_Project)
        self.menuFile.addAction(self.actionLoad_Project)
        self.menuFile.addAction(self.actionSave_Project)
        self.menuFile.addSeparator()
        self.menuFile.addAction(self.actionAbout_2)
        self.menuFile.addAction(self.actionQuit_2)
        self.menuDTOcean.addAction(self.actionGenerate_array_hydrodynamic)
        self.menuOpen_Window.addAction(self.actionHydrodynamic)
        self.menuOpen_Window.addAction(self.actionPerformance_Fit)
        self.menuOpen_Window.addAction(self.actionData_Visualisation)
        self.menuWindows.addAction(self.actionCascade)
        self.menuWindows.addAction(self.actionTiled)
        self.menuWindows.addSeparator()
        self.menuWindows.addAction(self.menuOpen_Window.menuAction())
        self.menubar.addAction(self.menuFile.menuAction())
        self.menubar.addAction(self.menuDTOcean.menuAction())
        self.menubar.addAction(self.menuWindows.menuAction())

        self.retranslateUi(MainWindow)
        QtCore.QObject.connect(self.actionQuit_2, QtCore.SIGNAL(_fromUtf8("triggered()")), MainWindow.close)
        QtCore.QMetaObject.connectSlotsByName(MainWindow)

    def retranslateUi(self, MainWindow):
        MainWindow.setWindowTitle(_translate("MainWindow", "WEC Analysis", None))
        self.menuFile.setTitle(_translate("MainWindow", "File", None))
        self.menuDTOcean.setTitle(_translate("MainWindow", "DTOcean", None))
        self.menuWindows.setTitle(_translate("MainWindow", "Windows", None))
        self.menuOpen_Window.setTitle(_translate("MainWindow", "Open Window", None))
        self.actionNew_Project.setText(_translate("MainWindow", "New Project", None))
        self.actionLoad_Project.setText(_translate("MainWindow", "Load Project", None))
        self.actionAbout.setText(_translate("MainWindow", "About", None))
        self.actionQuit.setText(_translate("MainWindow", "Quit", None))
        self.actionAbout_2.setText(_translate("MainWindow", "About", None))
        self.actionQuit_2.setText(_translate("MainWindow", "Quit", None))
        self.actionGenerate_array_hydrodynamic.setText(_translate("MainWindow", "Generate array hydrodynamic", None))
        self.actionCascade.setText(_translate("MainWindow", "Cascade", None))
        self.actionTiled.setText(_translate("MainWindow", "Tiled", None))
        self.actionMinimise_All.setText(_translate("MainWindow", "Minimise All", None))
        self.actionSave_Project.setText(_translate("MainWindow", "Save Project", None))
        self.actionHydrodynamic.setText(_translate("MainWindow", "Hydrodynamic", None))
        self.actionPerformance_Fit.setText(_translate("MainWindow", "Performance Fit", None))
        self.actionData_Visualisation.setText(_translate("MainWindow", "Data Visualisation", None))

