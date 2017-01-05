# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file '../designer/new_selection.ui'
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

class Ui_Form(object):
    def setupUi(self, Form):
        Form.setObjectName(_fromUtf8("Form"))
        Form.resize(424, 524)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(Form.sizePolicy().hasHeightForWidth())
        Form.setSizePolicy(sizePolicy)
        Form.setMinimumSize(QtCore.QSize(424, 524))
        Form.setMaximumSize(QtCore.QSize(424, 524))
        self.gridLayout = QtGui.QGridLayout(Form)
        self.gridLayout.setObjectName(_fromUtf8("gridLayout"))
        self.btn_prj = QtGui.QPushButton(Form)
        self.btn_prj.setObjectName(_fromUtf8("btn_prj"))
        self.gridLayout.addWidget(self.btn_prj, 1, 2, 1, 1)
        self.btn_t4 = QtGui.QPushButton(Form)
        self.btn_t4.setMinimumSize(QtCore.QSize(200, 200))
        self.btn_t4.setObjectName(_fromUtf8("btn_t4"))
        self.gridLayout.addWidget(self.btn_t4, 4, 2, 1, 1)
        self.btn_t2 = QtGui.QPushButton(Form)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.btn_t2.sizePolicy().hasHeightForWidth())
        self.btn_t2.setSizePolicy(sizePolicy)
        self.btn_t2.setMinimumSize(QtCore.QSize(200, 200))
        self.btn_t2.setObjectName(_fromUtf8("btn_t2"))
        self.gridLayout.addWidget(self.btn_t2, 3, 2, 1, 1)
        self.btn_t1 = QtGui.QPushButton(Form)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.btn_t1.sizePolicy().hasHeightForWidth())
        self.btn_t1.setSizePolicy(sizePolicy)
        self.btn_t1.setMinimumSize(QtCore.QSize(200, 200))
        self.btn_t1.setObjectName(_fromUtf8("btn_t1"))
        self.gridLayout.addWidget(self.btn_t1, 3, 1, 1, 1)
        self.btn_t3 = QtGui.QPushButton(Form)
        sizePolicy = QtGui.QSizePolicy(QtGui.QSizePolicy.Fixed, QtGui.QSizePolicy.Fixed)
        sizePolicy.setHorizontalStretch(0)
        sizePolicy.setVerticalStretch(0)
        sizePolicy.setHeightForWidth(self.btn_t3.sizePolicy().hasHeightForWidth())
        self.btn_t3.setSizePolicy(sizePolicy)
        self.btn_t3.setMinimumSize(QtCore.QSize(200, 200))
        self.btn_t3.setObjectName(_fromUtf8("btn_t3"))
        self.gridLayout.addWidget(self.btn_t3, 4, 1, 1, 1)
        self.label = QtGui.QLabel(Form)
        self.label.setObjectName(_fromUtf8("label"))
        self.gridLayout.addWidget(self.label, 0, 1, 1, 1)
        self.le_prj = QtGui.QLineEdit(Form)
        self.le_prj.setObjectName(_fromUtf8("le_prj"))
        self.gridLayout.addWidget(self.le_prj, 1, 1, 1, 1)
        self.label_2 = QtGui.QLabel(Form)
        self.label_2.setObjectName(_fromUtf8("label_2"))
        self.gridLayout.addWidget(self.label_2, 2, 1, 1, 1)

        self.retranslateUi(Form)
        QtCore.QMetaObject.connectSlotsByName(Form)

    def retranslateUi(self, Form):
        Form.setWindowTitle(_translate("Form", "WEC - Analysis / New Project Selection", None))
        self.btn_prj.setText(_translate("Form", "Browse", None))
        self.btn_t4.setText(_translate("Form", "Load WAMIT", None))
        self.btn_t2.setText(_translate("Form", "Run Nemoh", None))
        self.btn_t1.setText(_translate("Form", "WEC DB", None))
        self.btn_t3.setText(_translate("Form", "Load Nemoh", None))
        self.label.setText(_translate("Form", "Porject Folder", None))
        self.le_prj.setPlaceholderText(_translate("Form", "C:\\Users\\JohnDoe\\test_project", None))
        self.label_2.setText(_translate("Form", "Project Type", None))

