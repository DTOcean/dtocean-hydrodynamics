# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file '../designer/read_db_form.ui'
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
        Form.resize(488, 635)
        self.gridLayout = QtGui.QGridLayout(Form)
        self.gridLayout.setObjectName(_fromUtf8("gridLayout"))
        self.label = QtGui.QLabel(Form)
        self.label.setMinimumSize(QtCore.QSize(0, 50))
        self.label.setMaximumSize(QtCore.QSize(188888, 50))
        font = QtGui.QFont()
        font.setPointSize(40)
        self.label.setFont(font)
        self.label.setObjectName(_fromUtf8("label"))
        self.gridLayout.addWidget(self.label, 1, 0, 1, 1)
        self.btn_load_data_t1 = QtGui.QPushButton(Form)
        self.btn_load_data_t1.setObjectName(_fromUtf8("btn_load_data_t1"))
        self.gridLayout.addWidget(self.btn_load_data_t1, 0, 1, 1, 1)
        self.cb_wec_t1 = QtGui.QComboBox(Form)
        self.cb_wec_t1.setObjectName(_fromUtf8("cb_wec_t1"))
        self.cb_wec_t1.addItem(_fromUtf8(""))
        self.cb_wec_t1.addItem(_fromUtf8(""))
        self.gridLayout.addWidget(self.cb_wec_t1, 0, 0, 1, 1)
        self.ww_wec = QtWebKit.QWebView(Form)
        self.ww_wec.setUrl(QtCore.QUrl(_fromUtf8("about:blank")))
        self.ww_wec.setObjectName(_fromUtf8("ww_wec"))
        self.gridLayout.addWidget(self.ww_wec, 3, 0, 1, 3)
        self.line = QtGui.QFrame(Form)
        self.line.setFrameShape(QtGui.QFrame.HLine)
        self.line.setFrameShadow(QtGui.QFrame.Sunken)
        self.line.setObjectName(_fromUtf8("line"))
        self.gridLayout.addWidget(self.line, 2, 0, 1, 3)

        self.retranslateUi(Form)
        QtCore.QMetaObject.connectSlotsByName(Form)

    def retranslateUi(self, Form):
        Form.setWindowTitle(_translate("Form", "Hydrodynamic - Module", None))
        self.label.setText(_translate("Form", "Case description", None))
        self.btn_load_data_t1.setText(_translate("Form", "load", None))
        self.cb_wec_t1.setItemText(0, _translate("Form", "heaving cylinder", None))
        self.cb_wec_t1.setItemText(1, _translate("Form", "surging barge", None))

from PyQt4 import QtWebKit
