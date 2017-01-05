# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file '../designer/read_wamit_form.ui'
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
        Form.resize(552, 294)
        self.gridLayout = QtGui.QGridLayout(Form)
        self.gridLayout.setObjectName(_fromUtf8("gridLayout"))
        self.le_data_t4 = QtGui.QLineEdit(Form)
        self.le_data_t4.setObjectName(_fromUtf8("le_data_t4"))
        self.gridLayout.addWidget(self.le_data_t4, 0, 0, 1, 4)
        self.btn_browse_t4 = QtGui.QPushButton(Form)
        self.btn_browse_t4.setObjectName(_fromUtf8("btn_browse_t4"))
        self.gridLayout.addWidget(self.btn_browse_t4, 0, 4, 1, 1)
        self.mesh_f_t4 = QtGui.QLineEdit(Form)
        self.mesh_f_t4.setObjectName(_fromUtf8("mesh_f_t4"))
        self.gridLayout.addWidget(self.mesh_f_t4, 2, 0, 1, 2)
        self.btn_mesh_f_t4 = QtGui.QPushButton(Form)
        self.btn_mesh_f_t4.setObjectName(_fromUtf8("btn_mesh_f_t4"))
        self.gridLayout.addWidget(self.btn_mesh_f_t4, 2, 2, 1, 1)
        self.btn_view_mesh_t4 = QtGui.QPushButton(Form)
        self.btn_view_mesh_t4.setObjectName(_fromUtf8("btn_view_mesh_t4"))
        self.gridLayout.addWidget(self.btn_view_mesh_t4, 2, 3, 1, 2)
        self.btn_submit_t4 = QtGui.QPushButton(Form)
        self.btn_submit_t4.setObjectName(_fromUtf8("btn_submit_t4"))
        self.gridLayout.addWidget(self.btn_submit_t4, 3, 0, 1, 1)
        self.btn_load_data_t4 = QtGui.QPushButton(Form)
        self.btn_load_data_t4.setObjectName(_fromUtf8("btn_load_data_t4"))
        self.gridLayout.addWidget(self.btn_load_data_t4, 3, 1, 1, 4)
        self.gen_input_3 = QtGui.QGroupBox(Form)
        self.gen_input_3.setObjectName(_fromUtf8("gen_input_3"))
        self.gridLayout_18 = QtGui.QGridLayout(self.gen_input_3)
        self.gridLayout_18.setObjectName(_fromUtf8("gridLayout_18"))
        self.label_34 = QtGui.QLabel(self.gen_input_3)
        self.label_34.setObjectName(_fromUtf8("label_34"))
        self.gridLayout_18.addWidget(self.label_34, 2, 2, 1, 2)
        self.moor_dof_t4 = QtGui.QLineEdit(self.gen_input_3)
        self.moor_dof_t4.setText(_fromUtf8(""))
        self.moor_dof_t4.setObjectName(_fromUtf8("moor_dof_t4"))
        self.gridLayout_18.addWidget(self.moor_dof_t4, 2, 0, 1, 2)
        self.label_35 = QtGui.QLabel(self.gen_input_3)
        self.label_35.setObjectName(_fromUtf8("label_35"))
        self.gridLayout_18.addWidget(self.label_35, 0, 2, 1, 2)
        self.ndof_t4 = QtGui.QLineEdit(self.gen_input_3)
        self.ndof_t4.setText(_fromUtf8(""))
        self.ndof_t4.setObjectName(_fromUtf8("ndof_t4"))
        self.gridLayout_18.addWidget(self.ndof_t4, 0, 0, 1, 2)
        self.label_36 = QtGui.QLabel(self.gen_input_3)
        self.label_36.setObjectName(_fromUtf8("label_36"))
        self.gridLayout_18.addWidget(self.label_36, 1, 2, 1, 2)
        self.pto_dof_t4 = QtGui.QLineEdit(self.gen_input_3)
        self.pto_dof_t4.setText(_fromUtf8(""))
        self.pto_dof_t4.setObjectName(_fromUtf8("pto_dof_t4"))
        self.gridLayout_18.addWidget(self.pto_dof_t4, 1, 0, 1, 2)
        self.cb_gen_array_mat_t4 = QtGui.QCheckBox(self.gen_input_3)
        self.cb_gen_array_mat_t4.setObjectName(_fromUtf8("cb_gen_array_mat_t4"))
        self.gridLayout_18.addWidget(self.cb_gen_array_mat_t4, 3, 0, 1, 1)
        self.gridLayout.addWidget(self.gen_input_3, 1, 0, 1, 5)

        self.retranslateUi(Form)
        QtCore.QMetaObject.connectSlotsByName(Form)

    def retranslateUi(self, Form):
        Form.setWindowTitle(_translate("Form", "Hydrodynamic - Module", None))
        self.le_data_t4.setPlaceholderText(_translate("Form", "C:\\users\\johndoe\\desktop\\data_WAMIT", None))
        self.btn_browse_t4.setText(_translate("Form", "browse", None))
        self.mesh_f_t4.setPlaceholderText(_translate("Form", "C:\\user\\johndoe\\desktop\\cylinder.gdf", None))
        self.btn_mesh_f_t4.setText(_translate("Form", "browse", None))
        self.btn_view_mesh_t4.setText(_translate("Form", "view", None))
        self.btn_submit_t4.setText(_translate("Form", "Submit Inputs", None))
        self.btn_load_data_t4.setText(_translate("Form", "Load Data", None))
        self.gen_input_3.setTitle(_translate("Form", "general inputs", None))
        self.label_34.setText(_translate("Form", "ID of the independent degree of freedoms connected to the mooring ", None))
        self.moor_dof_t4.setPlaceholderText(_translate("Form", "1", None))
        self.label_35.setText(_translate("Form", "number of independent degree of freedoms", None))
        self.ndof_t4.setPlaceholderText(_translate("Form", "6", None))
        self.label_36.setText(_translate("Form", "ID of the independent degree of freedoms connected to the PTO ", None))
        self.pto_dof_t4.setPlaceholderText(_translate("Form", "1, 3", None))
        self.cb_gen_array_mat_t4.setText(_translate("Form", "load array interaction matrices", None))

