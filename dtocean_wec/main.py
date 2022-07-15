# -*- coding: utf-8 -*-

#    Copyright (C) 2016 Francesco Ferri, Mathew Topper
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
Created on Fri Jun 10 11:22:25 2016

.. moduleauthor:: Francesco Ferri <ff@civil.aau.dk>
.. moduleauthor:: Mathew Topper <mathew.topper@dataonlygreater.com>
"""

import os
import sys
import logging

from PyQt4.QtCore import *
from PyQt4.QtGui import *

from dtocean_hydro.configure import get_install_paths
import dtocean_wave.utils.hdf5_interface as h5i

from .form_utils import *
from .load_new import Ui_Form as Ui_LN
from .mdi_layout import Ui_MainWindow
from .new_selection import Ui_Form as Ui_NP
from .pfit_form import PowerPerformance
from .plot_form import Plotter
from .tab1 import ReadDb
from .tab2 import RunNemoh
from .tab3 import ReadNemoh
from .tab4 import ReadWamit
from .utils.mesh_plotter import PythonQtOpenGLMeshViewer


class QtHandler(logging.Handler):
    def __init__(self):
        logging.Handler.__init__(self)
    def emit(self, record):
        record = self.format(record)
        if record: XStream.stdout().write('%s\n'%record)
        # originally: XStream.stdout().write("{}\n".format(record))


logger = logging.getLogger(__name__)
handler = QtHandler()
handler.setFormatter(logging.Formatter("%(levelname)s: %(message)s"))
logger.addHandler(handler)
logger.setLevel(logging.DEBUG)


class XStream(QObject):
    _stdout = None
    _stderr = None
    messageWritten = pyqtSignal(str)
    def flush( self ):
        pass
    def fileno( self ):
        return -1
    def write( self, msg ):
        if ( not self.signalsBlocked() ):
            self.messageWritten.emit(unicode(msg))
    @staticmethod
    def stdout():
        if ( not XStream._stdout ):
            XStream._stdout = XStream()
            sys.stdout = XStream._stdout
        return XStream._stdout
    @staticmethod
    def stderr():
        if ( not XStream._stderr ):
            XStream._stderr = XStream()
            sys.stderr = XStream._stderr
        return XStream._stderr


class NewProject(QWidget, Ui_NP):
    def __init__(self, parent=None):
        QWidget.__init__(self, parent)
        self.setupUi(self)
 
class Entrance(QWidget, Ui_LN):
    def __init__(self, parent=None):
        QWidget.__init__(self, parent)
        self.setupUi(self)

class MainWindow(QMainWindow, Ui_MainWindow):
    count = 0

    def __init__(self):
        QMainWindow.__init__(self)
        self.setupUi(self)
        self.menubar.triggered[QAction].connect(self.windowaction)
        self._project_folder = r'C:\Users\francesco\Desktop\test_prj'
        self._project_name = 'test_prj'
        self.lw_area.hide()
        self.run_entrance()
        self._data = None
        self.setAttribute(Qt.WA_DeleteOnClose, True)
        
        XStream.stdout().messageWritten.connect( self.console.insertPlainText )
        XStream.stderr().messageWritten.connect( self.console.insertPlainText )
        
        self.lw_area.clicked.connect(self.on_treeView_clicked)
        
        # Get path to data dirs through configuration
        path_dict = get_install_paths()
        
        self.bin_path = path_dict["bin_path"]
        self.wec_share_path = path_dict["wec_share_path"]
        
        return
        
    @pyqtSlot(QModelIndex)
    def on_treeView_clicked(self, index):
        print 'selected item index found at %s with data: %s' % (index.row(), index.data())
        if index.data() == 'Hydrodynamic':
            try:
                self.form_hyd.isVisible()
            except:
                self.reopen_form_hyd()
                
        elif index.data() == 'Performance Fit':
            try:
                self.form_power.isVisible()
            except:
                self.reopen_form_power()
        elif index.data() == 'Data Visualisation':
            try:
                self.form_plot.isVisible()
            except:
                self.reopen_form_plot()
        else:
            pass
     
    def task_show_mesh(self, path):
        self.mesh_view = PythonQtOpenGLMeshViewer(path['path'], path['f_n'])
        self.mesh_view.show()
        
    def task_reset_forms(self):
        # TODO: implement the clear method on data changed
        print("Not implemented yet")
    
    def task_clear_plot(self):
        self.form_plot.setEnabled(False)
        self.list_model.item(2).setForeground(QBrush(QColor(0,0,0)))
        
    def task_clear_power(self):
        if 'p_fit' in self._data.keys():
            del self._data['p_fit'] 
        self.reopen_form_power()
        self.form_power.setEnabled(False)
        self.list_model.item(1).setForeground(QBrush(QColor(0,0,0)))
        
    def reopen_form_hyd(self):
        prj_id = self._data['inputs_hydrodynamic']['general_input']['input_type']
        del(self.form_hyd)
        if prj_id == 1:
            self.form_hyd = ReadDb(self)
        elif prj_id == 2:
            self.form_hyd = RunNemoh(self)
        elif prj_id == 3:
            self.form_hyd = ReadNemoh(self)
        elif prj_id == 4:
            self.form_hyd = ReadWamit(self)
        
        self.mdi_area.addSubWindow(self.form_hyd)
        self.form_hyd.show()
        self.form_hyd.populate_prj()
        
    def reopen_form_power(self ):
        del(self.form_power)
        self.form_power = PowerPerformance(self)
        self.mdi_area.addSubWindow(self.form_power)
        self.form_power.show()
        ret = self.form_power.set_data(self._data)
        if ret == 0:
            self.list_model.item(1).setForeground(QBrush(QColor(0,0,150)))
            self.list_model.item(2).setForeground(QBrush(QColor(0,0,150)))
            self.actionGenerate_array_hydrodynamic.setEnabled(False)
        elif ret == 1:
            self.list_model.item(1).setForeground(QBrush(QColor(0,150,0)))
            self.list_model.item(2).setForeground(QBrush(QColor(0,150,0)))
            if self.form_hyd._data['inputs_hydrodynamic']['general_input']['get_array_mat'] == 1:
                self.actionGenerate_array_hydrodynamic.setEnabled(True)
        else:
            self.list_model.item(1).setForeground(QBrush(QColor(0,0,0)))
            self.list_model.item(2).setForeground(QBrush(QColor(0,0,0)))
            self.actionGenerate_array_hydrodynamic.setEnabled(False)
        
        
    def reopen_form_plot(self):
        del(self.form_plot)
        self.form_plot = Plotter(self)
        self.mdi_area.addSubWindow(self.form_plot)
        self.form_plot.show()
        self.form_plot.set_data(self._data)
        

    def pull_data_from_child(self):
        try:
            self._data['inputs_hydrodynamic'] = self.form_hyd._data['inputs_hydrodynamic']
        except:
            print("No hydrodynamic inputs data")
        
        try:
            self._data['inputs_pm_fit'] = self.form_power._data['inputs_pm_fit']
        except:
            print("No performance fit inputs data")
        
        try:
            self._data['hyd'] = self.form_hyd._data['hyd']
        except:
            print("No hydrodynamic results")
        
        
        try:
            self._data['p_fit'] = self.form_power._data['p_fit']
        except:
            print("No performance fit results")
        

    def windowaction(self, action):
        print "triggered"
        if action.text() == "Quit":
            self.close()
		
        if action.text() == "New Project":
            if self._data is None:
                self.gen_new_project()
            else:
                self.new_project_choice()
    		
        if action.text() == "Cascade":
           self.mdi_area.cascadeSubWindows()
    		
        if action.text() == "Tiled":
           self.mdi_area.tileSubWindows()
           
        if action.text() == "Generate array hydrodynamic":
            self.save_dtocean_format()
            
        if action.text() == "Save Project":
            if not self._data is None:
                self.pull_data_from_child()
                self.save_project()
        
        if action.text() == "Load Project":
            if self._data is None:
                self.load_project()
            else:
                self.load_project_choice()
    
    def save_choice(self):
        msgBox = QMessageBox()
        msgBox.setText("The project contains unsaved data.")
        msgBox.setInformativeText("Do you want to save your changes before to close the project?")
        msgBox.setStandardButtons(QMessageBox.Save | QMessageBox.Discard | QMessageBox.Cancel)
        msgBox.setDefaultButton(QMessageBox.Save)
        ret = msgBox.exec_()
        return ret
        
    def load_project_choice(self):
        ret = self.save_choice()
        
        if ret == QMessageBox.Save:
            # should never be reached
            self.save_project()
            self.clear_project_data()
            self.load_project()
        elif ret == QMessageBox.Discard:
            self.clear_project_data()
            self.load_project()
        elif ret == QMessageBox.Cancel:
            print("No project has been loaded")
            pass
        else:
            print("No project has been loaded")
            pass
        
    def closeEvent(self, event):
        choice = QMessageBox.question(self, 'Quit', 'Do you want to exit the program?',
                                            QMessageBox.Yes | QMessageBox.No)
        if choice == QMessageBox.Yes:
            event.accept()
        else:
            event.ignore()       
        
    def new_project_choice(self):
        ret = self.save_choice()
        
        if ret == QMessageBox.Save:
            # should never be reached
            prj_folder = self._data['prj_folder']
            prj_fn = self._data['prj_filename']
            project = os.path.join(prj_folder, prj_fn)
            h5i.save_dict_to_hdf5(self._data, project)
            self.gen_new_project()
        elif ret == QMessageBox.Discard:
            self.gen_new_project()
        elif ret == QMessageBox.Cancel:
            print("No new project has been generated")
            pass
        else:
            print("No new project has been generated")
            pass
        
    def task_save_hydrodynamic(self, data):
        self._data['inputs_hydrodynamic'] = data
        self.save_project()
    
    def save_project(self):
        prj_folder = self._data['prj_folder']
        prj_fn = self._data['prj_filename']
        project = os.path.join(prj_folder, prj_fn)
        h5i.save_dict_to_hdf5(self._data, project)
        
    def task_end_pfit(self, data):
        self._data['p_fit'] = data
        self.save_project()

        try:        
            self.form_plot.set_data(self._data)
        except:
            self.reopen_form_plot()
            
        self.list_model.item(1).setForeground(QBrush(QColor(0,150,0)))
        self.list_model.item(2).setForeground(QBrush(QColor(0,150,0)))
        if self.form_hyd._data['inputs_hydrodynamic']['general_input']['get_array_mat'] == 1:
                self.actionGenerate_array_hydrodynamic.setEnabled(True)
            
    def set_pfit_data(self, data):
        self._data['inputs_pm_fit'] = data
        self.save_project()
        
    def set_wec_db_results(self, data):
        self._data = data
        self._data['inputs_pm_fit']['data_folder'] = self._data['inputs_hydrodynamic']['general_input']['data_folder']
        self.set_hydrodynamic_results(self._data['hyd'])

    def set_hydrodynamic_results(self, data):
        self._data['hyd'] = data
        self.save_project()
        
        try:        
            self.form_plot.setEnabled(True)
        except:
            self.reopen_form_plot()
        try:        
            self.form_power.setEnabled(True)
        except:
            self.reopen_form_power()
            
        self.list_model.item(0).setForeground(QBrush(QColor(0,150,0)))
        self.list_model.item(1).setForeground(QBrush(QColor(0,0,150)))
        self.list_model.item(2).setForeground(QBrush(QColor(0,0,150)))
        self.form_plot.set_data(self._data)
        self.form_power.set_data(self._data)
        
        
    def run_entrance(self):
        self.entrance = Entrance(self)
        self.mdi_area.addSubWindow(self.entrance)
        self.entrance.show()
        self.entrance.btn_load.clicked.connect(self.load_project)
        self.entrance.btn_new.clicked.connect(self.gen_new_project)
        
        
    def clear_project_data(self):
        if not self._data is None:
            self.list_model.item(0).setForeground(QBrush(QColor(0,0,150)))
            self.list_model.item(1).setForeground(QBrush(QColor(0,0,0)))
            self.list_model.item(2).setForeground(QBrush(QColor(0,0,0)))
        
        self._data = None
        
    def set_forms_data(self):
        self.form_hyd.set_data(self._data)
        
        #self.form_power.set_data(self._data)
        #self.form_plot.set_data(self._data)
    
    def load_project(self):
        self.mdi_area.closeAllSubWindows()
        self.clear_project_data()
        project_folder = QFileDialog.getExistingDirectory(self,
                                                          "Select a project folder")
        project_folder = str(project_folder)        
        prj_name = os.path.basename(os.path.normpath(project_folder))
        prj_filename = "".join([prj_name, '_data_collection.hdf5'])
        project = os.path.join(project_folder, prj_filename)
        
        if os.path.isfile(project):
            dictionary = h5i.load_dict_from_hdf5(project)
            dictionary['prj_folder'] = project_folder
            dictionary['prj_filename'] = prj_filename
            dictionary['prj_name'] = prj_name
            self._data = dictionary
            self.populate_project()
        else:
            choice = QMessageBox().question(self, 'Wrong Project Folder',
                                                    'The selected folder does not contains a valid project data collection.\n Do you want to select another folder?',
                                                    QMessageBox.Yes | QMessageBox.No) 
            if choice == QMessageBox.Yes:
                self.load_project()
            else:
                return -1
    
    def populate_project(self):
        prj_id = self._data['inputs_hydrodynamic']['general_input']['input_type']
        if prj_id == 1:
            self.create_project(1)
        elif prj_id == 2:
            self.create_project(2)
        elif prj_id == 3:
            self.create_project(3)
        elif prj_id == 4:
            self.create_project(4)
        else:
            print("Input type not understood")
            
        # self.form_hyd.populate_prj()
        
      
    
    def gen_new_project(self):
        self.mdi_area.closeAllSubWindows()
        self.clear_project_data()
        self.new_project = NewProject(self)
        self.mdi_area.addSubWindow(self.new_project)
        self.new_project.show()
        
        
        self.new_project.btn_t1.clicked.connect(lambda: self.create_project(1))
        self.new_project.btn_t2.clicked.connect(lambda: self.create_project(2))
        self.new_project.btn_t3.clicked.connect(lambda: self.create_project(3))
        self.new_project.btn_t4.clicked.connect(lambda: self.create_project(4))
        self.new_project.btn_prj.clicked.connect(self.browse_prj)
     
    def browse_prj(self):
        folder = QFileDialog.getExistingDirectory(self, "Select a project folder")
        folder = str(folder)
        prj_name = os.path.basename(os.path.normpath(folder))
        if os.path.isdir(folder):
            if os.listdir(folder):
                msgBox = QMessageBox()
                msgBox.setText('The project folder is not empty \n Select another folder or clear the actual folder?')
                msgBox.addButton(QPushButton('Select'), QMessageBox.YesRole)
                msgBox.addButton(QPushButton('Clear'), QMessageBox.NoRole)
                msgBox.addButton(QPushButton('Cancel'), QMessageBox.RejectRole)
                ret = msgBox.exec_();
                if ret == 0:
                    self.browse_prj()
                elif ret == 1:
                    try:
                        clean_prj_folder(folder)
                    except:
                        return -1
                else:
                    return -1
        else:
            return -1
        self.new_project.le_prj.setText(folder)
        self._data = h5i.create_empty_project(folder, prj_name)
    
    def populate_listview(self):
        self.lw_area.show()
        model = QStandardItemModel()
        ls_item = [
            'Hydrodynamic', # Must be store-bought
            'Performance Fit', # Must be homemade
            'Data Visualisation', # Must be saucy
        ]
         
        for it in ls_item:
            # Create an item with a caption
            item = QStandardItem(it)
            # Add a checkbox to it
            item.setCheckable(False)
            item.setEditable(False)
            item.setForeground(QBrush(QColor(0,0,0)))
            # Add the item to the model
            model.appendRow(item)
        
        model.item(0).setForeground(QBrush(QColor(0,0,150)))
        self.lw_area.setModel(model)
        self.list_model = model
        
    
        
    def set_full_project(self):
        self.form_plot = Plotter(self)
        self.form_power = PowerPerformance(self)
        self.mdi_area.addSubWindow(self.form_power)
        self.mdi_area.addSubWindow(self.form_plot)
        self.form_plot.show()
        self.form_power.show()
        
        self.form_plot.setEnabled(False)
        self.form_power.setEnabled(False)
        
    def create_project(self, prj_id):
        if self._data is None:
             self.new_project.le_prj.setStyleSheet("QLineEdit {border: 1px solid red}")
        else:
            self.mdi_area.closeActiveSubWindow()
            self.populate_listview()
            self.set_full_project()
            MainWindow.count = MainWindow.count+1
            if prj_id == 1:
                self.form_hyd = ReadDb(self)
            elif prj_id == 2:
                self.form_hyd = RunNemoh(self)
            elif prj_id == 3:
                self.form_hyd = ReadNemoh(self)
            elif prj_id == 4:
                self.form_hyd = ReadWamit(self)
            else:
                print("Input type not understood")
            
            self.mdi_area.addSubWindow(self.form_hyd)
            self.form_hyd.show()
            
            self.mdi_area.tileSubWindows()
            self.actionSave_Project.setEnabled(True)
            self.set_forms_data()

    def save_dtocean_format(self):
        status = final_data_check(self._data)
        if status[0]:
            x = self._data['hyd']
            y = self._data['p_fit']
            
            data2dtocean = dict(x.items() + y.items() + [(k, x[k] + y[k]) for k in set(x) & set(y)])
            #update_wec_power_matrix(data2dtocean)
            h5i.save_dict_to_hdf5(data2dtocean, os.path.join(self._data['prj_folder'], 'wec_solution.h5'))
        else:
            print('The following fields are missing from the data attribute')
            print(' , \n'.join(status[1]))

    
    
def main():
   app = QApplication(sys.argv)
   ex = MainWindow()
   ex.show()
   sys.exit(app.exec_())
	
if __name__ == '__main__':
   main()