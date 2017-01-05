import os
import sys
import argparse

from PyQt4 import QtGui
from main import MainWindow

module_path = os.path.realpath(__file__)

def main(debug=False):

    app = QtGui.QApplication(sys.argv)
    ex = MainWindow()
    ex.show()
    sys.exit(app.exec_())
    
    return
    
def gui_interface():

    '''Command line interface for dtocean-gui.
    
    Example:
    
        To get help::
        
            $ dtocean-wec -h
            
    '''
    
    epiStr = ('''AAU (c) 2016.''')
              
    desStr = "Run the WEC standalone GUI."

    parser = argparse.ArgumentParser(description=desStr,
                                     epilog=epiStr)

    parser.add_argument("--debug",
                        help=("disable stream redirection"),
                        action='store_true')
                                     
    args = parser.parse_args()
    debug = args.debug
        
    main(debug)

    return

