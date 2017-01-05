import sys
from PyQt4 import QtGui, QtCore


class MainWindow(QtGui.QMainWindow):

    def __init__(self):
        super(MainWindow, self).__init__()

        self.do_something() #sanity check
        self.cw = ChildWidget(self)
        self.setCentralWidget(self.cw)
        self.show()

    def do_something(self):

        print 'doing something!'


class ChildWidget(QtGui.QWidget):

    def __init__(self, parent):
        super(ChildWidget, self).__init__(parent)

        self.button1 = QtGui.QPushButton()
        self.button1.clicked.connect(self.do_something_else)

        self.button2 = QtGui.QPushButton()
        self.button2.clicked.connect(self.parent().do_something)

        self.layout = QtGui.QVBoxLayout()
        self.layout.addWidget(self.button1)
        self.layout.addWidget(self.button2)
        self.setLayout(self.layout)
        self.show()

    def do_something_else(self):

        print 'doing something else!'


def main():
    app = QtGui.QApplication(sys.argv)
    ex = MainWindow()
    sys.exit(app.exec_())

if __name__ == '__main__':
    main()