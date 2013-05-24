import sys
from PyQt4 import QtGui, QtCore
from PyQt4.QtCore import *
from PyQt4.QtGui import *

class ModelViewer(QtGui.QWidget):
    def __init__(self, parent=None):
        QtGui.QWidget.__init__(self, parent)

        self.setGeometry(1000, 1000, 200, 80)
        self.setWindowTitle('Data Table')

        # quit = QtGui.QPushButton('Close', self)
        # quit.setGeometry(10, 10, 60, 35)
        #
        # self.connect(quit, QtCore.SIGNAL('clicked()'),
        #              self, QtCore.SLOT('close()'))

    def show_model(self, model):
        """
        :param Radial1DModel: object with attribute t_rads,
        :return: Affective method, opens window to display table
        """
        tablemodel = MyTableModel(self)
        tablemodel.add_data(model)
        tableview = QtGui.QTableView()
        tableview.setModel(tablemodel)

        layout = QtGui.QVBoxLayout(self)
        layout.addWidget(tableview)
        self.setLayout(layout)

        self.show()

class MyTableModel(QtCore.QAbstractTableModel):
    def __init__(self, parent=None, *args):
        QtCore.QAbstractTableModel.__init__(self, parent, *args)

    def add_data(self, datain):
        self.arraydata = []
        self.arraydata.append(datain)

    def rowCount(self, parent):
        return len(self.arraydata)

    def columnCount(self, parent):
        return len(self.arraydata[0])

    def data(self, index, role):
        if not index.isValid():
            return None
        elif role != QtCore.Qt.DisplayRole:
            return None
        return (self.arraydata[index.row()][index.column()])

if __name__ == '__main__':
    app = QtCore.QCoreApplication.instance()
    if app is None:
        app = QtGui.QApplication([])

    try:
        from IPython.lib.guisupport import start_event_loop_qt4
        start_event_loop_qt4(app)
    except ImportError:
        app.exec_()
