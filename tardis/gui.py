import sys
from PyQt4 import QtGui, QtCore
from PyQt4.QtCore import *
from PyQt4.QtGui import *
import numpy

class ModelViewer(QtGui.QWidget):
    def __init__(self, parent=None):

        # assumes that qt has already been initialized in ipython with "%gui qt"
        app = QtCore.QCoreApplication.instance()
        if app is None:
            app = QtGui.QApplication([])

        try:
            from IPython.lib.guisupport import start_event_loop_qt4
            start_event_loop_qt4(app)
        except ImportError:
            app.exec_()

        QtGui.QWidget.__init__(self, parent)

        self.setGeometry(300, 300, 400, 300)
        self.setWindowTitle('Data Table')
        self.tablemodel = ModelViewer(self)

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

        self.tablemodel.add_data(model)
        self.tableview = QtGui.QTableView()
        self.tableview.setModel(self.tablemodel)

        layout = QtGui.QVBoxLayout(self)
        layout.addWidget(self.tableview)
        self.setLayout(layout)

        self.show()

class MyTableModel(QtCore.QAbstractTableModel):
    def __init__(self, parent=None, *args):
        QtCore.QAbstractTableModel.__init__(self, parent, *args)

    def add_data(self, datain):
        self.arraydata = []
        data = list(datain)
        if type((data)[0]) == numpy.float64:
            for index, item in enumerate(data):
                data[index] = float(item)
        self.arraydata.append(data)

    def rowCount(self, parent):
        return len(self.arraydata[0])

    def columnCount(self, parent):
        return len(self.arraydata)

    def data(self, index, role):
        if not index.isValid():
            return None
        elif role != QtCore.Qt.DisplayRole:
            return None
        return (self.arraydata[index.column()][index.row()])

if __name__ == '__main__':
    mdl = ModelViewer()
    mdl.show_model([numpy.float64(17126371.28374891624), numpy.float64(21276.89398667946923698464694), numpy.float64(3893164.141592654)])
