import sys
from PyQt4 import QtGui, QtCore
from PyQt4.QtCore import *
from PyQt4.QtGui import *
import numpy as np

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
        self.tablemodel = MyTableModel(["t_rad", "Ws"])
        self.tableview = QtGui.QTableView()
        self.layout = QtGui.QVBoxLayout()

        # quit = QtGui.QPushButton('Close', self)
        # quit.setGeometry(10, 10, 60, 35)
        #
        # self.connect(quit, QtCore.SIGNAL('clicked()'),
        #              self, QtCore.SLOT('close()'))

    def show_model(self, model=None):
        """
        :param Radial1DModel: object with attribute t_rads, ws
        :return: Affective method, opens window to display table
        """
        if model is not None:
            self.tablemodel.arraydata = []
            self.add_data(model.t_rads.tolist())
            self.add_data(model.ws.tolist())
        self.tableview.setModel(self.tablemodel)
        #self.close_all_widgets(self.layout)
        self.layout.addWidget(self.tableview)
        self.setLayout(self.layout)
        self.show()

    def add_data(self, datain):
        self.tablemodel.add_data(datain)

    def close_all_widgets(self, layout):
        for i in range(layout.count()): layout.itemAt(i).widget().close()

class MyTableModel(QtCore.QAbstractTableModel):
    def __init__(self, headerdata, parent=None, *args):
        QtCore.QAbstractTableModel.__init__(self, parent, *args)
        self.arraydata = []
        self.headerdata = headerdata

    def add_data(self, datain):
        if datain is not None:
            # if type(data[0]) == np.float64:
            #     for index, item in enumerate(data):
            #         data[index] = float(item)
            #     self.arraydata.append(data)
            # elif type(data[0]) == list:
            #     for index, item in enumerate(data):
            #         self.a
            self.arraydata.append(datain)

    def rowCount(self, parent):
        return len(self.arraydata[0])

    def columnCount(self, parent):
        return len(self.arraydata)

    def headerData(self, section, orientation, role=QtCore.Qt.DisplayRole):
        if orientation == QtCore.Qt.Vertical and role == QtCore.Qt.DisplayRole:
            return QtCore.QVariant("Shell: %d" % (section + 1))
        elif orientation == QtCore.Qt.Horizontal and role == QtCore.Qt.DisplayRole:
            return QtCore.QVariant(self.headerdata[section])
        return QtCore.QVariant()

    def data(self, index, role):
        if not index.isValid():
            return None
        elif role != QtCore.Qt.DisplayRole:
            return None
        return (self.arraydata[index.column()][index.row()])
