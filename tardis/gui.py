import sys, os
import numpy as np
import matplotlib
import matplotlib.pylab as plt
from matplotlib.patches import Circle, Wedge
from matplotlib.collections import PatchCollection
from matplotlib.figure import *
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from PyQt4 import QtGui, QtCore

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

        self.setGeometry(100, 100, 965, 600)
        self.setWindowTitle('Shell Viewer')
        self.tablemodel = MyTableModel(["t_rad", "Ws"])
        self.tableview = QtGui.QTableView()
        self.graph = MatplotlibWidget()
        self.layout = QtGui.QHBoxLayout()

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
        self.layout.addWidget(self.graph)
        self.layout.addWidget(self.tableview)
        self.setLayout(self.layout)
        self.plot()
        self.show()

    def update_data(self, model=None):
        if model is not None:
            self.tablemodel.arraydata = []
            self.add_data(model.t_rads.tolist())
            self.add_data(model.ws.tolist())
        for r in range(self.tablemodel.rowCount()):
            for c in range(self.tablemodel.columnCount()):
                index = self.tablemodel.createIndex(r, c)
                self.tablemodel.setData(index, QtCore.QVariant(self.tablemodel.arraydata[c][r]))

    def add_data(self, datain):
        self.tablemodel.add_data(datain)

    def plot(self):
        for i in range(20):
            globals()["shell_" + str(i)] = Wedge((0,0), i + 1, 0, 360, width=0.1)
            # globals()["shell_" + str(i)].picker = True
            self.graph.ax.add_patch(globals()["shell_" + str(i)])

        #patches = []
        #self.graph.ax.add_patch(shell_0)
        self.graph.ax.axis(xmin=-20,xmax=20,ymin=-20,ymax=20)
        self.graph.draw()

class MyTableModel(QtCore.QAbstractTableModel):
    def __init__(self, headerdata, parent=None, *args):
        QtCore.QAbstractTableModel.__init__(self, parent, *args)
        self.arraydata = []
        self.headerdata = headerdata

    def add_data(self, datain):
        self.arraydata.append(datain)

    def rowCount(self, parent=QtCore.QModelIndex()):
        return len(self.arraydata[0])

    def columnCount(self, parent=QtCore.QModelIndex()):
        return len(self.arraydata)

    def headerData(self, section, orientation, role=QtCore.Qt.DisplayRole):
        if orientation == QtCore.Qt.Vertical and role == QtCore.Qt.DisplayRole:
            return QtCore.QVariant("Shell: %d" % (section + 1))
        elif orientation == QtCore.Qt.Horizontal and role == QtCore.Qt.DisplayRole:
            return QtCore.QVariant(self.headerdata[section])
        return QtCore.QVariant()

    def data(self, index, role=QtCore.Qt.DisplayRole):
        if not index.isValid():
            return None
        elif role != QtCore.Qt.DisplayRole:
            return None
        return (self.arraydata[index.column()][index.row()])

    def setData(self, index, value, role=QtCore.Qt.EditRole):
        if not index.isValid():
            return False
        elif role != QtCore.Qt.EditRole:
            return False
        self.arraydata[index.column()][index.row()] = value
        self.emit(SIGNAL('dataChanged(const QModelIndex &, const QModelIndex &)'), index, index)
        return True

class MatplotlibWidget(FigureCanvas):

    def __init__(self):
        self.figure = Figure()
        self.ax = self.figure.add_subplot(111)

        FigureCanvas.__init__(self, self.figure)
        FigureCanvas.setSizePolicy(self, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Expanding)
        FigureCanvas.updateGeometry(self)

    #     cid = self.figure.canvas.mpl_connect('button_press_event', self.onclick)
    #
    # def onclick(event):
    #     print 'button=%d, x=%d, y=%d, xdata=%f, ydata=%f'% (
    #     event.button, event.x, event.y, event.xdata, event.ydata)
    #
    # def pick_handler(event):
    #     mouseevent = event.mouseevent
    #     artist = event.artist
    #     print event
