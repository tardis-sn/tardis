import numpy as np
import matplotlib
import matplotlib.pylab as plt
from matplotlib.patches import Circle, Wedge
from matplotlib import colors
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
        self.graph = MatplotlibWidget(self)
        self.layout = QtGui.QHBoxLayout()
        self.model = None

    def show_model(self, model=None):
        """
        :param Radial1DModel: object with attribute t_rads, ws
        :return: Affective method, opens window to display table
        """
        if model:
            self.change_model(model)
        self.tableview.setModel(self.tablemodel)
        self.layout.addWidget(self.graph)
        self.layout.addWidget(self.tableview)
        self.setLayout(self.layout)
        self.plot()
        self.show()

    def update_data(self, model=None):
        if model:
            self.change_model(model)
        for r in range(self.tablemodel.rowCount()):
            for c in range(self.tablemodel.columnCount()):
                index = self.tablemodel.createIndex(r, c)
                self.tablemodel.setData(index, QtCore.QVariant(self.tablemodel.arraydata[c][r]))
        self.plot()

    def change_model(self, model):
        self.model = model
        self.tablemodel.arraydata = []
        self.add_data(model.t_rads.tolist())
        self.add_data(model.ws.tolist())

    def add_data(self, datain):
        self.tablemodel.add_data(datain)

    def plot(self):
        self.graph.ax.clear()
        # set title and axis labels
        self.shells = []
        t_rad_normalizer = colors.Normalize(vmin=self.model.t_rads.min(), vmax=self.model.t_rads.max())
        t_rad_color_map = plt.cm.ScalarMappable(norm=t_rad_normalizer, cmap=plt.cm.jet)
        t_rad_color_map.set_array(self.model.t_rads)
        if self.graph.cb:
            self.graph.cb.set_clim(vmin=self.model.t_rads.min(), vmax=self.model.t_rads.max())
            self.graph.cb.update_normal(t_rad_color_map)
            #self.graph.figure.delaxes(self.graph.figure.axes[1])
            #self.graph.figure.subplots_adjust(right=0.90)  #default right padding
        else:
            self.graph.cb = self.graph.figure.colorbar(t_rad_color_map)
            self.graph.cb.set_label('T (K)')
        #normalizing_factor = 0.2 * (self.model.r_outer[-1] - self.model.r_inner[0]) / self.model.r_inner[0]
        normalizing_factor = 0.8e-15
        print normalizing_factor
        for i, t_rad in enumerate(self.model.t_rads):
            r_inner = self.model.r_inner[i] * normalizing_factor
            r_outer = self.model.r_outer[i] * normalizing_factor
            self.shells.append(Shell(i, (0,0), r_inner, r_outer, facecolor=t_rad_color_map.to_rgba(t_rad),
                                picker=self.graph.shell_picker))
            # shells[i].set_out_radius(self.model.r_outer[i] / 1e15)
            # shells[i].set_in_radius(self.model.r_inner[i] / 1e15)
            self.graph.ax.add_patch(self.shells[i])
        self.graph.ax.set_xlim(0, 2)
        self.graph.ax.set_ylim(0, 2)
        #self.graph.ax.axis(xmin=0,xmax=20,ymin=0,ymax=20)
        print self.model.r_outer[-1] / 1e15
        #self.graph.ax.axis(xmin=0,xmax=self.model.r_outer[-1] / 1e15,ymin=0,ymax=self.model.r_outer[-1] / 1e15)
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
        self.emit(QtCore.SIGNAL('dataChanged(const QModelIndex &, const QModelIndex &)'), index, index)
        return True

class MatplotlibWidget(FigureCanvas):

    def __init__(self, parent):
        self.parent = parent
        self.figure = Figure()
        self.ax = self.figure.add_subplot(111)
        self.cb = None

        FigureCanvas.__init__(self, self.figure)
        FigureCanvas.setSizePolicy(self, QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Expanding)
        FigureCanvas.updateGeometry(self)

        cid = self.figure.canvas.mpl_connect('pick_event', self.onpick)

    def onpick(self, event):
        mouseevent = event.mouseevent
        shell = event.artist
        print 'Picked event:', event, mouseevent.xdata, mouseevent.ydata, shell.get_out_radius()
        self.parent.tableview.selectRow(shell.get_index())

    def shell_picker(self, shell, mouseevent):
        """
        find the points within a certain distance from the mouseclick in
        data coords and attach some extra attributes, pickx and picky
        which are the data points that were picked
        """
        tolerance = 0.001 # distance that click is allowed from shell, based on matplotlib data coordinates
        if mouseevent.xdata is None:
            return False, dict()
        mouse_r2 = mouseevent.xdata ** 2 + mouseevent.ydata ** 2
        if (shell.r_inner - tolerance) ** 2 < mouse_r2 < (shell.r_outer + tolerance) ** 2:
            return True, dict()
        return False, dict()

class Shell(matplotlib.patches.Wedge):

    def __init__(self, index, center, r_inner, r_outer, **kwargs):
        super(Shell, self).__init__(center, r_outer, 0, 90, width=r_outer - r_inner, **kwargs)
        self.index = index
        self.center = center
        self.r_outer = r_outer
        self.r_inner = r_inner
        self.width = r_outer - r_inner

    def get_index(self):
        return self.index

    def get_out_radius(self):
        return self.r_outer

    def get_in_radius(self):
        return self.radius - self.width

    def get_center(self):
        return self.center

    def get_width(self):
        return self.width

    def set_out_radius(self, r):
        if r > self.get_in_radius():
            self.radius = r
            return True
        return False

    def set_in_radius(self, r):
        if 0 < r < self.get_out_radius():
            self.width = self.get_out_radius() - r
            return True
        return False

    def set_center(self, center):
        self.center = center
        return True

    def set_width(self, width):
        if width > 0:
            self.width = width
            return True
        return False
