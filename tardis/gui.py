import numpy as np
import matplotlib
#matplotlib.use('KtAgg')
import matplotlib.pylab as plt
import matplotlib.gridspec as gridspec
from matplotlib import colors
from matplotlib.patches import Circle, Wedge
from matplotlib.figure import *
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt4 import NavigationToolbar2QT as NavigationToolbar
from PyQt4 import QtGui, QtCore

class ModelViewer(QtGui.QWidget):
    def __init__(self, parent=None):
        # assumes that qt has already been initialized by starting IPython with the flag "--pylab=qt"
        app = QtCore.QCoreApplication.instance()
        if app is None:
            app = QtGui.QApplication([])
        try:
            from IPython.lib.guisupport import start_event_loop_qt4
            start_event_loop_qt4(app)
        except ImportError:
            app.exec_()

        super(ModelViewer, self).__init__(parent)
        self.model = None
        self.setGeometry(20, 35, 1200, 500)
        self.setWindowTitle('Shells Viewer')
        self.tablemodel = MyTableModel(["t_rad", "Ws"])
        self.tableview = QtGui.QTableView()
        self.graph = MatplotlibWidget(self, 'model')
        self.graph_label = QtGui.QLabel('Select Property to Plot:')
        self.graph_button = QtGui.QToolButton()
        self.spectrum = MatplotlibWidget(self)
        self.spectrum_label = QtGui.QLabel('Select Spectrum to Plot:')
        self.spectrum_button = QtGui.QToolButton()
        self.layout = QtGui.QHBoxLayout()
        self.graph_sublayout = QtGui.QVBoxLayout()
        self.graph_subsublayout = QtGui.QHBoxLayout()
        self.spectrum_sublayout = QtGui.QVBoxLayout()
        self.spectrum_subsublayout = QtGui.QHBoxLayout()

    def show_model(self, model=None):
        if model:
            self.change_model(model)
        self.tableview.setModel(self.tablemodel)
        self.tableview.setMinimumWidth(200)
        self.tableview.connect(self.tableview.verticalHeader(), QtCore.SIGNAL('sectionDoubleClicked(int)'),
                               self.on_header_double_clicked)
        self.graph_button.setText('t_rads')
        self.spectrum_button.setText('spec_flux_angstrom')
        self.graph_button.setPopupMode(QtGui.QToolButton.MenuButtonPopup)
        self.spectrum_button.setPopupMode(QtGui.QToolButton.MenuButtonPopup)
        self.graph_button.setMenu(QtGui.QMenu(self.graph_button))
        self.spectrum_button.setMenu(QtGui.QMenu(self.spectrum_button))
        self.graph_button.menu().addAction('t_rads').triggered.connect(self.change_graph_to_t_rads)
        self.graph_button.menu().addAction('Ws').triggered.connect(self.change_graph_to_ws)
        self.spectrum_button.menu().addAction('spec_flux_angstrom').triggered.connect(self.change_spectrum_to_spec_flux_angstrom)
        self.spectrum_button.menu().addAction('spec_virtual_flux_angstrom').triggered.connect(self.change_spectrum_to_spec_virtual_flux_angstrom)
        self.layout.addWidget(self.tableview)
        self.graph_subsublayout.addWidget(self.graph_label)
        self.graph_subsublayout.addWidget(self.graph_button)
        self.graph_sublayout.addLayout(self.graph_subsublayout)
        self.graph_sublayout.addWidget(self.graph)
        self.layout.addLayout(self.graph_sublayout)
        self.spectrum_subsublayout.addWidget(self.spectrum_label)
        self.spectrum_subsublayout.addWidget(self.spectrum_button)
        self.spectrum_sublayout.addLayout(self.spectrum_subsublayout)
        self.spectrum_sublayout.addWidget(self.spectrum)
        self.spectrum_sublayout.addWidget(self.spectrum.toolbar)
        self.layout.addLayout(self.spectrum_sublayout)
        self.setLayout(self.layout)
        self.plot_model()
        self.plot_spectrum()
        self.show()

    def update_data(self, model=None):
        if model:
            self.change_model(model)
        for r in range(self.tablemodel.rowCount()):
            for c in range(self.tablemodel.columnCount()):
                index = self.tablemodel.createIndex(r, c)
                self.tablemodel.setData(index, QtCore.QVariant(self.tablemodel.arraydata[c][r]))
        self.plot_model()
        if self.graph_button.text == 'Ws':
            self.change_graph_to_ws()
        self.plot_spectrum()
        if self.spectrum_button.text == 'spec_virtual_flux_angstrom':
            self.change_spectrum_to_spec_virtual_flux_angstrom()
        self.show()

    def change_model(self, model):
        self.model = model
        self.tablemodel.arraydata = []
        self.add_data(model.t_rads.tolist())
        self.add_data(model.ws.tolist())

    def add_data(self, datain):
        self.tablemodel.add_data(datain)

    def change_spectrum_to_spec_virtual_flux_angstrom(self):
        self.change_spectrum(self.model.spec_virtual_flux_angstrom, 'spec_virtual_flux_angstrom')

    def change_spectrum_to_spec_flux_angstrom(self):
        self.change_spectrum(self.model.spec_flux_angstrom, 'spec_flux_angstrom')

    def change_spectrum(self, data, name):
        self.spectrum_button.setText(name)
        self.spectrum.dataplot[0].set_ydata(data)
        self.spectrum.ax.relim()
        self.spectrum.ax.autoscale()
        self.spectrum.draw()

    def plot_spectrum(self):
        self.spectrum.ax.clear()
        self.spectrum.ax.set_title('Spectrum')
        self.spectrum.ax.set_xlabel('Wavelength (A)')
        self.spectrum.ax.set_ylabel('Intensity')
        self.spectrum.dataplot = self.spectrum.ax.plot(self.model.spec_angstrom, self.model.spec_flux_angstrom, label='b')
        self.spectrum.draw()

    def change_graph_to_ws(self):
        self.change_graph(self.model.ws, 'Ws', '')

    def change_graph_to_t_rads(self):
        self.change_graph(self.model.t_rads, 't_rads', '(K)')

    def change_graph(self, data, name, unit):
        self.graph_button.setText(name)
        self.graph.dataplot[0].set_ydata(data)
        self.graph.ax1.relim()
        self.graph.ax1.autoscale()
        self.graph.ax1.set_title(name + ' vs Shell')
        self.graph.ax1.set_ylabel(name + ' ' + unit)
        normalizer = colors.Normalize(vmin=data.min(), vmax=data.max())
        color_map = plt.cm.ScalarMappable(norm=normalizer, cmap=plt.cm.jet)
        color_map.set_array(data)
        self.graph.cb.set_clim(vmin=data.min(), vmax=data.max())
        self.graph.cb.update_normal(color_map)
        if unit == '(K)':
            unit = 'T (K)'
        self.graph.cb.set_label(unit)
        for i, item in enumerate(data):
            self.shells[i].set_facecolor(color_map.to_rgba(item))
        self.graph.draw()

    def plot_model(self):
        self.graph.ax1.clear()
        self.graph.ax1.set_title('t_rads vs Shell')
        self.graph.ax1.set_xlabel('Shell Number')
        self.graph.ax1.set_ylabel('t_rads (K)')
        self.graph.ax1.yaxis.get_major_formatter().set_powerlimits((0, 1))
        self.graph.dataplot = self.graph.ax1.plot(range(len(self.model.t_rads)), self.model.t_rads)
        self.graph.ax2.clear()
        self.graph.ax2.set_title('Shell View')
        self.graph.ax2.set_xlabel('Distance from Center (cm)')
        self.graph.ax2.set_ylabel('Distance from Center (cm)')
        self.shells = []
        t_rad_normalizer = colors.Normalize(vmin=self.model.t_rads.min(), vmax=self.model.t_rads.max())
        t_rad_color_map = plt.cm.ScalarMappable(norm=t_rad_normalizer, cmap=plt.cm.jet)
        t_rad_color_map.set_array(self.model.t_rads)
        if self.graph.cb:
            self.graph.cb.set_clim(vmin=self.model.t_rads.min(), vmax=self.model.t_rads.max())
            self.graph.cb.update_normal(t_rad_color_map)
        else:
            self.graph.cb = self.graph.figure.colorbar(t_rad_color_map)
            self.graph.cb.set_label('T (K)')
        self.graph.normalizing_factor = 0.2 * (self.model.r_outer[-1] - self.model.r_inner[0]) / self.model.r_inner[0]
        #self.graph.normalizing_factor = 8e-16
        for i, t_rad in enumerate(self.model.t_rads):
            r_inner = self.model.r_inner[i] * self.graph.normalizing_factor
            r_outer = self.model.r_outer[i] * self.graph.normalizing_factor
            self.shells.append(Shell(i, (0,0), r_inner, r_outer, facecolor=t_rad_color_map.to_rgba(t_rad),
                                     picker=self.graph.shell_picker))
            self.graph.ax2.add_patch(self.shells[i])
        self.graph.ax2.set_xlim(0, self.model.r_outer[-1] * self.graph.normalizing_factor)
        self.graph.ax2.set_ylim(0, self.model.r_outer[-1] * self.graph.normalizing_factor)
        #self.graph.gs.tight_layout(self.graph.figure)
        self.graph.draw()

    def on_header_double_clicked(self, index):
        shell_info = ShellInfo(index, self)

class ShellInfo(QtGui.QDialog):

    def __init__(self, index, parent=None):
        super(ShellInfo, self).__init__(parent)
        self.index = index
        self.setGeometry(500, 150, 650, 650)
        self.setWindowTitle('Shell %d Info' % (self.index + 1))
        #self.graph = MatplotlibWidget(self)
        self.atoms = MyTableModel(['Atoms'])
        self.layout = QtGui.QVBoxLayout()
        self.layout.addWidget(self.graph)
        self.layout.addWidget(self.graph.toolbar)
        self.setLayout(self.layout)
        self.show()

class MyTableModel(QtCore.QAbstractTableModel):
    def __init__(self, headerdata, parent=None, *args):
        super(MyTableModel, self).__init__(parent, *args)
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

    def __init__(self, parent, fig=None):
        self.parent = parent
        self.figure = Figure()
        if fig != 'model':
            self.ax = self.figure.add_subplot(111)
        else:
            self.gs = gridspec.GridSpec(2, 1, height_ratios=[1, 3])
            self.ax1 = self.figure.add_subplot(self.gs[0])
            self.ax2 = self.figure.add_subplot(self.gs[1])
        self.cb = None

        super(MatplotlibWidget, self).__init__(self.figure)
        super(MatplotlibWidget, self).setSizePolicy(QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Expanding)
        super(MatplotlibWidget, self).updateGeometry()
        if fig != 'model':
            self.toolbar = NavigationToolbar(self, parent)
        else:
            cid = self.figure.canvas.mpl_connect('pick_event', self.onpick)

    def onpick(self, event):
        self.parent.tableview.selectRow(event.artist.index)
        for i in range(len(self.parent.shells)):
            if i != event.artist.index and i != event.artist.index + 1:
                self.parent.shells[i].set_edgecolor('k')
            else:
                self.parent.shells[i].set_edgecolor('w')
        self.draw()

    def shell_picker(self, shell, mouseevent):
        if mouseevent.xdata is None:
            return False, dict()
        mouse_r2 = mouseevent.xdata ** 2 + mouseevent.ydata ** 2
        if shell.r_inner ** 2 < mouse_r2 < shell.r_outer ** 2:
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
