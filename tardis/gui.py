import numpy as np
import pdb
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
from astropy import units as u
import analysis

# def current_ion_index(index, index_list):
#     if not index in index_list:
#         return None
#     if not (index - 1) in index_list:
#         return 0
#     else:
#         return current_ion_index(index - 1, index_list) + 1
#
# def current_ion_index(index, duplicate_list):
#     if duplicate_list[index - 1] != duplicate_list[index]:
#         return 0
#     else:
#         return current_ion_index(index - 1, duplicate_list) + 1

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
        self.shell_info = {}
        self.line_info = []
        self.setGeometry(20, 35, 1250, 500)
        self.setWindowTitle('Shells Viewer')
        self.tablemodel = MyTableModel([['Shell: '], ["t_rad", "Ws"]], (1, 0))
        self.tableview = QtGui.QTableView()
        self.graph = MatplotlibWidget(self, 'model')
        self.graph_label = QtGui.QLabel('Select Property:')
        self.graph_button = QtGui.QToolButton()
        self.spectrum = MatplotlibWidget(self)
        self.spectrum_label = QtGui.QLabel('Select Spectrum:')
        self.spectrum_button = QtGui.QToolButton()
        self.spectrum_span_button = QtGui.QPushButton('Show Wavelength Range')
        self.layout = QtGui.QHBoxLayout()
        self.graph_sublayout = QtGui.QVBoxLayout()
        self.graph_subsublayout = QtGui.QHBoxLayout()
        self.spectrum_sublayout = QtGui.QVBoxLayout()
        self.spectrum_subsublayout = QtGui.QHBoxLayout()

        self.tableview.setMinimumWidth(200)
        self.tableview.connect(self.tableview.verticalHeader(), QtCore.SIGNAL('sectionClicked(int)'), self.graph.highlight_shell)
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
        self.spectrum_span_button.clicked.connect(self.spectrum.show_span)
        self.layout.addWidget(self.tableview)
        self.graph_subsublayout.addWidget(self.graph_label)
        self.graph_subsublayout.addWidget(self.graph_button)
        self.graph_sublayout.addLayout(self.graph_subsublayout)
        self.graph_sublayout.addWidget(self.graph)
        self.layout.addLayout(self.graph_sublayout)
        self.spectrum_subsublayout.addWidget(self.spectrum_span_button)
        self.spectrum_subsublayout.addWidget(self.spectrum_label)
        self.spectrum_subsublayout.addWidget(self.spectrum_button)
        self.spectrum_sublayout.addLayout(self.spectrum_subsublayout)
        self.spectrum_sublayout.addWidget(self.spectrum)
        self.spectrum_sublayout.addWidget(self.spectrum.toolbar)
        self.layout.addLayout(self.spectrum_sublayout)
        self.setLayout(self.layout)

    def show_model(self, model=None):
        if model:
            self.change_model(model)
        self.tableview.setModel(self.tablemodel)
        self.plot_model()
        self.plot_spectrum()
        self.show()

    def update_data(self, model=None):
        if model:
            self.change_model(model)
        self.tablemodel.updateTable()
        for index in self.shell_info.keys():
            self.shell_info[index].update_tables()
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
        self.tablemodel.addData(model.t_rads.tolist())
        self.tablemodel.addData(model.ws.tolist())

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
        self.shell_info[index] = ShellInfo(index, self)

class ShellInfo(QtGui.QDialog):

    def __init__(self, index, parent=None):
        super(ShellInfo, self).__init__(parent)
        self.parent = parent
        self.shell_index = index
        self.setGeometry(400, 150, 200, 400)
        self.setWindowTitle('Shell %d Abundances' % (self.shell_index + 1))
        self.atomstable = QtGui.QTableView()
        self.ionstable = QtGui.QTableView()
        self.levelstable = QtGui.QTableView()
        self.atomstable.connect(self.atomstable.verticalHeader(), QtCore.SIGNAL('sectionDoubleClicked(int)'),
                               self.on_atom_header_double_clicked)

        self.plasma = self.parent.model.plasmas[self.shell_index]
        self.table1_data = self.plasma.number_density
        self.atomsdata = MyTableModel([['Z = '], ['Count (Shell %d)' % (self.shell_index + 1)]], iterate_header=(2, 0), index_info=self.table1_data.index.values.tolist())
        self.ionsdata = None
        self.levelsdata = None
        self.atomsdata.addData(self.table1_data.values.tolist())
        self.atomstable.setModel(self.atomsdata)

        self.layout = QtGui.QHBoxLayout()
        self.layout.addWidget(self.atomstable)
        self.layout.addWidget(self.ionstable)
        self.layout.addWidget(self.levelstable)
        self.setLayout(self.layout)
        self.ionstable.hide()
        self.levelstable.hide()
        self.show()

    def on_atom_header_double_clicked(self, index):
        self.current_atom_index = self.table1_data.index.values.tolist()[index]
        self.table2_data = self.plasma.ion_populations.ix[self.current_atom_index]
        self.ionsdata = MyTableModel([['Ion: '], ['Count (Z = %d)' % self.current_atom_index]], iterate_header=(2, 0), index_info=self.table2_data.index.values.tolist())
        normalized_data = []
        for item in self.table2_data.values.tolist():
            normalized_data.append(float(item / self.table1_data.ix[self.current_atom_index]))
        self.ionsdata.addData(normalized_data)
        self.ionstable.setModel(self.ionsdata)
        self.ionstable.connect(self.ionstable.verticalHeader(), QtCore.SIGNAL('sectionDoubleClicked(int)'),
                               self.on_ion_header_double_clicked)
        self.levelstable.hide()
        self.ionstable.setColumnWidth(0, 120)
        self.ionstable.show()
        self.setGeometry(400, 150, 380, 400)
        self.show()

    def on_ion_header_double_clicked(self, index):
        self.current_ion_index = self.table2_data.index.values.tolist()[index]
        self.table3_data = self.plasma.level_populations.ix[self.current_atom_index, self.current_ion_index]
        self.levelsdata = MyTableModel([['Level: '], ['Count (Ion %d)' % self.current_ion_index]], iterate_header=(2, 0), index_info=self.table3_data.index.values.tolist())
        normalized_data = []
        for item in self.table3_data.values.tolist():
            normalized_data.append(float(item / self.table2_data.ix[self.current_ion_index]))
        self.levelsdata.addData(normalized_data)
        self.levelstable.setModel(self.levelsdata)
        self.levelstable.setColumnWidth(0, 120)
        self.levelstable.show()
        self.setGeometry(400, 150, 600, 400)
        self.show()

    def update_tables(self):
        self.plasma = self.parent.model.plasmas[self.shell_index]
        self.table1_data = self.plasma.number_density
        self.atomsdata.index_info=self.table1_data.index.values.tolist()
        self.atomsdata.arraydata = []
        self.atomsdata.addData(self.table1_data.values.tolist())
        self.atomsdata.updateTable()
        self.ionstable.hide()
        self.levelstable.hide()
        self.setGeometry(400, 150, 200, 400)
        self.show()

class LineInfo(QtGui.QDialog):

    def __init__(self, parent, wavelength_start, wavelength_end):
        super(LineInfo, self).__init__(parent)
        self.parent = parent
        self.setGeometry(200 + len(self.parent.line_info) * 20, 150, 250, 400)
        #self.setWindowTitle('Last Line Interaction: %f - %f (A)' % (self.wavelength_start, self.wavelength_end))
        self.setWindowTitle('Line Interaction')
        self.atomstable = QtGui.QTableView()
        self.transitionsintable = QtGui.QTableView()
        self.transitionsouttable = QtGui.QTableView()
        #self.atomstable.connect(self.atomstable.verticalHeader(), QtCore.SIGNAL('sectionClicked(int)'), self.on_atom_header_clicked)
        self.atomstable.connect(self.atomstable, QtCore.SIGNAL("clicked(QModelIndex)"), self.on_atom_clicked)
        self.get_data(wavelength_start, wavelength_end)
        #self.last_line_out.groupby('atomic_number').wavelength.count().astype(float) / self.last_line_out.groupby('atomic_number').wavelength.count().sum()
        self.atomsdata = MyTableModel([self.header_list, ['Abundance']])
        #iterate_header=(2, 0), index_info=self.atom_values.index.values.tolist())
        self.atomsdata.addData(self.ion_table)
        #self.transitionsdata = MyTableModel([])
        #self.atomsdata.addData(self.level_transitions_in)
        #self.atomsdata.addData(self.last_line_out_atom_table.values.tolist())
        self.atomstable.setModel(self.atomsdata)
        # for ions in self.ions_index:
        #     self.atomstable.hideRow(ions)
        self.layout = QtGui.QHBoxLayout()
        self.layout.addWidget(self.atomstable)
        self.layout.addWidget(self.transitionsintable)
        self.layout.addWidget(self.transitionsouttable)
        self.setLayout(self.layout)
        self.transitionsintable.hide()
        self.transitionsouttable.hide()
        self.show()

    def get_data(self, wavelength_start, wavelength_end):
        self.wavelength_start = wavelength_start * u.angstrom
        self.wavelength_end = wavelength_end * u.angstrom
        last_line_in_ids, last_line_out_ids = analysis.get_last_line_interaction(self.wavelength_start, self.wavelength_end, self.parent.model)
        self.last_line_in, self.last_line_out = self.parent.model.atom_data.lines.ix[last_line_in_ids], self.parent.model.atom_data.lines.ix[last_line_out_ids]
        self.grouped_lines_in, self.grouped_lines_out = self.last_line_in.groupby(['atomic_number', 'ion_number']), self.last_line_out.groupby(['atomic_number', 'ion_number'])
        #self.atom_values = self.last_line_in.groupby('atomic_number').wavelength.count().astype(float) / self.last_line_in.groupby('atomic_number').wavelength.count().sum()
        self.ions_in, self.ions_out = self.grouped_lines_in.groups.keys(), self.grouped_lines_out.groups.keys()
        self.ions_in.sort()
        self.ions_out.sort()
        #self.current_atom_index = []
        #self.ions_index = []
        self.header_list = []
        self.ion_table = (self.grouped_lines_in.wavelength.count().astype(float) / self.grouped_lines_in.wavelength.count().sum()).values.tolist()
        for z, ion in self.ions_in:
            self.header_list.append('Z = %d: Ion %d' % (z, ion))
        # for index, item in enumerate(self.atom_values.values.tolist()):
        #     #self.full_table_data.append(item)
        #     current_atom_index = self.atom_values.index.values.tolist()[index]
        #     for i_index, item in enumerate((self.last_line_in[self.last_line_in.atomic_number == current_atom_index].groupby('ion_number').wavelength.count().astype(float) / self.last_line_in.groupby('ion_number').wavelength.count().sum()).values.tolist()):
        #         self.full_table_data.append(item)
        #         self.current_atom_index.append(current_atom_index)
        #         self.header_list.append('Z = ' + str(current_atom_index) + ': Ion ' + str(i_index))
        #         #self.ions_index.append(len(self.full_table_data) - 1)

    def update(self, wavelength_start, wavelength_end):
        self.get_data(wavelength_start, wavelength_end)
        #self.atomsdata.headerdata = [self.header_list, ['Percent']]
        #self.emit(QtCore.SIGNAL("LayoutAboutToBeChanged()"))
        self.atomsdata = MyTableModel([self.header_list, ['Percent']])
        self.atomsdata.arraydata = []
        self.atomsdata.addData(self.full_table_data)
        self.atomsdata.updateTable()
        #self.emit(QtCore.SIGNAL("LayoutChanged()"))
        self.atomstable = QtGui.QTableView()
        self.atomstable.connect(self.atomstable.verticalHeader(), QtCore.SIGNAL('sectionClicked(int)'), self.on_atom_header_clicked)
        self.atomstable.setModel(self.atomsdata)
        for row in range(self.atomsdata.rowCount()):
            if not (row in self.ions_index):
                self.atomstable.showRow(row)
            else:
                self.atomstable.hideRow(row)
        self.show()

    def get_transition_table(self, lines, atom, ion):
        grouped = lines.groupby(['atomic_number', 'ion_number'])
        transitions_with_duplicates = lines.ix[grouped.groups[(atom, ion)]].groupby(['level_number_lower', 'level_number_upper']).groups
        transitions = lines.ix[grouped.groups[(atom, ion)]].drop_duplicates().groupby(['level_number_lower', 'level_number_upper']).groups
        transitions_count = []
        transitions_parsed = []
        for item in transitions.values():
            c = 0
            for ditem in transitions_with_duplicates.values():
                c += ditem.count(item[0])
            transitions_count.append(c)
        s = 0
        for item in transitions_count:
            s += item
        for index in range(len(transitions_count)):
            transitions_count[index] /= float(s)
        for key, value in transitions.items():
            #pdb.set_trace()
            transitions_parsed.append("%d-%d (%.2f A)" % (key[0], key[1], self.parent.model.atom_data.lines.ix[value[0]]['wavelength']))
        return transitions_parsed, transitions_count

    def on_atom_clicked(self, qindex):
        index = qindex.row()
        # gindex = index
        # if gindex in self.ions_index:
        #     hack = 1
        # else:
        #     hack = 0
        # while not gindex in self.ions_index:
        #     gindex += 1     # Last item in table is always an ion
        #current_atom_index = self.atom_values.index.values.tolist()[index - self.ions_index.index(gindex) - hack]
        # if not index in self.ions_index:
        #     self.transitions_level_lower = self.last_line_in[(self.last_line_in.atomic_number == current_atom_index)].level_number_lower.values.tolist()
        #     self.transitions_level_upper = self.last_line_in[(self.last_line_in.atomic_number == current_atom_index)].level_number_upper.values.tolist()
        #     self.transitions_wavelength = self.last_line_in[(self.last_line_in.atomic_number == current_atom_index)].wavelength.values.tolist()
        # else:
        #i_index = current_ion_index(index, self.ions_index)
        # current_atom_index = self.[index]
        # i_index = current_ion_index(index, self.current_atom_index)
        # self.transitions_level_lower = self.last_line_in[(self.last_line_in.atomic_number == current_atom_index) & (self.last_line_in.ion_number == i_index)].level_number_lower.values.tolist()
        # self.transitions_level_upper = self.last_line_in[(self.last_line_in.atomic_number == current_atom_index) & (self.last_line_in.ion_number == i_index)].level_number_upper.values.tolist()
        # self.transitions_wavelength = self.last_line_in[(self.last_line_in.atomic_number == current_atom_index) & (self.last_line_in.ion_number == i_index)].wavelength.values.tolist()
        self.transitionsin_parsed, self.transitionsin_count = self.get_transition_table(self.last_line_in, self.ions_in[index][0], self.ions_in[index][1])
        self.transitionsout_parsed, self.transitionsout_count = self.get_transition_table(self.last_line_out, self.ions_out[index][0], self.ions_out[index][1])
        self.transitionsindata = MyTableModel([self.transitionsin_parsed, ['Lines In']])
        self.transitionsoutdata = MyTableModel([self.transitionsout_parsed, ['Lines Out']])
        self.transitionsindata.addData(self.transitionsin_count)
        self.transitionsoutdata.addData(self.transitionsout_count)
        self.transitionsintable.setModel(self.transitionsindata)
        self.transitionsouttable.setModel(self.transitionsoutdata)
        self.transitionsintable.show()
        self.transitionsouttable.show()
        self.setGeometry(400, 150, 750, 400)
        self.show()

    def on_atom_header_clicked(self, index):
        #self.current_atom_index = self.atom_values.index.values.tolist()[index]
        #self.last_line_out[self.last_line_out.atomic_number == self.current_atom_index].groupby('ion_number').wavelength.count().astype(float) / self.last_line_out[self.last_line_out.atomic_number == self.current_atom_index].groupby('ion_number').wavelength.count().sum()
        #self.ionsdata = MyTableModel([['Ion: '], ['Lines In']], iterate_header=(2, 0), index_info=self.last_line_in_ion_table.index.values.tolist())
        #self.ionsdata.addData(self.last_line_in_ion_table.values.tolist())
        #self.ionsdata.addData(self.last_line_out_ion_table.values.tolist())
        #self.ionstable.setModel(self.ionsdata)
        #self.ionstable.setColumnWidth(0, 120)
        #self.ionstable.show()
        #self.setGeometry(400, 150, 400, 400)
        i_index = index + 1
        if not index in self.ions_index:
            while i_index in self.ions_index:
                if self.atomstable.isRowHidden(i_index):
                    self.atomstable.showRow(i_index)
                else:
                    self.atomstable.hideRow(i_index)
                i_index += 1
        self.show()

class MyTableModel(QtCore.QAbstractTableModel):
    def __init__(self, headerdata=None, iterate_header=(0, 0), index_info=None, parent=None, *args):
        super(MyTableModel, self).__init__(parent, *args)
        self.headerdata = headerdata
        self.arraydata = []
        self.iterate_header = iterate_header
        self.index_info = index_info

    def addData(self, datain):
        self.arraydata.append(datain)

    def rowCount(self, parent=QtCore.QModelIndex()):
        return len(self.arraydata[0])

    def columnCount(self, parent=QtCore.QModelIndex()):
        return len(self.arraydata)

    def headerData(self, section, orientation, role=QtCore.Qt.DisplayRole):
        if orientation == QtCore.Qt.Vertical and role == QtCore.Qt.DisplayRole:
            if self.iterate_header[0] == 1:
                return QtCore.QVariant(self.headerdata[0][0] + str(section + 1))
            elif self.iterate_header[0] == 2:
                if self.index_info:
                    return QtCore.QVariant(self.headerdata[0][0] + str(self.index_info[section]))
                else:
                    return QtCore.QVariant(self.headerdata[0][0] + str(section + 1))
            else:
                return QtCore.QVariant(self.headerdata[0][section])
        elif orientation == QtCore.Qt.Horizontal and role == QtCore.Qt.DisplayRole:
            if self.iterate_header[1] == 1:
                return QtCore.QVariant(self.headerdata[1][0] + str(section + 1))
            elif self.iterate_header[1] == 2:
                if self.index_info:
                    return QtCore.QVariant(self.headerdata[1][0] + str(self.index_info[section]))
            else:
                return QtCore.QVariant(self.headerdata[1][section])
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

    def updateTable(self):
        for r in range(self.rowCount()):
            for c in range(self.columnCount()):
                index = self.createIndex(r, c)
                self.setData(index, QtCore.QVariant(self.arraydata[c][r]))

class MatplotlibWidget(FigureCanvas):

    def __init__(self, parent, fig=None):
        self.parent = parent
        self.figure = Figure()
        self.cid = {}
        if fig != 'model':
            self.ax = self.figure.add_subplot(111)
        else:
            self.gs = gridspec.GridSpec(2, 1, height_ratios=[1, 3])
            self.ax1 = self.figure.add_subplot(self.gs[0])
            self.ax2 = self.figure.add_subplot(self.gs[1])
        self.cb = None
        self.span = None

        super(MatplotlibWidget, self).__init__(self.figure)
        super(MatplotlibWidget, self).setSizePolicy(QtGui.QSizePolicy.Expanding, QtGui.QSizePolicy.Expanding)
        super(MatplotlibWidget, self).updateGeometry()
        if fig != 'model':
            self.toolbar = NavigationToolbar(self, parent)
            self.cid[0] = self.figure.canvas.mpl_connect('pick_event', self.on_span_pick)
        else:
            self.cid[0] = self.figure.canvas.mpl_connect('pick_event', self.on_shell_pick)

    def show_span(self, garbage=0, left=5000, right=10000):
        if self.parent.spectrum_span_button.text() == 'Show Wavelength Range':
            if not self.span:
                self.span = self.ax.axvspan(left, right, color='r', alpha=0.3, picker=self.span_picker)
            else:
                self.span.set_visible(True)
            self.parent.spectrum_span_button.setText('Hide Wavelength Range')
        else:
            self.span.set_visible(False)
            self.parent.spectrum_span_button.setText('Show Wavelength Range')
        self.draw()

    def on_span_pick(self, event):
        self.figure.canvas.mpl_disconnect(self.cid[0])
        self.span.set_edgecolor('m')
        self.span.set_linewidth(5)
        self.draw()
        if event.edge == 'left':
            self.cid[1] = self.figure.canvas.mpl_connect('motion_notify_event', self.on_span_left_motion)
        elif event.edge == 'right':
            self.cid[1] = self.figure.canvas.mpl_connect('motion_notify_event', self.on_span_right_motion)
        self.cid[2] = self.figure.canvas.mpl_connect('button_press_event', self.on_span_resized)

    def on_span_left_motion(self, mouseevent):
        self.span.xy[0][0] = mouseevent.xdata
        self.span.xy[1][0] = mouseevent.xdata
        self.span.xy[4][0] = mouseevent.xdata
        self.draw()
        #self.parent.line_info[-1].update(self.span.xy[0][0], self.span.xy[2][0])

    def on_span_right_motion(self, mouseevent):
        self.span.xy[2][0] = mouseevent.xdata
        self.span.xy[3][0] = mouseevent.xdata
        self.draw()
        #self.parent.line_info[-1].update(self.span.xy[0][0], self.span.xy[2][0])

    def on_span_resized(self, mouseevent):
        self.figure.canvas.mpl_disconnect(self.cid[1])
        self.figure.canvas.mpl_disconnect(self.cid[2])
        self.cid[0] = self.figure.canvas.mpl_connect('pick_event', self.on_span_pick)
        self.span.set_edgecolor('r')
        self.span.set_linewidth(1)
        self.draw()
        self.parent.line_info.append(LineInfo(self.parent, self.span.xy[0][0], self.span.xy[2][0]))
        #self.parent.line_info[-1].show()

    def on_shell_pick(self, event):
        self.highlight_shell(event.artist.index)

    def highlight_shell(self, index):
        self.parent.tableview.selectRow(index)
        for i in range(len(self.parent.shells)):
            if i != index and i != index + 1:
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

    def span_picker(self, span, mouseevent, tolerance=5):
        left = float(span.xy[0][0])
        right = float(span.xy[2][0])
        tolerance = span.axes.transData.inverted().transform((tolerance, 0))[0] - span.axes.transData.inverted().transform((0, 0))[0]
        event_attributes = {'edge': None}
        if mouseevent.xdata is None:
            return False, event_attributes
        if left - tolerance <= mouseevent.xdata <= left + tolerance:
            event_attributes['edge'] = 'left'
            return True, event_attributes
        elif right - tolerance <= mouseevent.xdata <= right + tolerance:
            event_attributes['edge'] = 'right'
            return True, event_attributes
        return False, event_attributes

class Shell(matplotlib.patches.Wedge):

    def __init__(self, index, center, r_inner, r_outer, **kwargs):
        super(Shell, self).__init__(center, r_outer, 0, 90, width=r_outer - r_inner, **kwargs)
        self.index = index
        self.center = center
        self.r_outer = r_outer
        self.r_inner = r_inner
        self.width = r_outer - r_inner
