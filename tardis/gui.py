import numpy as np
import matplotlib
#matplotlib.use('KtAgg')
import matplotlib.pylab as plt
import matplotlib.gridspec as gridspec
from matplotlib import colors
from matplotlib.patches import Circle
from matplotlib.figure import *
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt4 import NavigationToolbar2QT as NavigationToolbar
import os
if os.environ.get('QT_API', None) is None:
    from PyQt4 import QtGui, QtCore
else:
    from PySide import QtGui, QtCore
    
from astropy import units as u
from tardis import analysis, util

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
        self.tablemodel = SimpleTableModel([['Shell: '], ["Rad. temp", "Ws"]], (1, 0))
        self.tableview = QtGui.QTableView()
        self.graph = MatplotlibWidget(self, 'model')
        self.graph_label = QtGui.QLabel('Select Property:')
        self.graph_button = QtGui.QToolButton()
        self.spectrum = MatplotlibWidget(self)
        self.spectrum_label = QtGui.QLabel('Select Spectrum:')
        self.spectrum_button = QtGui.QToolButton()
        self.spectrum_span_button = QtGui.QPushButton('Show Wavelength Range')
        self.spectrum_line_info_button = QtGui.QPushButton('Show Line Info')
        self.layout = QtGui.QHBoxLayout()
        self.graph_sublayout = QtGui.QVBoxLayout()
        self.graph_subsublayout = QtGui.QHBoxLayout()
        self.spectrum_sublayout = QtGui.QVBoxLayout()
        self.spectrum_subsublayout = QtGui.QHBoxLayout()

        self.tableview.setMinimumWidth(200)
        self.tableview.connect(self.tableview.verticalHeader(), QtCore.SIGNAL('sectionClicked(int)'), self.graph.highlight_shell)
        self.tableview.connect(self.tableview.verticalHeader(), QtCore.SIGNAL('sectionDoubleClicked(int)'),
                               self.on_header_double_clicked)
        self.graph_button.setText('Rad. temp')
        self.spectrum_button.setText('spec_flux_angstrom')
        self.graph_button.setPopupMode(QtGui.QToolButton.MenuButtonPopup)
        self.spectrum_button.setPopupMode(QtGui.QToolButton.MenuButtonPopup)
        self.graph_button.setMenu(QtGui.QMenu(self.graph_button))
        self.spectrum_button.setMenu(QtGui.QMenu(self.spectrum_button))
        self.graph_button.menu().addAction('Rad. temp').triggered.connect(self.change_graph_to_t_rads)
        self.graph_button.menu().addAction('Ws').triggered.connect(self.change_graph_to_ws)
        self.spectrum_button.menu().addAction('spec_flux_angstrom').triggered.connect(self.change_spectrum_to_spec_flux_angstrom)
        self.spectrum_button.menu().addAction('spec_virtual_flux_angstrom').triggered.connect(self.change_spectrum_to_spec_virtual_flux_angstrom)
        self.spectrum_span_button.clicked.connect(self.spectrum.show_span)
        self.spectrum_line_info_button.clicked.connect(self.spectrum.show_line_info)
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
        self.spectrum_sublayout.addWidget(self.spectrum_line_info_button)
        self.spectrum_sublayout.addWidget(self.spectrum)
        self.spectrum_sublayout.addWidget(self.spectrum.toolbar)
        self.layout.addLayout(self.spectrum_sublayout)
        self.spectrum_line_info_button.hide()
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
        self.tablemodel.addData(model.t_rads.value.tolist())
        self.tablemodel.addData(model.ws.tolist())

    def change_spectrum_to_spec_virtual_flux_angstrom(self):
        if self.model.spectrum_virtual.luminosity_density_lambda is None:
            luminosity_density_lambda = np.zeros_like(self.model.spectrum_virtual.wavelength)
        else:
            luminosity_density_lambda = self.model.spectrum_virtual.luminosity_density_lambda.value

        self.change_spectrum(luminosity_density_lambda, 'spec_flux_angstrom')

    def change_spectrum_to_spec_flux_angstrom(self):
        if self.model.spectrum.luminosity_density_lambda is None:
            luminosity_density_lambda = np.zeros_like(self.model.spectrum.wavelength)
        else:
            luminosity_density_lambda = self.model.spectrum.luminosity_density_lambda.value

        self.change_spectrum(luminosity_density_lambda, 'spec_flux_angstrom')

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
        wavelength = self.model.spectrum.wavelength.value
        if self.model.spectrum.luminosity_density_lambda is None:
            luminosity_density_lambda = np.zeros_like(wavelength)
        else:
            luminosity_density_lambda = self.model.spectrum.luminosity_density_lambda.value

        self.spectrum.dataplot = self.spectrum.ax.plot(wavelength, luminosity_density_lambda, label='b')
        self.spectrum.draw()

    def change_graph_to_ws(self):
        self.change_graph(self.model.ws, 'Ws', '')

    def change_graph_to_t_rads(self):
        self.change_graph(self.model.t_rads.value, 't_rads', '(K)')

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
        self.graph.ax1.set_title('Rad. Temp vs Shell')
        self.graph.ax1.set_xlabel('Shell Number')
        self.graph.ax1.set_ylabel('Rad. Temp (K)')
        self.graph.ax1.yaxis.get_major_formatter().set_powerlimits((0, 1))
        self.graph.dataplot = self.graph.ax1.plot(range(len(self.model.t_rads.value)), self.model.t_rads.value)
        self.graph.ax2.clear()
        self.graph.ax2.set_title('Shell View')
        self.graph.ax2.set_xlabel('Arbitrary')
        self.graph.ax2.set_ylabel('Arbitrary')
        self.shells = []
        t_rad_normalizer = colors.Normalize(vmin=self.model.t_rads.value.min(), vmax=self.model.t_rads.value.max())
        t_rad_color_map = plt.cm.ScalarMappable(norm=t_rad_normalizer, cmap=plt.cm.jet)
        t_rad_color_map.set_array(self.model.t_rads.value)
        if self.graph.cb:
            self.graph.cb.set_clim(vmin=self.model.t_rads.value.min(), vmax=self.model.t_rads.value.max())
            self.graph.cb.update_normal(t_rad_color_map)
        else:
            self.graph.cb = self.graph.figure.colorbar(t_rad_color_map)
            self.graph.cb.set_label('T (K)')
        self.graph.normalizing_factor = 0.2 * (self.model.tardis_config.structure.r_outer.value[-1] - self.model.tardis_config.structure.r_inner.value[0]) / self.model.tardis_config.structure.r_inner.value[0]
        #self.graph.normalizing_factor = 8e-16
        for i, t_rad in enumerate(self.model.t_rads.value):
            r_inner = self.model.tardis_config.structure.r_inner.value[i] * self.graph.normalizing_factor
            r_outer = self.model.tardis_config.structure.r_outer.value[i] * self.graph.normalizing_factor
            # This is the portion that has to be altered for the GUI to change from ellipses to circles.
            # However this needs a constant factor as seen in the below function for circles.
            #def circle (x,y,r) where x is r_inner and y is r_outer, I am not sure what constant has to be added amd 
                #random_object = self.shells.append(Shell(i, (0 + r,0 + r), r_inner -r, r_outer - r))
                #return random object
            self.shells.append(Shell(i, (0,0), r_inner, r_outer, facecolor=t_rad_color_map.to_rgba(t_rad),
                                     picker=self.graph.shell_picker))
            self.graph.ax2.add_patch(self.shells[i])
        self.graph.ax2.set_xlim(0, self.model.tardis_config.structure.r_outer.value[-1] * self.graph.normalizing_factor)
        self.graph.ax2.set_ylim(0, self.model.tardis_config.structure.r_outer.value[-1] * self.graph.normalizing_factor)
        self.graph.figure.tight_layout()
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
        self.atomstable.connect(self.atomstable.verticalHeader(), QtCore.SIGNAL('sectionClicked(int)'),
                               self.on_atom_header_double_clicked)


        self.table1_data = self.parent.model.tardis_config.abundances[self.shell_index]
        self.atomsdata = SimpleTableModel([['Z = '], ['Count (Shell %d)' % (self.shell_index + 1)]], iterate_header=(2, 0), index_info=self.table1_data.index.values.tolist())
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
        self.table2_data = self.parent.model.plasma_array.ion_populations[self.shell_index].ix[self.current_atom_index]
        self.ionsdata = SimpleTableModel([['Ion: '], ['Count (Z = %d)' % self.current_atom_index]], iterate_header=(2, 0), index_info=self.table2_data.index.values.tolist())
        normalized_data = []
        for item in self.table2_data.values:
            normalized_data.append(float(item /
                                   self.parent.model.tardis_config.number_densities[self.shell_index]
                                   .ix[self.current_atom_index]))


        self.ionsdata.addData(normalized_data)
        self.ionstable.setModel(self.ionsdata)
        self.ionstable.connect(self.ionstable.verticalHeader(), QtCore.SIGNAL('sectionClicked(int)'),
                               self.on_ion_header_double_clicked)
        self.levelstable.hide()
        self.ionstable.setColumnWidth(0, 120)
        self.ionstable.show()
        self.setGeometry(400, 150, 380, 400)
        self.show()

    def on_ion_header_double_clicked(self, index):
        self.current_ion_index = self.table2_data.index.values.tolist()[index]
        self.table3_data = self.parent.model.plasma_array.level_populations[self.shell_index].ix[self.current_atom_index,
                                                                                                 self.current_ion_index]
        self.levelsdata = SimpleTableModel([['Level: '], ['Count (Ion %d)' % self.current_ion_index]], iterate_header=(2, 0), index_info=self.table3_data.index.values.tolist())
        normalized_data = []
        for item in self.table3_data.values.tolist():
            normalized_data.append(float(item / self.table2_data.ix[self.current_ion_index]))
        self.levelsdata.addData(normalized_data)
        self.levelstable.setModel(self.levelsdata)
        self.levelstable.setColumnWidth(0, 120)
        self.levelstable.show()
        self.setGeometry(400, 150, 580, 400)
        self.show()

    def update_tables(self):
        self.table1_data = self.parent.model.plasma_array[self.shell_index].number_densities
        self.atomsdata.index_info=self.table1_data.index.values.tolist()
        self.atomsdata.arraydata = []
        self.atomsdata.addData(self.table1_data.values.tolist())
        self.atomsdata.updateTable()
        self.ionstable.hide()
        self.levelstable.hide()
        self.setGeometry(400, 150, 200, 400)
        self.show()


class LineInteractionTables(QtGui.QWidget):

    def __init__(self, line_interaction_analysis, atom_data, description):
        super(LineInteractionTables, self).__init__()
        self.text_description = QtGui.QLabel(str(description))
        self.species_table = QtGui.QTableView()
        self.transitions_table = QtGui.QTableView()
        self.layout = QtGui.QHBoxLayout()
        self.line_interaction_analysis = line_interaction_analysis
        self.atom_data = atom_data
        line_interaction_species_group = line_interaction_analysis.last_line_in.groupby(['atomic_number', 'ion_number'])
        self.species_selected = sorted(line_interaction_species_group.groups.keys())
        species_symbols = [util.species_tuple_to_string(item, atom_data) for item in self.species_selected]
        species_table_model = SimpleTableModel([species_symbols, ['Species']])
        species_abundances = (line_interaction_species_group.wavelength.count().astype(float) /
                                    line_interaction_analysis.last_line_in.wavelength.count()).astype(float).tolist()
        species_abundances = map(float, species_abundances)
        species_table_model.addData(species_abundances)
        self.species_table.setModel(species_table_model)

        line_interaction_species_group.wavelength.count()
        self.layout.addWidget(self.text_description)
        self.layout.addWidget(self.species_table)
        self.species_table.connect(self.species_table.verticalHeader(), QtCore.SIGNAL('sectionClicked(int)'),
                               self.on_species_clicked)
        self.layout.addWidget(self.transitions_table)

        self.setLayout(self.layout)
        self.show()

    def on_species_clicked(self, index):
        current_species = self.species_selected[index]
        last_line_in = self.line_interaction_analysis.last_line_in
        last_line_out = self.line_interaction_analysis.last_line_out

        last_line_in_filter = (last_line_in.atomic_number == current_species[0]).values & \
                              (last_line_in.ion_number == current_species[1]).values

        current_last_line_in = last_line_in[last_line_in_filter].reset_index()
        current_last_line_out = last_line_out[last_line_in_filter].reset_index()

        current_last_line_in['line_id_out'] = current_last_line_out['line_id']


        last_line_in_string = []
        last_line_count = []
        grouped_line_interactions = current_last_line_in.groupby(['line_id', 'line_id_out'])
        exc_deexc_string = 'exc. %d-%d (%.2f A) de-exc. %d-%d (%.2f A)'

        for line_id, row in grouped_line_interactions.wavelength.count().iteritems():
            current_line_in = self.atom_data.lines.ix[line_id[0]]
            current_line_out = self.atom_data.lines.ix[line_id[1]]
            last_line_in_string.append(exc_deexc_string % (current_line_in['level_number_lower'],
                                                           current_line_in['level_number_upper'],
                                                           current_line_in['wavelength'],
                                                           current_line_out['level_number_upper'],
                                                           current_line_out['level_number_lower'],
                                                           current_line_out['wavelength']))
            last_line_count.append(int(row))


        last_line_in_model = SimpleTableModel([last_line_in_string, ['Num. pkts %d' %
                                                                     current_last_line_in.wavelength.count()]])
        last_line_in_model.addData(last_line_count)
        self.transitions_table.setModel(last_line_in_model)



class LineInfo(QtGui.QDialog):

    def __init__(self, parent, wavelength_start, wavelength_end):
        super(LineInfo, self).__init__(parent)
        self.parent = parent
        self.setGeometry(180 + len(self.parent.line_info) * 20, 150, 250, 400)
        self.setWindowTitle('Line Interaction: %.2f - %.2f (A) ' % (wavelength_start, wavelength_end,
        ))
        self.layout = QtGui.QVBoxLayout()
        packet_nu_line_interaction = analysis.LastLineInteraction.from_model(self.parent.model)
        packet_nu_line_interaction.packet_filter_mode = 'packet_nu'
        packet_nu_line_interaction.wavelength_start = wavelength_start * u.angstrom
        packet_nu_line_interaction.wavelength_end = wavelength_end * u.angstrom
        
        line_in_nu_line_interaction = analysis.LastLineInteraction.from_model(self.parent.model)
        line_in_nu_line_interaction.packet_filter_mode = 'line_in_nu'
        line_in_nu_line_interaction.wavelength_start = wavelength_start * u.angstrom
        line_in_nu_line_interaction.wavelength_end = wavelength_end * u.angstrom


        self.layout.addWidget(LineInteractionTables(packet_nu_line_interaction, self.parent.model.atom_data, 'filtered by frequency of packet'))
        self.layout.addWidget(LineInteractionTables(line_in_nu_line_interaction, self.parent.model.atom_data, 'filtered by frequency of line interaction'))


        self.setLayout(self.layout)
        self.show()

    def get_data(self, wavelength_start, wavelength_end):
        self.wavelength_start = wavelength_start * u.angstrom
        self.wavelength_end = wavelength_end * u.angstrom
        last_line_in_ids, last_line_out_ids = analysis.get_last_line_interaction(self.wavelength_start, self.wavelength_end, self.parent.model)
        self.last_line_in, self.last_line_out = self.parent.model.atom_data.lines.ix[last_line_in_ids], self.parent.model.atom_data.lines.ix[last_line_out_ids]
        self.grouped_lines_in, self.grouped_lines_out = self.last_line_in.groupby(['atomic_number', 'ion_number']), self.last_line_out.groupby(['atomic_number', 'ion_number'])
        self.ions_in, self.ions_out = self.grouped_lines_in.groups.keys(), self.grouped_lines_out.groups.keys()
        self.ions_in.sort()
        self.ions_out.sort()
        self.header_list = []
        self.ion_table = (self.grouped_lines_in.wavelength.count().astype(float) / self.grouped_lines_in.wavelength.count().sum()).values.tolist()
        for z, ion in self.ions_in:
            self.header_list.append('Z = %d: Ion %d' % (z, ion))

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
            transitions_parsed.append("%d-%d (%.2f A)" % (key[0], key[1], self.parent.model.atom_data.lines.ix[value[0]]['wavelength']))
        return transitions_parsed, transitions_count

    def on_atom_clicked(self, index):
        self.transitionsin_parsed, self.transitionsin_count = self.get_transition_table(self.last_line_in, self.ions_in[index][0], self.ions_in[index][1])
        self.transitionsout_parsed, self.transitionsout_count = self.get_transition_table(self.last_line_out, self.ions_out[index][0], self.ions_out[index][1])
        self.transitionsindata = SimpleTableModel([self.transitionsin_parsed, ['Lines In']])
        self.transitionsoutdata = SimpleTableModel([self.transitionsout_parsed, ['Lines Out']])
        self.transitionsindata.addData(self.transitionsin_count)
        self.transitionsoutdata.addData(self.transitionsout_count)
        self.transitionsintable.setModel(self.transitionsindata)
        self.transitionsouttable.setModel(self.transitionsoutdata)
        self.transitionsintable.show()
        self.transitionsouttable.show()
        self.setGeometry(180 + len(self.parent.line_info) * 20, 150, 750, 400)
        self.show()

    def on_atom_clicked2(self, index):
        self.transitionsin_parsed, self.transitionsin_count = self.get_transition_table(self.last_line_in, self.ions_in[index][0], self.ions_in[index][1])
        self.transitionsout_parsed, self.transitionsout_count = self.get_transition_table(self.last_line_out, self.ions_out[index][0], self.ions_out[index][1])
        self.transitionsindata = SimpleTableModel([self.transitionsin_parsed, ['Lines In']])
        self.transitionsoutdata = SimpleTableModel([self.transitionsout_parsed, ['Lines Out']])
        self.transitionsindata.addData(self.transitionsin_count)
        self.transitionsoutdata.addData(self.transitionsout_count)
        self.transitionsintable2.setModel(self.transitionsindata)
        self.transitionsouttable2.setModel(self.transitionsoutdata)
        self.transitionsintable2.show()
        self.transitionsouttable2.show()
        self.setGeometry(180 + len(self.parent.line_info) * 20, 150, 750, 400)
        self.show()

class SimpleTableModel(QtCore.QAbstractTableModel):
    def __init__(self, headerdata=None, iterate_header=(0, 0), index_info=None, parent=None, *args):
        super(SimpleTableModel, self).__init__(parent, *args)
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
                return self.headerdata[0][0] + str(section + 1)
            elif self.iterate_header[0] == 2:
                if self.index_info:
                    return self.headerdata[0][0] + str(self.index_info[section])
                else:
                    return self.headerdata[0][0] + str(section + 1)
            else:
                return self.headerdata[0][section]
        elif orientation == QtCore.Qt.Horizontal and role == QtCore.Qt.DisplayRole:
            if self.iterate_header[1] == 1:
                return self.headerdata[1][0] + str(section + 1)
            elif self.iterate_header[1] == 2:
                if self.index_info:
                    return self.headerdata[1][0] + str(self.index_info[section])
            else:
                return self.headerdata[1][section]
        return None  #returning none instead of an empty string enables the gui to properly display the header files      

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
                self.setData(index, self.arraydata[c][r])

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

    def show_line_info(self):
        self.parent.line_info.append(LineInfo(self.parent, self.span.xy[0][0], self.span.xy[2][0]))

    def show_span(self, garbage=0, left=5000, right=10000):
        if self.parent.spectrum_span_button.text() == 'Show Wavelength Range':
            if not self.span:
                self.span = self.ax.axvspan(left, right, color='r', alpha=0.3, picker=self.span_picker)
            else:
                self.span.set_visible(True)
            self.parent.spectrum_line_info_button.show()
            self.parent.spectrum_span_button.setText('Hide Wavelength Range')
        else:
            self.span.set_visible(False)
            self.parent.spectrum_line_info_button.hide()
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
        if mouseevent.xdata < self.span.xy[2][0]:
            self.span.xy[0][0] = mouseevent.xdata
            self.span.xy[1][0] = mouseevent.xdata
            self.span.xy[4][0] = mouseevent.xdata
            self.draw()

    def on_span_right_motion(self, mouseevent):
        if mouseevent.xdata > self.span.xy[0][0]:
            self.span.xy[2][0] = mouseevent.xdata
            self.span.xy[3][0] = mouseevent.xdata
            self.draw()

    def on_span_resized(self, mouseevent):
        self.figure.canvas.mpl_disconnect(self.cid[1])
        self.figure.canvas.mpl_disconnect(self.cid[2])
        self.cid[0] = self.figure.canvas.mpl_connect('pick_event', self.on_span_pick)
        self.span.set_edgecolor('r')
        self.span.set_linewidth(1)
        self.draw()

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
