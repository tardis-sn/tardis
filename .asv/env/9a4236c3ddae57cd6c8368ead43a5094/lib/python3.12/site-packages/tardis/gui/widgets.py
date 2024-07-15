import os

import tardis.util.base

if os.environ.get("QT_API", None) == "pyqt":
    from PyQt5 import QtGui, QtCore, QtWidgets
elif os.environ.get("QT_API", None) == "pyside":
    from PySide2 import QtGui, QtCore, QtWidgets
else:
    raise ImportError(
        """QT_API was not set! Please exit the IPython console\n
         and at the bash prompt use : \n\n export QT_API=pyside \n or\n
         export QT_API=pyqt \n\n For more information refer to user guide."""
    )

import matplotlib
from matplotlib.figure import *
import matplotlib.gridspec as gridspec
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5 import (
    NavigationToolbar2QT as NavigationToolbar,
)
from matplotlib import colors
from matplotlib.patches import Circle
import matplotlib.pylab as plt
from astropy import units as u

import tardis
from tardis import analysis, util


class MatplotlibWidget(FigureCanvas):
    """Canvas to draw graphs on."""

    def __init__(self, tablecreator, parent, fig=None):
        """Create the canvas. Add toolbar depending on the parent."""

        # Force-deactivate LaTeX
        plt.rcParams["text.usetex"] = False

        self.tablecreator = tablecreator
        self.parent = parent
        self.figure = Figure()  # (frameon=False,facecolor=(1,1,1))
        self.cid = {}
        if fig != "model":
            self.ax = self.figure.add_subplot(111)
        else:
            self.gs = gridspec.GridSpec(2, 1, height_ratios=[1, 3])
            self.ax1 = self.figure.add_subplot(self.gs[0])
            self.ax2 = self.figure.add_subplot(self.gs[1])  # , aspect='equal')
        self.cb = None
        self.span = None

        super(MatplotlibWidget, self).__init__(self.figure)
        super(MatplotlibWidget, self).setSizePolicy(
            QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Expanding
        )
        super(MatplotlibWidget, self).updateGeometry()
        if fig != "model":
            self.toolbar = NavigationToolbar(self, parent)
            self.cid[0] = self.figure.canvas.mpl_connect(
                "pick_event", self.on_span_pick
            )
        else:
            self.cid[0] = self.figure.canvas.mpl_connect(
                "pick_event", self.on_shell_pick
            )

    def show_line_info(self):
        """Show line info for span selected region."""
        self.parent.line_info.append(
            LineInfo(
                self.parent,
                self.span.xy[0][0],
                self.span.xy[2][0],
                self.tablecreator,
            )
        )

    def show_span(self, garbage=0, left=5000, right=10000):
        """Hide/Show/Change the buttons that show line info
        in spectrum plot widget.
        """
        if self.parent.spectrum_span_button.text() == "Show Wavelength Range":
            if not self.span:
                self.span = self.ax.axvspan(
                    left, right, color="r", alpha=0.3, picker=self.span_picker
                )
            else:
                self.span.set_visible(True)
            self.parent.spectrum_line_info_button.show()
            self.parent.spectrum_span_button.setText("Hide Wavelength Range")
        else:
            self.span.set_visible(False)
            self.parent.spectrum_line_info_button.hide()
            self.parent.spectrum_span_button.setText("Show Wavelength Range")
        self.draw()

    def on_span_pick(self, event):
        """Callback to 'pick'(grab with mouse) the span selector tool."""
        self.figure.canvas.mpl_disconnect(self.cid[0])
        self.span.set_edgecolor("m")
        self.span.set_linewidth(5)
        self.draw()
        if event.edge == "left":
            self.cid[1] = self.figure.canvas.mpl_connect(
                "motion_notify_event", self.on_span_left_motion
            )
        elif event.edge == "right":
            self.cid[1] = self.figure.canvas.mpl_connect(
                "motion_notify_event", self.on_span_right_motion
            )
        self.cid[2] = self.figure.canvas.mpl_connect(
            "button_press_event", self.on_span_resized
        )

    def on_span_left_motion(self, mouseevent):
        """Update data of span selector tool on left movement of mouse and
        redraw.
        """
        if mouseevent.xdata < self.span.xy[2][0]:
            self.span.xy[0][0] = mouseevent.xdata
            self.span.xy[1][0] = mouseevent.xdata
            self.span.xy[4][0] = mouseevent.xdata
            self.draw()

    def on_span_right_motion(self, mouseevent):
        """Update data of span selector tool on right movement of mouse and
        redraw.
        """
        if mouseevent.xdata > self.span.xy[0][0]:
            self.span.xy[2][0] = mouseevent.xdata
            self.span.xy[3][0] = mouseevent.xdata
            self.draw()

    def on_span_resized(self, mouseevent):
        """Redraw the red rectangle to currently selected span."""
        self.figure.canvas.mpl_disconnect(self.cid[1])
        self.figure.canvas.mpl_disconnect(self.cid[2])
        self.cid[0] = self.figure.canvas.mpl_connect(
            "pick_event", self.on_span_pick
        )
        self.span.set_edgecolor("r")
        self.span.set_linewidth(1)
        self.draw()

    def on_shell_pick(self, event):
        """Highlight the shell that was picked."""
        self.highlight_shell(event.artist.index)

    def highlight_shell(self, index):
        """Change edgecolor of highlighted shell."""
        self.parent.tableview.selectRow(index)
        for i in range(len(self.parent.shells)):
            if i != index and i != index + 1:
                self.parent.shells[i].set_edgecolor("k")
            else:
                self.parent.shells[i].set_edgecolor("w")
        self.draw()

    def shell_picker(self, shell, mouseevent):
        """Enable picking shells in the shell plot."""
        if mouseevent.xdata is None:
            return False, dict()
        mouse_r2 = mouseevent.xdata**2 + mouseevent.ydata**2
        if shell.r_inner**2 < mouse_r2 < shell.r_outer**2:
            return True, dict()
        return False, dict()

    def span_picker(self, span, mouseevent, tolerance=5):
        """Detect mouseclicks inside tolerance region of the span selector
        tool and pick it.
        """
        left = float(span.xy[0][0])
        right = float(span.xy[2][0])
        tolerance = (
            span.axes.transData.inverted().transform((tolerance, 0))[0]
            - span.axes.transData.inverted().transform((0, 0))[0]
        )
        event_attributes = {"edge": None}
        if mouseevent.xdata is None:
            return False, event_attributes
        if left - tolerance <= mouseevent.xdata <= left + tolerance:
            event_attributes["edge"] = "left"
            return True, event_attributes
        elif right - tolerance <= mouseevent.xdata <= right + tolerance:
            event_attributes["edge"] = "right"
            return True, event_attributes
        return False, event_attributes


class Shell(matplotlib.patches.Wedge):
    """A data holder to store measurements of shells that will be drawn in
    the plot.
    """

    def __init__(self, index, center, r_inner, r_outer, **kwargs):
        super(Shell, self).__init__(
            center, r_outer, 0, 90, width=r_outer - r_inner, **kwargs
        )
        self.index = index
        self.center = center
        self.r_outer = r_outer
        self.r_inner = r_inner
        self.width = r_outer - r_inner


class ConfigEditor(QtWidgets.QWidget):
    """The configuration editor widget.

    This widget is added to the stacked widget that is the central widget of
    the main top level window created by Tardis.
    """

    def __init__(self, yamlconfigfile, parent=None):
        """Create and return the configuration widget.

        Parameters
        ----------
        yamlconfigfile : string
            File name of the yaml configuration file.
        parent : None
            Set to None. The parent is changed when the widget is
            appended to the layout of its parent.
        """
        super(ConfigEditor, self).__init__(parent)

        # Configurations from the input and template
        configDict = yaml.load(open(yamlconfigfile), Loader=yaml.CLoader)
        templatedictionary = {
            "tardis_config_version": [True, "v1.0"],
            "supernova": {
                "luminosity_requested": [True, "1 solLum"],
                "time_explosion": [True, None],
                "distance": [False, None],
                "luminosity_wavelength_start": [False, "0 angstrom"],
                "luminosity_wavelength_end": [False, "inf angstrom"],
            },
            "atom_data": [True, "File Browser"],
            "plasma": {
                "initial_t_inner": [False, "-1K"],
                "initial_t_rad": [False, "10000K"],
                "disable_electron_scattering": [False, False],
                "ionization": [True, None],
                "excitation": [True, None],
                "radiative_rates_type": [True, None],
                "line_interaction_type": [True, None],
                "w_epsilon": [False, 1e-10],
                "delta_treatment": [False, None],
                "nlte": {
                    "species": [False, []],
                    "coronal_approximation": [False, False],
                    "classical_nebular": [False, False],
                },
            },
            "model": {
                "structure": {
                    "type": [
                        True,
                        [
                            "file|_:_|filename|_:_|"
                            "filetype|_:_|v_inner_boundary|_:_|v_outer_boundary",
                            "specific|_:_|velocity|_:_|density",
                        ],
                    ],
                    "filename": [True, None],
                    "filetype": [True, None],
                    "v_inner_boundary": [False, "0 km/s"],
                    "v_outer_boundary": [False, "inf km/s"],
                    "velocity": [True, None],
                    "density": {
                        "type": [
                            True,
                            [
                                "branch85_w7|_:_|w7_time_0"
                                "|_:_|w7_time_0|_:_|w7_time_0",
                                "exponential|_:_|time_0|_:_|rho_0|_:_|" "v_0",
                                "power_law|_:_|time_0|_:_|rho_0"
                                "|_:_|v_0|_:_|exponent",
                                "uniform|_:_|value",
                            ],
                        ],
                        "w7_time_0": [False, "0.000231481 day"],
                        "w7_rho_0": [False, "3e29 g/cm^3"],
                        "w7_v_0": [False, "1 km/s"],
                        "time_0": [True, None],
                        "rho_0": [True, None],
                        "v_0": [True, None],
                        "exponent": [True, None],
                        "value": [True, None],
                    },
                },
                "abundances": {
                    "type": [
                        True,
                        ["file|_:_|filetype|_:_|" "filename", "uniform"],
                    ],
                    "filename": [True, None],
                    "filetype": [False, None],
                },
            },
            "montecarlo": {
                "seed": [False, 23111963],
                "no_of_packets": [True, None],
                "iterations": [True, None],
                "black_body_sampling": {
                    "start": "1 angstrom",
                    "stop": "1000000 angstrom",
                    "num": "1.e+6",
                },
                "last_no_of_packets": [False, -1],
                "no_of_virtual_packets": [False, 0],
                "convergence_strategy": {
                    "type": [
                        True,
                        [
                            "damped|_:_|damping_constant|_:_|t_inner|_:_|"
                            "t_rad|_:_|w|_:_|lock_t_inner_cycles|_:_|"
                            "t_inner_update_exponent",
                            "specific|_:_|threshold"
                            "|_:_|fraction|_:_|hold_iterations|_:_|t_inner"
                            "|_:_|t_rad|_:_|w|_:_|lock_t_inner_cycles|_:_|"
                            "damping_constant|_:_|t_inner_update_exponent",
                        ],
                    ],
                    "t_inner_update_exponent": [False, -0.5],
                    "lock_t_inner_cycles": [False, 1],
                    "hold_iterations": [True, 3],
                    "fraction": [True, 0.8],
                    "damping_constant": [False, 0.5],
                    "threshold": [True, None],
                    "t_inner": {
                        "damping_constant": [False, 0.5],
                        "threshold": [False, None],
                    },
                    "t_rad": {
                        "damping_constant": [False, 0.5],
                        "threshold": [True, None],
                    },
                    "w": {
                        "damping_constant": [False, 0.5],
                        "threshold": [True, None],
                    },
                },
            },
            "spectrum": [True, None],
        }
        self.match_dicts(configDict, templatedictionary)

        self.layout = QtWidgets.QVBoxLayout()

        # Make tree
        self.trmodel = TreeModel(templatedictionary)
        self.colView = QtWidgets.QColumnView()
        self.colView.setModel(self.trmodel)
        # Five columns of width 256 each can be visible at once
        self.colView.setFixedWidth(256 * 5)
        self.colView.setItemDelegate(TreeDelegate(self))
        self.layout.addWidget(self.colView)

        # Recalculate button
        button = QtWidgets.QPushButton("Recalculate")
        button.setFixedWidth(90)
        self.layout.addWidget(button)
        button.clicked.connect(self.recalculate)

        # Finally put them all in
        self.setLayout(self.layout)

    def match_dicts(self, dict1, dict2):  # dict1<=dict2
        """Compare and combine two dictionaries.

        If there are new keys in `dict1` then they are appended to `dict2`.
        The `dict2` stores the values for the keys in `dict1` but it
        first modifies them by taking the value and appending to a
        list whose first item is either True or False, indicating if
        that key is mandatory or not. The goal of this method is to
        perform validation by inserting user provided values into the
        template dictionary. After inserting user given values into the
        template dictionary, all the keys have either default or user
        provided values. Then it can be used to build a tree to be
        shown in the ConfigEditor.

        Parameters
        ----------
        dict1 : dictionary
            The dictionary of user provided configuration.
        dict2 : dictionary
            The template dictionary with all default values set. This
            one may have some keys missing that are present in the
            `dict1`. Such keys will be appended.

        Raises
        ------
        IOError
            If the configuration file has an invalid value for a
            key that can only take values from a predefined list,
            then this error is raised.
        """
        for key in dict1:
            if key in dict2:
                if isinstance(dict2[key], dict):
                    self.match_dicts(dict1[key], dict2[key])

                elif isinstance(dict2[key], list):
                    if isinstance(dict2[key][1], list):
                        # options = dict2[key][1] #This is passed by reference.
                        # So copy the list manually.
                        options = [
                            dict2[key][1][i] for i in range(len(dict2[key][1]))
                        ]

                        for i in range(len(options)):
                            options[i] = options[i].split("|_:_|")[0]

                        optionselected = dict1[key]

                        if optionselected in options:
                            indexofselected = options.index(optionselected)
                            temp = dict2[key][1][0]

                            dict2[key][1][0] = dict2[key][1][indexofselected]
                            dict2[key][1][indexofselected] = temp

                        else:
                            print("The selected and available options")
                            print(optionselected)
                            print(options)
                            raise IOError(
                                "An invalid option was"
                                " provided in the input file"
                            )

                else:
                    dict2[key] = dict1[key]
            else:
                toappend = [False, dict1[key]]
                dict2[key] = toappend

    def recalculate(self):
        """Recalculate and display the model from the modified data in
        the ConfigEditor.
        """
        pass


class ModelViewer(QtWidgets.QWidget):
    """The widget that holds all the plots and tables that visualize
    the data in the tardis model. This is also appended to the stacked
    widget in the top level window.
    """

    def __init__(self, tablecreator, parent=None):
        """Create all widgets that are children of ModelViewer."""
        QtWidgets.QWidget.__init__(self, parent)

        # Data structures
        self.model = None
        self.shell_info = {}
        self.line_info = []

        # functions
        self.createTable = tablecreator

        # Shells widget
        self.shellWidget = self.make_shell_widget()

        # Spectrum widget
        self.spectrumWidget = self.make_spectrum_widget()

        # Plot tab widget
        self.plotTabWidget = QtWidgets.QTabWidget()
        self.plotTabWidget.addTab(self.shellWidget, "&Shells")
        self.plotTabWidget.addTab(self.spectrumWidget, "S&pectrum")

        # Table widget
        self.tablemodel = self.createTable(
            [["Shell: "], ["Rad. temp", "Ws", "V"]], (1, 0)
        )
        self.tableview = QtWidgets.QTableView()
        self.tableview.setMinimumWidth(200)

        self.sectionClicked = QtCore.Signal(int)
        self.tableview.verticalHeader().sectionClicked.connect(
            self.graph.highlight_shell
        )

        self.sectionDoubleClicked = QtCore.Signal(int)
        self.tableview.verticalHeader().sectionDoubleClicked.connect(
            self.on_header_double_clicked
        )

        # Label for text output
        self.outputLabel = QtWidgets.QLabel()
        self.outputLabel.setFrameStyle(
            QtWidgets.QFrame.StyledPanel | QtWidgets.QFrame.Sunken
        )
        self.outputLabel.setStyleSheet("QLabel{background-color:white;}")

        # Group boxes
        graphsBox = QtWidgets.QGroupBox("Visualized results")
        textsBox = QtWidgets.QGroupBox("Model parameters")
        tableBox = QtWidgets.QGroupBox("Tabulated results")

        # For textbox
        textlayout = QtWidgets.QHBoxLayout()
        textlayout.addWidget(self.outputLabel)

        tableslayout = QtWidgets.QVBoxLayout()
        tableslayout.addWidget(self.tableview)
        tableBox.setLayout(tableslayout)

        visualayout = QtWidgets.QVBoxLayout()
        visualayout.addWidget(self.plotTabWidget)
        graphsBox.setLayout(visualayout)

        self.layout = QtWidgets.QHBoxLayout()
        self.layout.addWidget(graphsBox)
        textntablelayout = QtWidgets.QVBoxLayout()
        textsBox.setLayout(textlayout)
        textntablelayout.addWidget(textsBox)
        textntablelayout.addWidget(tableBox)

        self.layout.addLayout(textntablelayout)
        self.setLayout(self.layout)

    def fill_output_label(self):
        """Read some data from tardis model and display on the label for
        quick user access.
        """

        model_converged = (
            '<font color="green"><b>True</b></font>'
            if self.model.converged
            else '<font color="red"><b>False</b></font>'
        )

        labeltext = f"Iterations requested: {self.model.iterations} <br/> Iterations executed:  {self.model.iterations_executed}<br/>\
                     Model converged     : {model_converged} <br/> Simulation Time    :  {self.model.transport.time_of_simulation.value} s <br/>\
                     Inner Temperature   : {self.model.model.t_inner.value} K <br/> Number of packets  :  {self.model.last_no_of_packets}<br/>\
                     Inner Luminosity    : {self.model.transport.calculate_luminosity_inner(self.model.model)}"

        self.outputLabel.setText(labeltext)

    def make_shell_widget(self):
        """Create the plot of the the shells and place it inside a
        container widget. Return the container widget.
        """
        # Widgets for plot of shells
        self.graph = MatplotlibWidget(self.createTable, self, "model")
        self.graph_label = QtWidgets.QLabel("Select Property:")
        self.graph_button = QtWidgets.QToolButton()
        self.graph_button.setText("Rad. temp")
        self.graph_button.setPopupMode(QtWidgets.QToolButton.MenuButtonPopup)
        self.graph_button.setMenu(QtWidgets.QMenu(self.graph_button))
        self.graph_button.menu().addAction("Rad. temp").triggered.connect(
            self.change_graph_to_t_rads
        )
        self.graph_button.menu().addAction("Ws").triggered.connect(
            self.change_graph_to_ws
        )

        # Layouts: bottom up
        self.graph_subsublayout = QtWidgets.QHBoxLayout()
        self.graph_subsublayout.addWidget(self.graph_label)
        self.graph_subsublayout.addWidget(self.graph_button)

        self.graph_sublayout = QtWidgets.QVBoxLayout()
        self.graph_sublayout.addLayout(self.graph_subsublayout)
        self.graph_sublayout.addWidget(self.graph)

        containerWidget = QtWidgets.QWidget()
        containerWidget.setLayout(self.graph_sublayout)
        return containerWidget

    def make_spectrum_widget(self):
        """Create the spectrum plot and associated buttons and append to
        a container widget. Return the container widget.
        """
        self.spectrum = MatplotlibWidget(self.createTable, self)
        self.spectrum_label = QtWidgets.QLabel("Select Spectrum:")
        self.spectrum_button = QtWidgets.QToolButton()
        self.spectrum_button.setText("spec_flux_angstrom")
        self.spectrum_button.setPopupMode(QtWidgets.QToolButton.MenuButtonPopup)
        self.spectrum_button.setMenu(QtWidgets.QMenu(self.spectrum_button))
        self.spectrum_button.menu().addAction(
            "spec_flux_angstrom"
        ).triggered.connect(self.change_spectrum_to_spec_flux_angstrom)
        self.spectrum_button.menu().addAction(
            "spec_virtual_flux_angstrom"
        ).triggered.connect(self.change_spectrum_to_spec_virtual_flux_angstrom)
        self.spectrum_span_button = QtWidgets.QPushButton(
            "Show Wavelength Range"
        )
        self.spectrum_span_button.clicked.connect(self.spectrum.show_span)
        self.spectrum_line_info_button = QtWidgets.QPushButton("Show Line Info")
        self.spectrum_line_info_button.hide()
        self.spectrum_line_info_button.clicked.connect(
            self.spectrum.show_line_info
        )

        self.spectrum_subsublayout = QtWidgets.QHBoxLayout()
        self.spectrum_subsublayout.addWidget(self.spectrum_span_button)
        self.spectrum_subsublayout.addWidget(self.spectrum_label)
        self.spectrum_subsublayout.addWidget(self.spectrum_button)

        self.spectrum_sublayout = QtWidgets.QVBoxLayout()
        self.spectrum_sublayout.addLayout(self.spectrum_subsublayout)
        self.spectrum_sublayout.addWidget(self.spectrum_line_info_button)
        self.spectrum_sublayout.addWidget(self.spectrum)
        self.spectrum_sublayout.addWidget(self.spectrum.toolbar)

        containerWidget = QtWidgets.QWidget()
        containerWidget.setLayout(self.spectrum_sublayout)
        return containerWidget

    def update_data(self, model=None):
        """Associate the given model with the GUI and display results."""
        if model:
            self.change_model(model)
        self.tablemodel.update_table()
        for index in self.shell_info.keys():
            self.shell_info[index].update_tables()
        self.plot_model()
        if self.graph_button.text == "Ws":
            self.change_graph_to_ws()
        self.plot_spectrum()
        if self.spectrum_button.text == "spec_virtual_flux_angstrom":
            self.change_spectrum_to_spec_virtual_flux_angstrom()
        self.show()

    def change_model(self, model):
        """Reset the model set in the GUI."""
        self.model = model
        self.tablemodel.arraydata = []
        self.tablemodel.add_data(model.model.t_rad.value.tolist())
        self.tablemodel.add_data(model.model.w.tolist())
        self.tablemodel.add_data(model.model.velocity.value.tolist())

    def change_spectrum_to_spec_virtual_flux_angstrom(self):
        """Change the spectrum data to the virtual spectrum."""
        if (
            self.model.transport.spectrum_virtual.luminosity_density_lambda
            is None
        ):
            luminosity_density_lambda = np.zeros_like(
                self.model.transport.spectrum_virtual.wavelength
            )
        else:
            luminosity_density_lambda = (
                self.model.transport.spectrum_virtual.luminosity_density_lambda.value
            )

        self.change_spectrum(luminosity_density_lambda, "spec_flux_angstrom")

    def change_spectrum_to_spec_flux_angstrom(self):
        """Change spectrum data back from virtual spectrum. (See the
        method above)."""
        if self.model.transport.spectrum.luminosity_density_lambda is None:
            luminosity_density_lambda = np.zeros_like(
                self.model.transport.spectrum.wavelength
            )
        else:
            luminosity_density_lambda = (
                self.model.transport.spectrum.luminosity_density_lambda.value
            )

        self.change_spectrum(luminosity_density_lambda, "spec_flux_angstrom")

    def change_spectrum(self, data, name):
        """Replot the spectrum plot using the data provided. Called
        when changing spectrum types. See the two methods above.
        """
        self.spectrum_button.setText(name)
        self.spectrum.dataplot[0].set_ydata(data)
        self.spectrum.ax.relim()
        self.spectrum.ax.autoscale()
        self.spectrum.draw()

    def plot_spectrum(self):
        """Plot the spectrum and add labels to the graph."""
        self.spectrum.ax.clear()
        self.spectrum.ax.set_title("Spectrum")
        self.spectrum.ax.set_xlabel("Wavelength (A)")
        self.spectrum.ax.set_ylabel("Intensity")
        wavelength = self.model.transport.spectrum.wavelength.value
        if self.model.transport.spectrum.luminosity_density_lambda is None:
            luminosity_density_lambda = np.zeros_like(wavelength)
        else:
            luminosity_density_lambda = (
                self.model.transport.spectrum.luminosity_density_lambda.value
            )

        self.spectrum.dataplot = self.spectrum.ax.plot(
            wavelength, luminosity_density_lambda, label="b"
        )
        self.spectrum.draw()

    def change_graph_to_ws(self):
        """Change the shell plot to show dilution factor."""
        self.change_graph(self.model.model.w, "Ws", "")

    def change_graph_to_t_rads(self):
        """Change the graph back to radiation Temperature."""
        self.change_graph(self.model.model.t_rad.value, "t_rad", "(K)")

    def change_graph(self, data, name, unit):
        """Called to change the shell plot by the two methods above."""
        self.graph_button.setText(name)
        self.graph.dataplot[0].set_ydata(data)
        self.graph.ax1.relim()
        self.graph.ax1.autoscale()
        self.graph.ax1.set_title(name + " vs Shell")
        self.graph.ax1.set_ylabel(name + " " + unit)
        normalizer = colors.Normalize(vmin=data.min(), vmax=data.max())
        color_map = plt.cm.ScalarMappable(norm=normalizer, cmap=plt.cm.jet)
        color_map.set_array(data)
        self.graph.cb.set_clim(vmin=data.min(), vmax=data.max())
        self.graph.cb.update_normal(color_map)
        if unit == "(K)":
            unit = "T (K)"
        self.graph.cb.set_label(unit)
        for i, item in enumerate(data):
            self.shells[i].set_facecolor(color_map.to_rgba(item))
        self.graph.draw()

    def plot_model(self):
        """Plot the two graphs, the shell model and the line plot
        both showing the radiation temperature and set labels.
        """
        self.graph.ax1.clear()
        self.graph.ax1.set_title("Rad. Temp vs Shell")
        self.graph.ax1.set_xlabel("Shell Number")
        self.graph.ax1.set_ylabel("Rad. Temp (K)")
        self.graph.ax1.yaxis.get_major_formatter().set_powerlimits((0, 1))
        self.graph.dataplot = self.graph.ax1.plot(
            range(len(self.model.model.t_rad.value)),
            self.model.model.t_rad.value,
        )
        self.graph.ax2.clear()
        self.graph.ax2.set_title("Shell View")
        self.graph.ax2.set_xticklabels([])
        self.graph.ax2.set_yticklabels([])
        self.graph.ax2.grid = True

        self.shells = []
        t_rad_normalizer = colors.Normalize(
            vmin=self.model.model.t_rad.value.min(),
            vmax=self.model.model.t_rad.value.max(),
        )
        t_rad_color_map = plt.cm.ScalarMappable(
            norm=t_rad_normalizer, cmap=plt.cm.jet
        )
        t_rad_color_map.set_array(self.model.model.t_rad.value)
        if self.graph.cb:
            self.graph.cb.set_clim(
                vmin=self.model.model.t_rad.value.min(),
                vmax=self.model.model.t_rad.value.max(),
            )
            self.graph.cb.update_normal(t_rad_color_map)
        else:
            self.graph.cb = self.graph.figure.colorbar(t_rad_color_map)
            self.graph.cb.set_label("T (K)")
        self.graph.normalizing_factor = (
            0.2
            * (
                self.model.model.r_outer.value[-1]
                - self.model.model.r_inner.value[0]
            )
            / (self.model.model.r_inner.value[0])
        )

        # self.graph.normalizing_factor = 8e-16
        for i, t_rad in enumerate(self.model.model.t_rad.value):
            r_inner = (
                self.model.model.r_inner.value[i]
                * self.graph.normalizing_factor
            )
            r_outer = (
                self.model.model.r_outer.value[i]
                * self.graph.normalizing_factor
            )
            self.shells.append(
                Shell(
                    i,
                    (0, 0),
                    r_inner,
                    r_outer,
                    facecolor=t_rad_color_map.to_rgba(t_rad),
                    picker=self.graph.shell_picker,
                )
            )
            self.graph.ax2.add_patch(self.shells[i])
        self.graph.ax2.set_xlim(
            0,
            self.model.model.r_outer.value[-1] * self.graph.normalizing_factor,
        )
        self.graph.ax2.set_ylim(
            0,
            self.model.model.r_outer.value[-1] * self.graph.normalizing_factor,
        )
        self.graph.figure.tight_layout()
        self.graph.draw()

    def on_header_double_clicked(self, index):
        """Callback to get counts for different Z from table."""
        self.shell_info[index] = ShellInfo(index, self.createTable, self)


class ShellInfo(QtWidgets.QDialog):
    """Dialog to display Shell abundances."""

    def __init__(self, index, tablecreator, parent=None):
        """Create the widget to display shell info and set data."""
        super(ShellInfo, self).__init__(parent)

        self.createTable = tablecreator
        self.parent = parent
        self.shell_index = index
        self.setGeometry(400, 150, 200, 400)
        self.setWindowTitle(f"Shell {self.shell_index + 1} Abundances")
        self.atomstable = QtWidgets.QTableView()
        self.ionstable = QtWidgets.QTableView()
        self.levelstable = QtWidgets.QTableView()
        self.sectionClicked = QtCore.Signal(int)
        self.atomstable.verticalHeader().sectionClicked.connect(
            self.on_atom_header_double_clicked
        )

        self.table1_data = self.parent.model.plasma.abundance[self.shell_index]
        self.atomsdata = self.createTable(
            [["Z = "], [f"Count (Shell {self.shell_index + 1})"]],
            iterate_header=(2, 0),
            index_info=self.table1_data.index.values.tolist(),
        )
        self.ionsdata = None
        self.levelsdata = None
        self.atomsdata.add_data(self.table1_data.values.tolist())
        self.atomstable.setModel(self.atomsdata)

        self.layout = QtWidgets.QHBoxLayout()
        self.layout.addWidget(self.atomstable)
        self.layout.addWidget(self.ionstable)
        self.layout.addWidget(self.levelstable)
        self.setLayout(self.layout)
        self.ionstable.hide()
        self.levelstable.hide()
        self.show()

    def on_atom_header_double_clicked(self, index):
        """Called when a header in the first column is clicked to show
        ion populations."""
        self.current_atom_index = self.table1_data.index.values.tolist()[index]
        self.table2_data = self.parent.model.plasma.ion_number_density[
            self.shell_index
        ].ix[self.current_atom_index]
        self.ionsdata = self.createTable(
            [["Ion: "], [f"Count (Z = {self.current_atom_index})"]],
            iterate_header=(2, 0),
            index_info=self.table2_data.index.values.tolist(),
        )
        normalized_data = []
        for item in self.table2_data.values:
            normalized_data.append(
                float(
                    item
                    / self.parent.model.plasma.number_density[
                        self.shell_index
                    ].ix[self.current_atom_index]
                )
            )

        self.ionsdata.add_data(normalized_data)
        self.ionstable.setModel(self.ionsdata)
        self.sectionClicked = QtCore.Signal(int)
        self.ionstable.verticalHeader().sectionClicked.connect(
            self.on_ion_header_double_clicked
        )
        self.levelstable.hide()
        self.ionstable.setColumnWidth(0, 120)
        self.ionstable.show()
        self.setGeometry(400, 150, 380, 400)
        self.show()

    def on_ion_header_double_clicked(self, index):
        """Called on double click of ion headers to show level populations."""
        self.current_ion_index = self.table2_data.index.values.tolist()[index]
        self.table3_data = self.parent.model.plasma.level_number_density[
            self.shell_index
        ].ix[self.current_atom_index, self.current_ion_index]
        self.levelsdata = self.createTable(
            [["Level: "], [f"Count (Ion {self.current_ion_index})"]],
            iterate_header=(2, 0),
            index_info=self.table3_data.index.values.tolist(),
        )
        normalized_data = []
        for item in self.table3_data.values.tolist():
            normalized_data.append(
                float(item / self.table2_data.ix[self.current_ion_index])
            )
        self.levelsdata.add_data(normalized_data)
        self.levelstable.setModel(self.levelsdata)
        self.levelstable.setColumnWidth(0, 120)
        self.levelstable.show()
        self.setGeometry(400, 150, 580, 400)
        self.show()

    def update_tables(self):
        """Update table data for shell info viewer."""
        self.table1_data = self.parent.model.plasma.number_density[
            self.shell_index
        ]
        self.atomsdata.index_info = self.table1_data.index.values.tolist()
        self.atomsdata.arraydata = []
        self.atomsdata.add_data(self.table1_data.values.tolist())
        self.atomsdata.update_table()
        self.ionstable.hide()
        self.levelstable.hide()
        self.setGeometry(400, 150, 200, 400)
        self.show()


class LineInfo(QtWidgets.QDialog):
    """Dialog to show the line info used by spectrum widget."""

    def __init__(self, parent, wavelength_start, wavelength_end, tablecreator):
        """Create the dialog and set data in it from the model.
        Show widget."""
        super(LineInfo, self).__init__(parent)
        self.createTable = tablecreator
        self.parent = parent
        self.setGeometry(180 + len(self.parent.line_info) * 20, 150, 250, 400)
        self.setWindowTitle(
            f"Line Interaction: {wavelength_start:.2f} - {wavelength_end:.2f} (A) "
        )
        self.layout = QtWidgets.QVBoxLayout()
        packet_nu_line_interaction = (
            analysis.LastLineInteraction.from_simulation(self.parent.model)
        )
        packet_nu_line_interaction.packet_filter_mode = "packet_nu"
        packet_nu_line_interaction.wavelength_start = (
            wavelength_start * u.angstrom
        )
        packet_nu_line_interaction.wavelength_end = wavelength_end * u.angstrom

        line_in_nu_line_interaction = (
            analysis.LastLineInteraction.from_simulation(self.parent.model)
        )
        line_in_nu_line_interaction.packet_filter_mode = "line_in_nu"
        line_in_nu_line_interaction.wavelength_start = (
            wavelength_start * u.angstrom
        )
        line_in_nu_line_interaction.wavelength_end = wavelength_end * u.angstrom

        self.layout.addWidget(
            LineInteractionTables(
                packet_nu_line_interaction,
                self.parent.model.plasma.atomic_data.atom_data,
                self.parent.model.plasma.lines,
                "filtered by frequency of packet",
                self.createTable,
            )
        )
        self.layout.addWidget(
            LineInteractionTables(
                line_in_nu_line_interaction,
                self.parent.model.plasma.atomic_data.atom_data,
                self.parent.model.plasma.lines,
                "filtered by frequency of line interaction",
                self.createTable,
            )
        )

        self.setLayout(self.layout)
        self.show()

    def get_data(self, wavelength_start, wavelength_end):
        """Fetch line info data for the specified wavelength range
        from the model and create ionstable.
        """
        self.wavelength_start = wavelength_start * u.angstrom
        self.wavelength_end = wavelength_end * u.angstrom
        (
            last_line_in_ids,
            last_line_out_ids,
        ) = analysis.get_last_line_interaction(
            self.wavelength_start, self.wavelength_end, self.parent.model
        )
        self.last_line_in, self.last_line_out = (
            self.parent.model.atom_data.lines.ix[last_line_in_ids],
            self.parent.model.atom_data.lines.ix[last_line_out_ids],
        )
        self.grouped_lines_in, self.grouped_lines_out = (
            self.last_line_in.groupby(["atomic_number", "ion_number"]),
            self.last_line_out.groupby(["atomic_number", "ion_number"]),
        )
        self.ions_in, self.ions_out = (
            self.grouped_lines_in.groups.keys(),
            self.grouped_lines_out.groups.keys(),
        )
        self.ions_in.sort()
        self.ions_out.sort()
        self.header_list = []
        self.ion_table = (
            self.grouped_lines_in.wavelength.count().astype(float)
            / self.grouped_lines_in.wavelength.count().sum()
        ).values.tolist()
        for z, ion in self.ions_in:
            self.header_list.append(f"Z = {z}: Ion {ion}")

    def get_transition_table(self, lines, atom, ion):
        """Called by the two methods below to get transition table for
        given lines, atom and ions.

        """
        grouped = lines.groupby(["atomic_number", "ion_number"])
        transitions_with_duplicates = (
            lines.ix[grouped.groups[(atom, ion)]]
            .groupby(["level_number_lower", "level_number_upper"])
            .groups
        )
        transitions = (
            lines.ix[grouped.groups[(atom, ion)]]
            .drop_duplicates()
            .groupby(["level_number_lower", "level_number_upper"])
            .groups
        )
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
            transitions_parsed.append(
                f"{key[0]}-{key[1]}"
                f"{self.parent.model.atom_data.lines.ix[value[0]]['wavelength']:.2f} A"
            )
        return transitions_parsed, transitions_count

    def on_atom_clicked(self, index):
        """Create and show transition table for the clicked item in the
        dialog created by the spectrum widget.
        """
        (
            self.transitionsin_parsed,
            self.transitionsin_count,
        ) = self.get_transition_table(
            self.last_line_in, self.ions_in[index][0], self.ions_in[index][1]
        )
        (
            self.transitionsout_parsed,
            self.transitionsout_count,
        ) = self.get_transition_table(
            self.last_line_out, self.ions_out[index][0], self.ions_out[index][1]
        )
        self.transitionsindata = self.createTable(
            [self.transitionsin_parsed, ["Lines In"]]
        )
        self.transitionsoutdata = self.createTable(
            [self.transitionsout_parsed, ["Lines Out"]]
        )
        self.transitionsindata.add_data(self.transitionsin_count)
        self.transitionsoutdata.add_data(self.transitionsout_count)
        self.transitionsintable.setModel(self.transitionsindata)
        self.transitionsouttable.setModel(self.transitionsoutdata)
        self.transitionsintable.show()
        self.transitionsouttable.show()
        self.setGeometry(180 + len(self.parent.line_info) * 20, 150, 750, 400)
        self.show()

    def on_atom_clicked2(self, index):
        """Create and show transition table for the clicked item in the
        dialog created by the spectrum widget.
        """
        (
            self.transitionsin_parsed,
            self.transitionsin_count,
        ) = self.get_transition_table(
            self.last_line_in, self.ions_in[index][0], self.ions_in[index][1]
        )
        (
            self.transitionsout_parsed,
            self.transitionsout_count,
        ) = self.get_transition_table(
            self.last_line_out, self.ions_out[index][0], self.ions_out[index][1]
        )
        self.transitionsindata = self.createTable(
            [self.transitionsin_parsed, ["Lines In"]]
        )
        self.transitionsoutdata = self.createTable(
            [self.transitionsout_parsed, ["Lines Out"]]
        )
        self.transitionsindata.add_data(self.transitionsin_count)
        self.transitionsoutdata.add_data(self.transitionsout_count)
        self.transitionsintable2.setModel(self.transitionsindata)
        self.transitionsouttable2.setModel(self.transitionsoutdata)
        self.transitionsintable2.show()
        self.transitionsouttable2.show()
        self.setGeometry(180 + len(self.parent.line_info) * 20, 150, 750, 400)
        self.show()


class LineInteractionTables(QtWidgets.QWidget):
    """Widget to hold the line interaction tables used by
    LineInfo which in turn is used by spectrum widget.
    """

    def __init__(
        self,
        line_interaction_analysis,
        atom_data,
        lines_data,
        description,
        tablecreator,
    ):
        """Create the widget and set data."""
        super(LineInteractionTables, self).__init__()
        self.createTable = tablecreator
        self.text_description = QtWidgets.QLabel(str(description))
        self.species_table = QtWidgets.QTableView()
        self.transitions_table = QtWidgets.QTableView()
        self.layout = QtWidgets.QHBoxLayout()
        self.line_interaction_analysis = line_interaction_analysis
        self.atom_data = atom_data
        self.lines_data = lines_data.reset_index().set_index("line_id")
        line_interaction_species_group = (
            line_interaction_analysis.last_line_in.groupby(
                ["atomic_number", "ion_number"]
            )
        )
        self.species_selected = sorted(
            line_interaction_species_group.groups.keys()
        )
        species_symbols = [
            tardis.util.base.species_tuple_to_string(item)
            for item in self.species_selected
        ]
        species_table_model = self.createTable([species_symbols, ["Species"]])
        species_abundances = (
            (
                line_interaction_species_group.wavelength.count().astype(float)
                / line_interaction_analysis.last_line_in.wavelength.count()
            )
            .astype(float)
            .tolist()
        )
        species_abundances = list(map(float, species_abundances))
        species_table_model.add_data(species_abundances)
        self.species_table.setModel(species_table_model)

        line_interaction_species_group.wavelength.count()
        self.layout.addWidget(self.text_description)
        self.layout.addWidget(self.species_table)
        self.sectionClicked = QtCore.Signal(int)
        self.species_table.verticalHeader().sectionClicked.connect(
            self.on_species_clicked
        )
        self.layout.addWidget(self.transitions_table)

        self.setLayout(self.layout)
        self.show()

    def on_species_clicked(self, index):
        """Needs docstring"""
        current_species = self.species_selected[index]
        last_line_in = self.line_interaction_analysis.last_line_in
        last_line_out = self.line_interaction_analysis.last_line_out

        current_last_line_in = last_line_in.xs(
            key=(current_species[0], current_species[1]),
            level=["atomic_number", "ion_number"],
            drop_level=False,
        ).reset_index()
        current_last_line_out = last_line_out.xs(
            key=(current_species[0], current_species[1]),
            level=["atomic_number", "ion_number"],
            drop_level=False,
        ).reset_index()

        current_last_line_in["line_id_out"] = current_last_line_out.line_id

        last_line_in_string = []
        last_line_count = []
        grouped_line_interactions = current_last_line_in.groupby(
            ["line_id", "line_id_out"]
        )
        exc_deexc_string = "exc. %d-%d (%.2f A) de-exc. %d-%d (%.2f A)"

        for (
            line_id,
            row,
        ) in grouped_line_interactions.wavelength.count().iteritems():
            current_line_in = self.lines_data.loc[line_id[0]]
            current_line_out = self.lines_data.loc[line_id[1]]
            last_line_in_string.append(
                exc_deexc_string
                % (
                    current_line_in["level_number_lower"],
                    current_line_in["level_number_upper"],
                    current_line_in["wavelength"],
                    current_line_out["level_number_upper"],
                    current_line_out["level_number_lower"],
                    current_line_out["wavelength"],
                )
            )
            last_line_count.append(int(row))

        last_line_in_model = self.createTable(
            [
                last_line_in_string,
                [f"Num. pkts {current_last_line_in.wavelength.count()}"],
            ]
        )
        last_line_in_model.add_data(last_line_count)
        self.transitions_table.setModel(last_line_in_model)


class Tardis(QtWidgets.QMainWindow):
    """Create the top level window for the GUI and wait for call to
    display data.
    """

    def __init__(self, tablemodel, config=None, atom_data=None, parent=None):
        """Create the top level window and all widgets it contains.

        When called with no arguments it initializes the GUI in passive
        mode. When a yaml config file and atom data are provided the
        GUI starts in the active mode.

        Parameters
        ----------
        parent : None
            Set to None by default and shouldn't be changed unless
            you are developing something new.
        config : string
            yaml file with configuration information for TARDIS.
        atom_data : string
            hdf file that has the atom data.

        Raises
        ------
        TemporarilyUnavaliable
            Raised when an attempt is made to start the active mode.
            This will be removed when active mode is developed.
        """

        # assumes that qt has already been initialized by starting IPython
        # with the flag "--pylab=qt"gut
        # app = QtCore.QCoreApplication.instance()
        # if app is None:
        #     app = QtGui.QApplication([])
        # try:
        #     from IPython.lib.guisupport import start_event_loop_qt5
        #     start_event_loop_qt5(app)
        # except ImportError:
        #     app.exec_()

        QtWidgets.QMainWindow.__init__(self, parent)

        # path to icons folder
        self.path = os.path.join(tardis.__path__[0], "gui", "images")

        # Check if configuration file was provided
        self.mode = "passive"
        if config is not None:
            self.mode = "active"

        # Statusbar
        statusbr = self.statusBar()
        lblstr = '<font color="red"><b>Calculation did not converge</b></font>'
        self.successLabel = QtWidgets.QLabel(lblstr)
        self.successLabel.setFrameStyle(
            QtWidgets.QFrame.StyledPanel | QtWidgets.QFrame.Sunken
        )
        statusbr.addPermanentWidget(self.successLabel)
        self.modeLabel = QtWidgets.QLabel("Passive mode")
        statusbr.addPermanentWidget(self.modeLabel)
        statusbr.showMessage(self.mode, 5000)
        statusbr.showMessage("Ready", 5000)

        # Actions
        quitAction = QtWidgets.QAction("&Quit", self)
        quitAction.setIcon(
            QtGui.QIcon(os.path.join(self.path, "closeicon.png"))
        )
        quitAction.triggered.connect(self.close)

        self.viewMdv = QtWidgets.QAction("View &Model", self)
        self.viewMdv.setIcon(
            QtGui.QIcon(os.path.join(self.path, "mdvswitch.png"))
        )
        self.viewMdv.setCheckable(True)
        self.viewMdv.setChecked(True)
        self.viewMdv.setEnabled(False)
        self.viewMdv.triggered.connect(self.switch_to_mdv)

        self.viewForm = QtWidgets.QAction("&Edit Model", self)
        self.viewForm.setIcon(
            QtGui.QIcon(os.path.join(self.path, "formswitch.png"))
        )
        self.viewForm.setCheckable(True)
        self.viewForm.setEnabled(False)
        self.viewForm.triggered.connect(self.switch_to_form)

        # Menubar
        self.fileMenu = self.menuBar().addMenu("&File")
        self.fileMenu.addAction(quitAction)
        self.viewMenu = self.menuBar().addMenu("&View")
        self.viewMenu.addAction(self.viewMdv)
        self.viewMenu.addAction(self.viewForm)
        self.helpMenu = self.menuBar().addMenu("&Help")

        # Toolbar
        fileToolbar = self.addToolBar("File")
        fileToolbar.setObjectName("FileToolBar")
        fileToolbar.addAction(quitAction)

        viewToolbar = self.addToolBar("View")
        viewToolbar.setObjectName("ViewToolBar")
        viewToolbar.addAction(self.viewMdv)
        viewToolbar.addAction(self.viewForm)

        # Central Widget
        self.stackedWidget = QtWidgets.QStackedWidget()
        self.mdv = ModelViewer(tablemodel)
        self.stackedWidget.addWidget(self.mdv)

        # In case of active mode
        if self.mode == "active":
            # Disabled currently
            # self.formWidget = ConfigEditor(config)
            # #scrollarea
            # scrollarea = QtGui.QScrollArea()
            # scrollarea.setWidget(self.formWidget)
            # self.stackedWidget.addWidget(scrollarea)
            # self.viewForm.setEnabled(True)
            # self.viewMdv.setEnabled(True)
            # model = run_tardis(config, atom_data)
            # self.show_model(model)
            raise TemporarilyUnavaliable(
                "The active mode is under"
                "development. Please use the passive mode for now."
            )

        self.setCentralWidget(self.stackedWidget)

    def show_model(self, model=None):
        """Set the provided model into the GUI and show the main window.

        Parameters
        ----------
        model : TARDIS model object
            A keyword argument that takes the tardis model object.
        """
        if model:
            self.mdv.change_model(model)
        if model.converged:
            self.successLabel.setText('<font color="green">converged</font>')
        if self.mode == "active":
            self.modeLabel.setText("Active Mode")

        self.mdv.fill_output_label()
        self.mdv.tableview.setModel(self.mdv.tablemodel)
        self.mdv.plot_model()
        self.mdv.plot_spectrum()
        self.showMaximized()

    def switch_to_mdv(self):
        """Switch the cental stacked widget to show the modelviewer."""
        self.stackedWidget.setCurrentIndex(0)
        self.viewForm.setChecked(False)

    def switch_to_form(self):
        """Switch the cental stacked widget to show the ConfigEditor."""
        self.stackedWidget.setCurrentIndex(1)
        self.viewMdv.setChecked(False)


class TemporarilyUnavaliable(Exception):
    """Exception raised when creation of active mode of tardis is attempted."""

    def __init__(self, value):
        self.value = value

    def __str__(self):
        return repr(self.value)
