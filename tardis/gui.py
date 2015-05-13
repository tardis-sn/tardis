import numpy as np
import matplotlib
matplotlib.style.use('fivethirtyeight')
matplotlib.rcParams['font.family']='serif'
matplotlib.rcParams['font.size']=10.0
matplotlib.rcParams['lines.linewidth']=1.0
matplotlib.rcParams['axes.formatter.use_mathtext']=True
matplotlib.rcParams['axes.edgecolor']=matplotlib.rcParams['grid.color']
matplotlib.rcParams['axes.linewidth']=matplotlib.rcParams['grid.linewidth']
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

from tardis import resource_rc
import yaml

import exceptions
#Run this command before importing resource_rc
#pyside-rcc resources.qrc -o resource_rc.py

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

#The main TARDIS window
class Tardis(QtGui.QMainWindow):

    def __init__(self, parent=None, config=None, atom_data=None):
        # assumes that qt has already been initialized by starting IPython with the flag "--pylab=qt"
        app = QtCore.QCoreApplication.instance()
        if app is None:
            app = QtGui.QApplication([])
        try:
            from IPython.lib.guisupport import start_event_loop_qt4
            start_event_loop_qt4(app)
        except ImportError:
            app.exec_()

        super(Tardis, self).__init__(parent)

        #Check if configuration file was provided
        self.mode = 'passive'
        if config is not None:
            from tardis import run_tardis
            self.mode = 'active'

        #Statusbar
        statusbr = self.statusBar()
        self.successLabel = QtGui.QLabel('<font color="red"><b>Calculation did not converge</b></font>')
        self.successLabel.setFrameStyle(QtGui.QFrame.StyledPanel | QtGui.QFrame.Sunken)
        statusbr.addPermanentWidget(self.successLabel)
        self.modeLabel = QtGui.QLabel('Passive mode')
        statusbr.addPermanentWidget(self.modeLabel)
        statusbr.showMessage(self.mode, 5000)
        statusbr.showMessage("Ready", 5000) 

        #Actions
        quitAction = QtGui.QAction("&Quit", self)
        quitAction.setIcon(QtGui.QIcon(":/closeicon.png"))
        quitAction.triggered.connect(self.close)
        
        self.viewMdv = QtGui.QAction("View &Model", self)
        self.viewMdv.setIcon(QtGui.QIcon(":/mdvswitch.png"))
        self.viewMdv.setCheckable(True)
        self.viewMdv.setChecked(True)
        self.viewMdv.setEnabled(False)
        self.viewMdv.triggered.connect(self.switchToMdv)
        
        self.viewForm = QtGui.QAction("&Edit Model", self)
        self.viewForm.setIcon(QtGui.QIcon(":/formswitch.png"))
        self.viewForm.setCheckable(True)
        self.viewForm.setEnabled(False)
        self.viewForm.triggered.connect(self.switchToForm)

        #Menubar
        self.fileMenu = self.menuBar().addMenu("&File")
        self.fileMenu.addAction(quitAction)
        self.viewMenu = self.menuBar().addMenu("&View")
        self.viewMenu.addAction(self.viewMdv)
        self.viewMenu.addAction(self.viewForm)
        self.helpMenu = self.menuBar().addMenu("&Help")

        #Toolbar
        fileToolbar = self.addToolBar("File")
        fileToolbar.setObjectName("FileToolBar")  
        fileToolbar.addAction(quitAction)

        viewToolbar = self.addToolBar("View")
        viewToolbar.setObjectName("ViewToolBar")
        viewToolbar.addAction(self.viewMdv)
        viewToolbar.addAction(self.viewForm)

        #Central Widget
        self.stackedWidget = QtGui.QStackedWidget()
        self.mdv = ModelViewer() 
        self.stackedWidget.addWidget(self.mdv)
        
        #In case of active mode
        if self.mode == 'active':
            self.formWidget = ConfigEditor(config)
            #scrollarea
            scrollarea = QtGui.QScrollArea()
            scrollarea.setWidget(self.formWidget)
            self.stackedWidget.addWidget(scrollarea)
            self.viewForm.setEnabled(True)
            self.viewMdv.setEnabled(True)
            model = run_tardis(config, atom_data)
            self.show_model(model)

        self.setCentralWidget(self.stackedWidget)

    def show_model(self, model=None):
        if model:
            self.mdv.change_model(model)
        if model.converged:
            self.successLabel.setText('<font color="green">converged</font>')
        if self.mode == 'active':
            self.modeLabel.setText('Active Mode')

        self.mdv.fillOutputLabel()
        self.mdv.tableview.setModel(self.mdv.tablemodel)
        self.mdv.plot_model()
        self.mdv.plot_spectrum()
        self.showMaximized()

    #Note to self: Add slot decorator decorator. And try to get all the
    #view changer slots into a single function
    def switchToMdv(self):
        self.stackedWidget.setCurrentIndex(0)
        self.viewForm.setChecked(False)
    def switchToForm(self):
        self.stackedWidget.setCurrentIndex(1)
        self.viewMdv.setChecked(False)


class ConfigEditor(QtGui.QWidget):
    
    def __init__(self, yamlconfigfile, parent=None):

        super(ConfigEditor, self).__init__(parent)
        
        #Configurations from the input and template    
        configDict = yaml.load(open(yamlconfigfile))
        templatedictionary ={'tardis_config_version':[True, 'v1.0'],
                                'supernova':{ 'luminosity_requested':[True, '1 solLum'],
                                              'time_explosion':[True, None],
                                              'distance':[False, None],
                                              'luminosity_wavelength_start':[False, '0 angstrom'],
                                              'luminosity_wavelength_end':[False, 'inf angstrom'],
                                            },
                                'atom_data':[True,'File Browser'],
                                'plasma':{ 'initial_t_inner':[False, '-1K'],
                                           'initial_t_rad':[False,'10000K'],
                                           'disable_electron_scattering':[False, False],
                                           'ionization':[True, None],
                                           'excitation':[True, None],
                                           'radiative_rates_type':[True, None],
                                           'line_interaction_type':[True, None],
                                           'w_epsilon':[False, 1e-10],
                                           'delta_treatment':[False, None],
                                           'nlte':{ 'species':[False, []],
                                                    'coronal_approximation':[False, False],
                                                    'classical_nebular':[False, False]
                                                  }
                                          },
                                'model':{ 'structure':{'type':[True, ['file|_:_|filename|_:_|filetype|_:_|v_inner_boundary|_:_|v_outer_boundary', 'specific|_:_|velocity|_:_|density']],
                                          'filename':[True, None],
                                          'filetype':[True, None],
                                          'v_inner_boundary':[False, '0 km/s'],
                                          'v_outer_boundary':[False, 'inf km/s'],
                                          'velocity':[True, None],
                                          'density':{ 'type':[True, ['branch85_w7|_:_|w7_time_0|_:_|w7_time_0|_:_|w7_time_0','exponential|_:_|time_0|_:_|rho_0|_:_|v_0','power_law|_:_|time_0|_:_|rho_0|_:_|v_0|_:_|exponent', 'uniform|_:_|value']],
                                                      'w7_time_0':[False, '0.000231481 day'],
                                                      'w7_rho_0':[False, '3e29 g/cm^3'],
                                                      'w7_v_0': [False, '1 km/s'],
                                                      'time_0':[True, None],
                                                      'rho_0':[True, None],
                                                      'v_0': [True, None], 
                                                      'exponent': [True, None],
                                                      'value':[True, None] 
                                                    }
                                                      },
                                          'abundances':{ 'type':[True, ['file|_:_|filetype|_:_|filename', 'uniform']],
                                                         'filename':[True, None],
                                                         'filetype':[False, None]
                                                        }
                                        },
                                'montecarlo':{'seed':[False, 23111963],
                                              'no_of_packets':[True, None],
                                              'iterations':[True, None],
                                              'black_body_sampling':{
                                                                        'start': '1 angstrom',
                                                                        'stop': '1000000 angstrom',
                                                                        'num': '1.e+6',
                                                                    },
                                              'last_no_of_packets':[False, -1],
                                              'no_of_virtual_packets':[False, 0],
                                              'enable_reflective_inner_boundary':[False, False],
                                              'inner_boundary_albedo':[False, 0.0],
                                              'convergence_strategy':{ 'type':[True, ['damped|_:_|damping_constant|_:_|t_inner|_:_|t_rad|_:_|w|_:_|lock_t_inner_cycles|_:_|t_inner_update_exponent', 'specific|_:_|threshold|_:_|fraction|_:_|hold_iterations|_:_|t_inner|_:_|t_rad|_:_|w|_:_|lock_t_inner_cycles|_:_|damping_constant|_:_|t_inner_update_exponent']],
                                                                       't_inner_update_exponent':[False, -0.5],
                                                                       'lock_t_inner_cycles':[False, 1],
                                                                       'hold_iterations':[True, 3],
                                                                       'fraction':[True, 0.8],
                                                                       'damping_constant':[False, 0.5],
                                                                       'threshold':[True, None],
                                                                       't_inner':{ 'damping_constant':[False, 0.5],
                                                                                   'threshold': [False, None]
                                                                                 },
                                                                       't_rad':{'damping_constant':[False, 0.5],
                                                                                'threshold':[True, None]
                                                                                },
                                                                        'w':{'damping_constant': [False, 0.5],
                                                                             'threshold': [True, None]
                                                                             }
                                                                    }
                                              },
                                'spectrum':[True, None]
                                }
        self.matchDicts(configDict, templatedictionary)

        self.layout = QtGui.QVBoxLayout()

        #Make tree
        self.trmodel = TreeModel(templatedictionary)
        self.colView = QtGui.QColumnView()
        self.colView.setModel(self.trmodel)
        self.colView.setFixedWidth(256*5) #Five columns of width 256 each can be visible at once
        self.colView.setItemDelegate(TreeDelegate(self))
        self.layout.addWidget(self.colView)

        #Recalculate button
        button = QtGui.QPushButton('Recalculate')
        button.setFixedWidth(90)
        self.layout.addWidget(button)
        button.clicked.connect(self.recalculate) 

        #Finally put them all in
        self.setLayout(self.layout)

    def matchDicts(self, dict1, dict2): #dict1<=dict2
        for key in dict1:
            if key in dict2:
                if isinstance(dict2[key], dict):
                    self.matchDicts(dict1[key], dict2[key])

                elif isinstance(dict2[key], list):
                    if isinstance(dict2[key][1], list):
                        
                        #options = dict2[key][1] #This is passed by reference. So copy the list manually.
                        options = [dict2[key][1][i] for i in range(len(dict2[key][1]))]

                        for i in range(len(options)):
                            options[i] = options[i].split('|_:_|')[0]

                        optionselected = dict1[key]

                        if optionselected in options:
                            indexofselected = options.index(optionselected)
                            temp = dict2[key][1][0]
                            
                            dict2[key][1][0] = dict2[key][1][indexofselected]
                            dict2[key][1][indexofselected] = temp
                            
                            
                        else:
                            print 'The selected and available options'
                            print optionselected
                            print options
                            raise exceptions.IOError("An invalid option was provided in the input file: ")
                                #+str(key)+'='+optionselected)+" options="+str(options))

                else:
                    dict2[key] = dict1[key]
            else:
                toappend = [False, dict1[key]]
                dict2[key] = toappend


    def recalculate(self):
        pass

class Node(object):
    """Object that is an item of the tree"""

    def __init__(self, data, parent=None):
        self.parent = parent
        self.children = []
        self.data = data
        self.siblings = {} #For 'type' fields. Will store the nodes to 
                           #enable disable on selection

    def appendChild(self, child):
        self.children.append(child)
        child.parent = self

    def getChild(self, i):
        if i < self.numChildren():
            return self.children[i]
        else:
            return None

    def numChildren(self):
        return len(self.children)

    def numColumns(self):
        return len(self.data)

    def getData(self, i):
        try:
            return self.data[i]
        except IndexError:
            return None

    def getParent(self):
        return self.parent

    def getIndexOfSelf(self):
        if self.parent:
            return self.parent.children.index(self)
        else:
            return 0

    #--------------------------------------------------------------------------------------------------------------
    #For editable model
    def insertChildren(self, position, count, columns):
        if (position < 0) or (position > self.numChildren()):
            return False

        for row in range(count):
            data = [None for v in range(columns)]
            item = Node(data, self)
            self.children.insert(position, item)

        return True 

    def removeChildren(self, position, count):
        if (position < 0) or ((position + count) > self.numChildren()):
            return False

        for row in range(count):
            self.children.pop(position)

        return True

    def insertColumns(self, position, columns):
        if (position < 0) or (position > self.numColumns()):
            return False

        for column in range(columns):
            self.data.insert(position, None)

        for child in self.children:
            child.insertColumns(position, columns)

        return True

    def removeColumns(self, position, columns):
        if (position < 0) or ((position + columns) > self.numColumns()):
            return False

        for column in range(columns):
            self.itemData.pop(position)

        for child in self.children:
            child.removeColumns(position, columns)

        return True

    def setData(self, column, value):
        if column < 0 or column >= self.numColumns():
            return False

        self.data[column] = value

        return True

#Note to self: For columnview headerdata seems to be unused. Try removing.
class TreeModel(QtCore.QAbstractItemModel):
    def __init__(self, dictionary, parent=None):
        super(TreeModel, self).__init__(parent)

        self.root = Node(["column A"])
        self.disabledNodes = []
        self.typenodes = []
        self.dictToTree(dictionary, self.root)

    def dictToTree(self, dictionary, root):
        #Construct tree with all nodes
        self.treeFromNode(dictionary, root)

        #Append siblings to type nodes
        for node in self.typenodes: #For every type node
            parent = node.getParent()
            sibsdict = {}
            for i in range(parent.numChildren()):
                sibsdict[parent.getChild(i).getData(0)] = parent.getChild(i)

            typesleaf = node.getChild(0)
            for i in range(typesleaf.numColumns()):
                sibstrings = typesleaf.getData(i).split('|_:_|')
            
                typesleaf.setData(i, sibstrings[0])
                sibslist = []
                for j in range(1, len(sibstrings)):
                    if sibstrings[j] in sibsdict:
                        sibslist.append(sibsdict[sibstrings[j]])

                typesleaf.siblings[sibstrings[0]] = sibslist
            
            #Then append siblings of current selection for all type nodes to
            #disabled nodes
            
            for i in range(1,typesleaf.numColumns()):
                key = typesleaf.getData(i)
                for nd in typesleaf.siblings[key]:
                    self.disabledNodes.append(nd)


    def treeFromNode(self, dictionary, root):
        for key in dictionary:
            child = Node([key])
            root.appendChild(child)
            if isinstance(dictionary[key], dict):
                self.treeFromNode(dictionary[key], child)
            elif isinstance(dictionary[key], list):
                #if dictionary[key][0]:
                #    data = child.getData(0)
                #    data = data + '*'
                #    child.setData(0, data)

                if isinstance(dictionary[key][1], list):
                    leaf = Node(dictionary[key][1])    
                else:
                    leaf = Node([dictionary[key][1]])

                child.appendChild(leaf)
                if key == 'type':
                    self.typenodes.append(child)
            #else: 
            #    leaf = Node([dictionary[key]])
            #    child.appendChild(leaf)

    def dictFromNode(self, node): 
        children = [node.getChild(i) for i in range(node.numChildren())]
        if len(children) > 1:
            dictionary = {}
            for nd in children:
                if nd in self.disabledNodes:
                    pass
                else:
                    dictionary[nd.getData(0)] = self.dictFromNode(nd)
            return dictionary
        elif len(children)==1:
            return children[0].getData(0)

    def columnCount(self, index):
        if index.isValid():
            return index.internalPointer().numColumns()
        else:
            return self.root.numColumns()

    def data(self, index, role):
        if not index.isValid():
            return None

        if role != QtCore.Qt.DisplayRole:
            return None

        item = index.internalPointer()

        return item.getData(index.column())

    def flags(self, index):
        if not index.isValid():
            return QtCore.Qt.NoItemFlags

        node = index.internalPointer()
        if (node.getParent() in self.disabledNodes) or (node in self.disabledNodes):
            return QtCore.Qt.NoItemFlags

        if node.numChildren()==0:
            return QtCore.Qt.ItemIsEditable | QtCore.Qt.ItemIsEnabled | QtCore.Qt.ItemIsSelectable

        return QtCore.Qt.ItemIsEnabled | QtCore.Qt.ItemIsSelectable        

    def getItem(self, index):
        if index.isValid():
            item = index.internalPointer()
            if item:
                return item

        return self.root

    def headerData(self, section, orientation, role):
        if orientation == QtCore.Qt.Horizontal and role == QtCore.Qt.DisplayRole:
            return self.root.getData(section)

        return None

    def index(self, row, column, parent=QtCore.QModelIndex()):
        if parent.isValid() and parent.column() != 0:
            return QtCore.QModelIndex()

        parentItem = self.getItem(parent)
        childItem = parentItem.getChild(row)
        if childItem:
            return self.createIndex(row, column, childItem)
        else:
            return QtCore.QModelIndex()

    def insertColumns(self, position, columns, parent=QtCore.QModelIndex()):
        self.beginInsertColumns(parent, position, position + columns - 1)
        success = self.root.insertColumns(position, columns)
        self.endInsertColumns()

        return success

    def insertRows(self, position, rows, parent=QtCore.QModelIndex()):
        parentItem = self.getItem(parent)
        self.beginInsertRows(parent, position, position + rows - 1)
        success = parentItem.insertChildren(position, rows,
                self.rootItem.columnCount())
        self.endInsertRows()

        return success

    def removeColumns(self, position, columns, parent=QtCore.QModelIndex()):
        self.beginRemoveColumns(parent, position, position + columns - 1)
        success = self.root.removeColumns(position, columns)
        self.endRemoveColumns()

        if self.root.numColumns() == 0:
            self.removeRows(0, rowCount())

        return success

    def removeRows(self, position, rows, parent=QtCore.QModelIndex()):
        parentItem = self.getItem(parent)

        self.beginRemoveRows(parent, position, position + rows - 1)
        success = parentItem.removeChildren(position, rows)
        self.endRemoveRows()

        return success

    def parent(self, index):
        if not index.isValid():
            return QtCore.QModelIndex()

        childItem = index.internalPointer()
        parentItem = childItem.getParent()

        if parentItem == self.root:
            return QtCore.QModelIndex()

        return self.createIndex(parentItem.getIndexOfSelf(), 0, parentItem)

    def rowCount(self, parent=QtCore.QModelIndex()):
        parentItem = self.getItem(parent)

        return parentItem.numChildren()

    def setData(self, index, value, role=QtCore.Qt.EditRole):
        if role != QtCore.Qt.EditRole:
            return False

        item = self.getItem(index)
        result = item.setData(index.column(), value)

        if result:
            self.dataChanged.emit(index, index)
            #print self.dictFromNode(self.root)
        return result

    def setHeaderData(self, section, orientation, value, role=QtCore.Qt.EditRole):
        if role != QtCore.Qt.EditRole or orientation != QtCore.Qt.Horizontal:
            return False

        result = self.root.setData(section, value)
        if result:
            self.headerDataChanged.emit(orientation, section, section)

        return result

class TreeDelegate(QtGui.QStyledItemDelegate):
    def __init__(self, parent=None):
        super(TreeDelegate, self).__init__(parent)

    def createEditor(self, parent, option, index):
        node = index.internalPointer()
        if node.numColumns()>1:
            combobox = QtGui.QComboBox(parent)
            combobox.addItems([node.getData(i) for i in range(node.numColumns())])
            combobox.setEditable(False)
            return combobox
        else:
            #return QtGui.QStyledItemDelegate.createEditor(self, parent, option, index)
            editor =  QtGui.QLineEdit(parent)
            editor.setText(str(node.getData(0)))
            editor.returnPressed.connect(self.closeAndCommit)
            return editor

    def closeAndCommit(self):
        editor = self.sender()
        if isinstance(editor, QtGui.QLineEdit):
            self.commitData.emit(editor)
            self.closeEditor.emit(editor, QtGui.QAbstractItemDelegate.NoHint)

    def setModelData(self, editor, model, index):
        node = index.internalPointer()

        if node.numColumns() > 1 and node.getParent().getData(0) != 'type':
            selectedIndex = editor.currentIndex()
            firstItem = node.getData(0)
            node.setData(0, str(editor.currentText()))
            node.setData(selectedIndex, str(firstItem))

        elif node.numColumns() > 1 and node.getParent().getData(0) == 'type':
            selectedIndex = editor.currentIndex()
            firstItem = node.getData(0)
            node.setData(0, str(editor.currentText()))
            node.setData(selectedIndex, str(firstItem))

            itemsToDisable = node.siblings[firstItem]
            itemsToEnable = node.siblings[str(editor.currentText())]

            for nd in itemsToDisable:
                model.disabledNodes.append(nd)

            for nd in itemsToEnable:
                if nd in model.disabledNodes:
                    model.disabledNodes.remove(nd) 

        elif isinstance(editor, QtGui.QLineEdit): 
            node.setData(0, str(editor.text()))
        else:
            QtGui.QStyledItemDelegate.setModelData(self, editor, model, index)
            
        print model.dictFromNode(model.root) 
        f = open('dictester.dat','w')
        f.write(yaml.dump(model.dictFromNode(model.root)))
        f.close()

class ModelViewer(QtGui.QWidget):

    def __init__(self, parent=None):

        super(ModelViewer, self).__init__(parent)
        
        #Data structures
        self.model = None
        self.shell_info = {}
        self.line_info = []

        #Shells widget
        self.shellWidget = self.makeShellWidget()
        
        #Spectrum widget
        self.spectrumWidget = self.makeSpectrumWidget()

        #Plot tab widget
        self.plotTabWidget = QtGui.QTabWidget()
        self.plotTabWidget.addTab(self.shellWidget,"&Shells")
        self.plotTabWidget.addTab(self.spectrumWidget, "S&pectrum")

        #Table widget
        self.tablemodel = SimpleTableModel([['Shell: '], ["Rad. temp", "Ws"]], (1, 0))
        self.tableview = QtGui.QTableView()
        self.tableview.setMinimumWidth(200)
        self.tableview.connect(self.tableview.verticalHeader(), QtCore.SIGNAL('sectionClicked(int)'), self.graph.highlight_shell)
        self.tableview.connect(self.tableview.verticalHeader(), QtCore.SIGNAL('sectionDoubleClicked(int)'),
                               self.on_header_double_clicked)

        #Label for text output
        self.outputLabel = QtGui.QLabel()
        self.outputLabel.setFrameStyle(QtGui.QFrame.StyledPanel | QtGui.QFrame.Sunken)
        self.outputLabel.setStyleSheet("QLabel{background-color:white;}")

        #Group boxes
        graphsBox = QtGui.QGroupBox("Visualized results")
        textsBox = QtGui.QGroupBox("Model parameters")
        tableBox = QtGui.QGroupBox("Tabulated results")

        #For textbox
        textlayout = QtGui.QHBoxLayout()
        textlayout.addWidget(self.outputLabel)

        tableslayout = QtGui.QVBoxLayout()
        tableslayout.addWidget(self.tableview)
        tableBox.setLayout(tableslayout)

        visualayout = QtGui.QVBoxLayout()
        visualayout.addWidget(self.plotTabWidget)
        graphsBox.setLayout(visualayout)

        self.layout = QtGui.QHBoxLayout()
        self.layout.addWidget(graphsBox)
        textntablelayout = QtGui.QVBoxLayout()
        textsBox.setLayout(textlayout)
        textntablelayout.addWidget(textsBox)
        textntablelayout.addWidget(tableBox)

        self.layout.addLayout(textntablelayout)                
        self.setLayout(self.layout)

    def fillOutputLabel(self):
        labeltext = 'Iterations requested: {} <br/> Iterations executed:  {}<br/>\
                     Model converged     : {} <br/> Simulation Time    :  {} s <br/>\
                     Inner Temperature   : {} K <br/> Number of packets  :  {}<br/>\
                     Inner Luminosity    : {}'\
                     .format(self.model.iterations_max_requested, self.model.iterations_executed,\
                        '<font color="green"><b>True</b></font>' if self.model.converged else \
                        '<font color="red"><b>False</b></font>', self.model.time_of_simulation.value,\
                        self.model.t_inner.value, self.model.current_no_of_packets,\
                        self.model.luminosity_inner)
        self.outputLabel.setText(labeltext)

    def makeShellWidget(self):
        #Widgets for plot of shells
        self.graph = MatplotlibWidget(self, 'model')
        self.graph_label = QtGui.QLabel('Select Property:')
        self.graph_button = QtGui.QToolButton()
        self.graph_button.setText('Rad. temp')
        self.graph_button.setPopupMode(QtGui.QToolButton.MenuButtonPopup)
        self.graph_button.setMenu(QtGui.QMenu(self.graph_button))
        self.graph_button.menu().addAction('Rad. temp').triggered.connect(self.change_graph_to_t_rads)
        self.graph_button.menu().addAction('Ws').triggered.connect(self.change_graph_to_ws)

        #Layouts: bottom up
        self.graph_subsublayout = QtGui.QHBoxLayout()
        self.graph_subsublayout.addWidget(self.graph_label)
        self.graph_subsublayout.addWidget(self.graph_button)

        self.graph_sublayout = QtGui.QVBoxLayout()
        self.graph_sublayout.addLayout(self.graph_subsublayout)
        self.graph_sublayout.addWidget(self.graph)

        containerWidget = QtGui.QWidget()
        containerWidget.setLayout(self.graph_sublayout)
        return containerWidget

    def makeSpectrumWidget(self):
        self.spectrum = MatplotlibWidget(self)
        self.spectrum_label = QtGui.QLabel('Select Spectrum:')
        self.spectrum_button = QtGui.QToolButton()
        self.spectrum_button.setText('spec_flux_angstrom')
        self.spectrum_button.setPopupMode(QtGui.QToolButton.MenuButtonPopup)
        self.spectrum_button.setMenu(QtGui.QMenu(self.spectrum_button))
        self.spectrum_button.menu().addAction('spec_flux_angstrom').triggered.connect(self.change_spectrum_to_spec_flux_angstrom)
        self.spectrum_button.menu().addAction('spec_virtual_flux_angstrom').triggered.connect(self.change_spectrum_to_spec_virtual_flux_angstrom)
        self.spectrum_span_button = QtGui.QPushButton('Show Wavelength Range')
        self.spectrum_span_button.clicked.connect(self.spectrum.show_span)
        self.spectrum_line_info_button = QtGui.QPushButton('Show Line Info')
        self.spectrum_line_info_button.hide()
        self.spectrum_line_info_button.clicked.connect(self.spectrum.show_line_info)

        self.spectrum_subsublayout = QtGui.QHBoxLayout()
        self.spectrum_subsublayout.addWidget(self.spectrum_span_button)
        self.spectrum_subsublayout.addWidget(self.spectrum_label)
        self.spectrum_subsublayout.addWidget(self.spectrum_button)

        self.spectrum_sublayout = QtGui.QVBoxLayout()
        self.spectrum_sublayout.addLayout(self.spectrum_subsublayout)
        self.spectrum_sublayout.addWidget(self.spectrum_line_info_button)
        self.spectrum_sublayout.addWidget(self.spectrum)
        self.spectrum_sublayout.addWidget(self.spectrum.toolbar)

        containerWidget = QtGui.QWidget()
        containerWidget.setLayout(self.spectrum_sublayout)
        return containerWidget


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
        self.graph.ax2.set_xticklabels([])
        self.graph.ax2.set_yticklabels([])
        self.graph.ax2.grid = True

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
        return None

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
        self.figure = Figure()#(frameon=False,facecolor=(1,1,1))
        self.cid = {}
        if fig != 'model':
            self.ax = self.figure.add_subplot(111)
        else:
            self.gs = gridspec.GridSpec(2, 1, height_ratios=[1, 3])
            self.ax1 = self.figure.add_subplot(self.gs[0])
            self.ax2 = self.figure.add_subplot(self.gs[1])#, aspect='equal')
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
