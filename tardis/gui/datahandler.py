import os
from pkg_resources import parse_version

import numpy as np
import matplotlib
import matplotlib.pylab as plt


if os.environ.get("QT_API", None) == "pyqt":
    from PyQt5 import QtGui, QtCore, QtWidgets
elif os.environ.get("QT_API", None) == "pyside":
    from PySide2 import QtGui, QtCore, QtWidgets
else:
    raise ImportError(
        "QT_API was not set! Please exit the IPython console\n"
        " and at the bash prompt use : \n\n export QT_API=pyside \n or\n"
        " export QT_API=pyqt \n\n For more information refer to user guide."
    )
import yaml

from tardis import run_tardis
from tardis.gui.widgets import MatplotlibWidget, ModelViewer, ShellInfo
from tardis.gui.widgets import LineInfo, LineInteractionTables

if parse_version(matplotlib.__version__) >= parse_version("1.4"):
    matplotlib.style.use("fivethirtyeight")
else:
    print("Please upgrade matplotlib to a version >=1.4 for best results!")
matplotlib.rcParams["font.family"] = "serif"
matplotlib.rcParams["font.size"] = 10.0
matplotlib.rcParams["lines.linewidth"] = 1.0
matplotlib.rcParams["axes.formatter.use_mathtext"] = True
matplotlib.rcParams["axes.edgecolor"] = matplotlib.rcParams["grid.color"]
matplotlib.rcParams["axes.linewidth"] = matplotlib.rcParams["grid.linewidth"]


class Node(object):
    """Object that serves as the nodes in the TreeModel.

    Attributes
    ----------
        parent: None/Node
            The parent of the node.
        children: list of Node
            The children of the node.
        data: list of string
            The data stored on the node. Can be a key or a value.
        siblings: dictionary
            A dictionary of nodes that are siblings of this node. The
            keys are the values of the nodes themselves. This is
            used to keep track of which value the user has selected
            if the parent of this node happens to be a key that can
            take values from a list.

    """

    def __init__(self, data, parent=None):
        """Create one node with the data and parent provided.

        Parameters
        ----------
            data: list of string
                The data that is intended to be stored on the node.
            parent: Node
                Another node which is the parent of this node. The
                root node has parent set to None.

        Note
        ----
            A leaf node is a node that is the only child of its parent.
            For this tree this will always be the case. This is because
            the tree stores the key at every node except the leaf node
            where it stores the value for the key. So if in the dictionary
            the value of a key is another dictionary, then it will
            be a node with no leafs. If the key has a value that is a
            value or a list then it will have one child that is a leaf.
            The leaf can have no children. For example-

            'a':val1
            'b':{'c':val2
                 'd':val3
                 'e':{'f':val4
                      'g':val5
                         :val6
                         :val7}} *In this case the key g can take values
                                  val7, val5 and val6 and is currently
                                  set to val5.

            In the tree shown above all quoted values are keys in the
            dictionary and are non-leaf nodes in the tree. All values
            of the form valx are leaf nodes and are not dictionaries
            themselves. If the keys have non-dictionary values then they
            have a leaf attached. And no leaf can have a child.

        """
        self.parent = parent
        self.children = []
        self.data = data
        self.siblings = {}  # For 'type' fields. Will store the nodes to
        # enable disable on selection

    def append_child(self, child):
        """Add a child to this node."""
        self.children.append(child)
        child.parent = self

    def get_child(self, i):
        """Get the ith child of this node.

        No error is raised if the cild requested doesn't exist. A
        None is returned in such cases.

        """
        if i < self.num_children():
            return self.children[i]
        else:
            return None

    def num_children(self):
        """Number of children this node has."""
        return len(self.children)

    def num_columns(self):
        """Returns the number of strings stored in the data attribute."""
        return len(self.data)

    def get_data(self, i):
        """Returns the ith string from the data list.

        No error is raised if the data list index is exceeded. A None is
        returned in such cases.

        """
        try:
            return self.data[i]
        except IndexError:
            return None

    def get_parent(self):
        """Return the parent of this node."""
        return self.parent

    def get_index_of_self(self):
        """Returns the number at which it comes in the list of its
        parent's children. For root the index 0 is returned.

        """
        if self.parent:
            return self.parent.children.index(self)
        else:
            return 0

    def set_data(self, column, value):
        """Set the data for the ith index to the provided value. Returns
        true if the data was set successfully.

        """
        if column < 0 or column >= self.num_columns():
            return False

        self.data[column] = value

        return True


class TreeModel(QtCore.QAbstractItemModel):
    """The class that defines the tree for ConfigEditor.

    Parameters
    ----------
        root: Node
            Root node of the tree.
        disabledNodes: list of Node
            List of leaf nodes that are not editable currently.
        typenodes: list of Node
            List of nodes that correspond to keys that set container
            types. Look at tardis configuration template. These are the
            nodes that have values that can be set from a list.

    """

    def __init__(self, dictionary, parent=None):
        """Create a tree of tardis configuration dictionary.

        Parameters
        ----------
            dictionary: dictionary
                The dictionary that needs to be converted to the tree.
            parent: None
                Used to instantiate the QAbstractItemModel

        """
        QtCore.QAbstractItemModel.__init__(self, parent)

        self.root = Node(["column A"])
        self.disabledNodes = []
        self.typenodes = []
        self.dict_to_tree(dictionary, self.root)

    # mandatory functions for subclasses
    def columnCount(self, index):
        """Return the number of columns in the node pointed to by
        the given model index.

        """
        if index.isValid():
            return index.internalPointer().num_columns()
        else:
            return self.root.num_columns()

    def data(self, index, role):
        """Returns the asked data for the node specified by the modeLabel
        index."""
        if not index.isValid():
            return None

        if role != QtCore.Qt.DisplayRole:
            return None

        item = index.internalPointer()

        return item.get_data(index.column())

    def flags(self, index):
        """Return flags for the items whose model index is provided."""
        if not index.isValid():
            return QtCore.Qt.NoItemFlags

        node = index.internalPointer()
        if (node.get_parent() in self.disabledNodes) or (
            node in self.disabledNodes
        ):
            return QtCore.Qt.NoItemFlags

        if node.num_children() == 0:
            return (
                QtCore.Qt.ItemIsEditable
                | QtCore.Qt.ItemIsEnabled
                | QtCore.Qt.ItemIsSelectable
            )

        return QtCore.Qt.ItemIsEnabled | QtCore.Qt.ItemIsSelectable

    def getItem(self, index):
        """Returns the node to which the model index is pointing. If the
        model index is invalid then the root node is returned.

        """
        if index.isValid():
            item = index.internalPointer()
            if item:
                return item

        return self.root

    def headerData(self, section, orientation, role):
        """Returns header data. This is not used in QColumnView. But will
        be needed for QTreeView.

        """
        if (
            orientation == QtCore.Qt.Horizontal
            and role == QtCore.Qt.DisplayRole
        ):
            return self.root.get_data(section)

        return None

    def index(self, row, column, parent=QtCore.QModelIndex()):
        """Create a model index for the given row and column. For a
        tree model, the row is the set of nodes with the same parents and
        the column indexes the data in the node.

        """
        if parent.isValid() and parent.column() != 0:
            return QtCore.QModelIndex()

        parentItem = self.getItem(parent)
        childItem = parentItem.get_child(row)
        if childItem:
            return self.createIndex(row, column, childItem)
        else:
            return QtCore.QModelIndex()

    def insertColumns(self, position, columns, parent=QtCore.QModelIndex()):
        """Insert columns in the tree model."""
        self.beginInsertColumns(parent, position, position + columns - 1)
        success = self.root.insertColumns(position, columns)
        self.endInsertColumns()

        return success

    def insertRows(self, position, rows, parent=QtCore.QModelIndex()):
        """Insert rows in the tree model."""
        parentItem = self.getItem(parent)
        self.beginInsertRows(parent, position, position + rows - 1)
        success = parentItem.insertChildren(
            position, rows, self.rootItem.columnCount()
        )
        self.endInsertRows()

        return success

    def parent(self, index):
        """Return the parent of the node to which the index points."""
        if not index.isValid():
            return QtCore.QModelIndex()

        childItem = index.internalPointer()
        parentItem = childItem.get_parent()

        if parentItem == self.root:
            return QtCore.QModelIndex()

        return self.createIndex(parentItem.get_index_of_self(), 0, parentItem)

    def rowCount(self, parent=QtCore.QModelIndex()):
        """The number of rows for a given node.

        (The number of rows is just the number of children for a node.)

        """
        parentItem = self.getItem(parent)

        return parentItem.num_hildren()

    def setData(self, index, value, role=QtCore.Qt.EditRole):
        """Set the value as the data at the location pointed by the
        index.

        """
        if role != QtCore.Qt.EditRole:
            return False

        item = self.getItem(index)
        result = item.setData(index.column(), value)

        if result:
            self.dataChanged.emit(index, index)
        return result

    def setHeaderData(
        self, section, orientation, value, role=QtCore.Qt.EditRole
    ):
        """Change header data. Unused in columnview."""
        if role != QtCore.Qt.EditRole or orientation != QtCore.Qt.Horizontal:
            return False

        result = self.root.setData(section, value)
        if result:
            self.headerDataChanged.emit(orientation, section, section)

        return result

    # Custom functions
    def dict_to_tree(self, dictionary, root):
        """Create the tree and append siblings to nodes that need them.

        Parameters
        ----------
            dictionary: dictionary
                The dictionary that is to be converted to the tree.
            root: Node
                The root node of the tree.

        """
        # Construct tree with all nodes
        self.tree_from_node(dictionary, root)

        # Append siblings to type nodes
        for node in self.typenodes:  # For every type node
            parent = node.get_parent()
            sibsdict = {}
            for i in range(parent.num_children()):
                sibsdict[parent.get_child(i).get_data(0)] = parent.get_child(i)

            typesleaf = node.get_child(0)
            for i in range(typesleaf.num_columns()):
                sibstrings = typesleaf.get_data(i).split("|_:_|")

                typesleaf.set_data(i, sibstrings[0])
                sibslist = []
                for j in range(1, len(sibstrings)):
                    if sibstrings[j] in sibsdict:
                        sibslist.append(sibsdict[sibstrings[j]])

                typesleaf.siblings[sibstrings[0]] = sibslist

            # Then append siblings of current selection for all type nodes to
            # disabled nodes
            for i in range(1, typesleaf.num_columns()):
                key = typesleaf.get_data(i)
                for nd in typesleaf.siblings[key]:
                    self.disabledNodes.append(nd)

    def tree_from_node(self, dictionary, root):
        """Convert dictionary to tree. Called by dict_to_tree."""
        for key in dictionary:
            child = Node([key])
            root.append_child(child)
            if isinstance(dictionary[key], dict):
                self.tree_from_node(dictionary[key], child)
            elif isinstance(dictionary[key], list):
                if isinstance(dictionary[key][1], list):
                    leaf = Node(dictionary[key][1])
                else:
                    leaf = Node([dictionary[key][1]])

                child.append_child(leaf)
                if key == "type":
                    self.typenodes.append(child)

    def dict_from_node(self, node):
        """Take a node and convert the whole subtree rooted at it into a
        dictionary.

        """
        children = [node.get_child(i) for i in range(node.num_children())]
        if len(children) > 1:
            dictionary = {}
            for nd in children:
                if nd in self.disabledNodes:
                    pass
                else:
                    dictionary[nd.get_data(0)] = self.dict_from_node(nd)
            return dictionary
        elif len(children) == 1:
            return children[0].get_data(0)


class TreeDelegate(QtWidgets.QStyledItemDelegate):
    """Create a custom delegate to modify the columnview that displays the
    TreeModel.

    """

    def __init__(self, parent=None):
        """Call the constructor of the superclass."""
        QtWidgets.QStyledItemDelegate.__init__(self, parent)

    # Mandatory methods for subclassing
    def createEditor(self, parent, option, index):
        """Create a lineEdit or combobox depending on the type of node."""
        node = index.internalPointer()
        if node.num_columns() > 1:
            combobox = QtGui.QComboBox(parent)
            combobox.addItems(
                [node.get_data(i) for i in range(node.num_columns())]
            )
            combobox.setEditable(False)
            return combobox
        else:
            editor = QtWidgets.QLineEdit(parent)
            editor.setText(str(node.get_data(0)))
            editor.returnPressed.connect(self.close_and_commit)
            return editor

    def setModelData(self, editor, model, index):
        """Called when new data id set in the model. This is where the
        siblings of type nodes are enabled or disabled according to the
        new choice made.

        """
        node = index.internalPointer()

        if node.num_columns() > 1 and node.get_parent().get_data(0) != "type":
            selectedIndex = editor.currentIndex()
            firstItem = node.get_data(0)
            node.setData(0, str(editor.currentText()))
            node.setData(selectedIndex, str(firstItem))

        elif node.num_columns() > 1 and node.get_parent().get_data(0) == "type":
            selectedIndex = editor.currentIndex()
            firstItem = node.get_data(0)
            node.setData(0, str(editor.currentText()))
            node.setData(selectedIndex, str(firstItem))

            itemsToDisable = node.siblings[firstItem]
            itemsToEnable = node.siblings[str(editor.currentText())]

            for nd in itemsToDisable:
                model.disabledNodes.append(nd)

            for nd in itemsToEnable:
                if nd in model.disabledNodes:
                    model.disabledNodes.remove(nd)

        elif isinstance(editor, QtWidgets.QLineEdit):
            node.setData(0, str(editor.text()))
        else:
            QtWidgets.QStyledItemDelegate.setModelData(
                self, editor, model, index
            )

    # Custom methods
    def close_and_commit(self):
        """Saver for the line edits."""
        editor = self.sender()
        if isinstance(editor, QtWidgets.QLineEdit):
            self.commitData.emit(editor)
            self.closeEditor.emit(
                editor, QtWidgets.QAbstractItemDelegate.NoHint
            )


class SimpleTableModel(QtCore.QAbstractTableModel):
    """Create a table data structure for the table widgets."""

    def __init__(
        self,
        headerdata=None,
        iterate_header=(0, 0),
        index_info=None,
        parent=None,
        *args,
    ):
        """Call constructor of the QAbstractTableModel and set parameters
        given by user.
        """
        super(SimpleTableModel, self).__init__(parent, *args)
        self.headerdata = headerdata
        self.arraydata = []
        self.iterate_header = iterate_header
        self.index_info = index_info

    # Implementing methods mandatory for subclassing QAbstractTableModel
    def rowCount(self, parent=QtCore.QModelIndex()):
        """Return number of rows."""
        return len(self.arraydata[0])

    def columnCount(self, parent=QtCore.QModelIndex()):
        """Return number of columns."""
        return len(self.arraydata)

    def headerData(self, section, orientation, role=QtCore.Qt.DisplayRole):
        """Set the header data."""
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
        elif (
            orientation == QtCore.Qt.Horizontal
            and role == QtCore.Qt.DisplayRole
        ):
            if self.iterate_header[1] == 1:
                return self.headerdata[1][0] + str(section + 1)
            elif self.iterate_header[1] == 2:
                if self.index_info:
                    return self.headerdata[1][0] + str(self.index_info[section])
            else:
                return self.headerdata[1][section]
        return None

    def data(self, index, role=QtCore.Qt.DisplayRole):
        """Return data of specified index and role."""
        if not index.isValid():
            return None
        elif role != QtCore.Qt.DisplayRole:
            return None
        return self.arraydata[index.column()][index.row()]

    def setData(self, index, value, role=QtCore.Qt.EditRole):
        """Change the data in the model for specified index and role
        to specified value."""
        if not index.isValid():
            return False
        elif role != QtCore.Qt.EditRole:
            return False
        self.arraydata[index.column()][index.row()] = value

        self.dataChanged = QtCore.Signal(
            QtGui.QModelIndex(), QtGui.QModelIndex()
        )
        self.dataChanged.emit(index, index)
        return True

    # Methods used to inderact with the SimpleTableModel
    def update_table(self):
        """Update table to set all the new data."""
        for r in range(self.rowCount()):
            for c in range(self.columnCount()):
                index = self.createIndex(r, c)
                self.setData(index, self.arraydata[c][r])

    def add_data(self, datain):
        """Add data to the model."""
        self.arraydata.append(datain)
