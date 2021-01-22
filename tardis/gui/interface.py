import os

if os.environ.get("QT_API", None) == "pyqt":
    from PyQt5 import QtCore, QtWidgets
elif os.environ.get("QT_API", None) == "pyside":
    from PySide2 import QtCore, QtWidgets
else:
    raise ImportError(
        "QT_API was not set! Please exit the IPython console\n"
        " and at the bash prompt use : \n\n export QT_API=pyside \n or\n"
        " export QT_API=pyqt \n\n For more information refer to user guide."
    )
import sys

try:
    from IPython.lib.guisupport import get_app_qt5, start_event_loop_qt5
    from IPython.lib.guisupport import is_event_loop_running_qt5

    importFailed = False
except ImportError:
    importFailed = True

from tardis.gui.widgets import Tardis
from tardis.gui.datahandler import SimpleTableModel
from tardis import run_tardis


def show(model):
    """Take an instance of tardis model and display it.

    If IPython functions were successfully imported then QApplication instance
    is created using get_app_qt4. This will only create the app instance if
    it doesn't already exists. Otherwise it is explicitly started.

    Then the mainwindow is created, the model is attached to it and its
    show method is called.

    Finally the eventloop is started using IPython functions (which start them
    consistently) if they were imported. Otherwise it is started explicitly.

    """
    if importFailed:
        app = QtWidgets.QApplication([])
    else:
        app = get_app_qt5()

    tablemodel = SimpleTableModel
    win = Tardis(tablemodel)
    win.show_model(model)

    if importFailed:
        app.exec_()
    else:
        start_event_loop_qt5(app)

        # If the IPython console is being used, this will evaluate to true.
        # In that case the window created will be garbage collected unless a
        # reference to it is maintained after this function exits. So the win is
        # returned.
        if is_event_loop_running_qt5(app):
            return win


if __name__ == "__main__":
    """When this module is executed as script, take arguments, calculate model
    and call the show function.

    """
    yamlfile = sys.argv[1]
    atomfile = sys.argv[2]
    mdl = run_tardis(yamlfile, atomfile)
    show(mdl)
