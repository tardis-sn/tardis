import os

if os.environ.get('QT_API', None)=='pyqt':
    from PyQt4 import QtCore, QtGui
elif os.environ.get('QT_API', None)=='pyside':
    from PySide import QtCore, QtGui
else:
    raise ImportError('QT_API was not set! Please exit the IPython console\n'
        ' and at the bash prompt use : \n\n export QT_API=pyside \n or\n'
        ' export QT_API=pyqt \n\n For more information refer to user guide.')

import sys
from tardis.gui.widgets import MatplotlibWidget, ModelViewer, ShellInfo
from tardis.gui.widgets import LineInfo, LineInteractionTables, Tardis 
from tardis.gui.datahandler import SimpleTableModel

try:
    from IPython.lib.guisupport import get_app_qt4
    app = get_app_qt4()
except ImportError:
    app = QtGui.QApplication([])

def show(model):
    tablemodel = SimpleTableModel

    win = Tardis(tablemodel)
    win.show_model(model)

    try:
        from IPython.lib.guisupport import start_event_loop_qt4, is_event_loop_running_qt4 
        start_event_loop_qt4(app)
        if is_event_loop_running_qt4(app):
            return win
    except ImportError:
        app.exec_()

if __name__=='__main__':
    from tardis import run_tardis
    
    yamlfile = sys.argv[1]
    atomfile = sys.argv[2]
    mdl = run_tardis(yamlfile, atomfile)
    show(mdl) 