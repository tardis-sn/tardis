************************
Graphical User Interface
************************

TARDIS uses the `PyQt4 framework <http://www.riverbankcomputing.com/software/pyqt/download>`_ for its cross-platform
interface.

The GUI runs through the `IPython Interpreter <http://ipython.org/install.html>`_ which should be started with the
command ``ipython --pylab=qt``, so that it has acess to pylab.

Creating an instance of the :class:`ModelViewer`-class requires that PyQt4/PySide has already been initialized in
IPython. The above command to start IPython accomplishes this.

gui.py contains all the classes used to create the GUI for Tardis.

This module must be imported inside IPython console started above. The console provides the event loop and the place
to create/calculate the tardis model. So the module is basically a tool to visualize results. 

Running instructions
--------------------
    1. Decide which Qt binding you want to use (PySide or PyQt) and 
    accordingly set QT_API in shell
            ```bash
            export QT_API=pyside 
            or
            export QT_API=pyqt
            ``` 
    2. Start the IPython console with eventloop integration 
            ```bash
            ipython --pylab=qt4
            ```
    3. Display your model
            ```python
            from tardis import gui 
            win = gui.Tardis()
            win.show_model(mdl)
            ```
Raises
------
    TemporarilyUnavaliable
        Raised when the currently disabled active mode is requested.

GUI Layout and Features
-----------------------
