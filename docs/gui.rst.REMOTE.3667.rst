Graphical User Interface
========================

TARDIS uses the `PyQt4 framework <http://www.riverbankcomputing.com/software/pyqt/download>`_ for its cross-platform
interface.

The GUI runs through the `IPython Interpreter <http://ipython.org/install.html>`_ which should be started with the
command ``ipython-2.7 --pylab``, so that it has acess to pylab.

Creating an instance of the :class:`ModelViewer`-class requires that PyQt4 has already been initialized in
IPython with the command ``%gui qt``.

GUI Layout and Features
-----------------------


