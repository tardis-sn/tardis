Graphical User Interface
========================

TARDIS uses the `PyQt4 framework <http://www.riverbankcomputing.com/software/pyqt/download>`_ for its cross-platform
interface.

The GUI runs through the `IPython Interpreter <http://ipython.org/install.html>`_ which should be started with the
<<<<<<< HEAD
command ``ipython-2.7 --pylab=qt``, so that it has acess to pylab.

Creating an instance of the :class:`ModelViewer`-class requires that PyQt4 has already been initialized in
IPython. The above command to start IPython accomplishes this.
=======
command ``ipython-2.7 --pylab``, so that it has acess to pylab.

Creating an instance of the :class:`ModelViewer`-class requires that PyQt4 has already been initialized in
IPython with the command ``%gui qt``.
>>>>>>> 7d434b790386d60258c5921c9080ee347d3b3722

GUI Layout and Features
-----------------------


