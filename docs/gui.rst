.. _gui_explanation:

************************
Graphical User Interface
************************

TARDIS uses the `PyQt4 framework <http://www.riverbankcomputing.com/software/pyqt/download>`_ for its cross-platform
interface.

The GUI runs through the `IPython Interpreter <http://ipython.org/install.html>`_ which should be started with the
command ``ipython --pylab=qt4``, so that it has acess to pylab.

Creating an instance of the :class:`ModelViewer`-class requires that PyQt4/PySide has already been initialized in
IPython. The above command to start IPython accomplishes this.

gui.py contains all the classes used to create the GUI for Tardis.

This module must be imported inside IPython console started above. The console provides the event loop and the place
to create/calculate the tardis model. So the module is basically a tool to visualize results. 

Running Instructions
--------------------
Decide which Qt binding you want to use (PySide or PyQt) and 
accordingly set QT_API in shell::

    export QT_API=pyside
    or:
    export QT_API=pyqt
        
Start the IPython console with eventloop integration::

    ipython --pylab=qt4
        
Display your model::
    
    from tardis import run_tardis
    mdl = run_tardis('yamlconfigfile.yml', 'atomdatafile.h5')
    from tardis.gui import interface  
    interface.show(mdl)

(WIP) You can also run a configuration from the shell. To do this first set 
QT_API in the shell. Then navigate to the tardis installation folder and go to
tardis/tardis/gui. Then use the command::
    
    python interface.py path-to-configuration-file path-to-atomic-data-file 

GUI Layout and Features
-----------------------
When you launch the gui the first screen that you see will be something like the snapshot below. 
You can inspect the plot of shells on this window and use the toggle button to change between the 
plot of dilution factor and the radiation temperature in shells.
    
.. image:: graphics/OpeningScreen.png
    :width: 900

Shell Info
~~~~~~~~~~
To see the abundances in a shell double-click on the shell number in the table. You can follow a 
similar strategy in the tables that appear to bring up the ion and level populations for each shell.
The snapshot below shows all the tables that you can bring up by successively clicking the horizontal
table headers. 

.. image:: graphics/ShellInfo.png
    :width: 900

Line Info
~~~~~~~~~
If you switch to the spectrum tab you can see the spectrum plot and a button to show line info.
Once again you can double-click the table headers to get further information. The snapshot below
shows all the tables that you can bring up in this window.

.. image:: graphics/LineInfo.png
    :width: 900
