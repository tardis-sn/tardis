.. _running:

**************
Running TARDIS
**************

To run TARDIS requires two files. The atomic database (for more info refer to :ref:`atomic-data-download`) and a
configuration file (more info at :ref:`config-file`).

Running TARDIS in the commandline
=================================

After installing TARDIS just download the example directory `<https://www.dropbox.com/s/svvyr5i7m8ouzdt/tardis_example.tar.gz>`_
and run TARDIS with:


.. code-block:: none

    tar zxvf tardis_example.tar.gz
    cd tardis_example
    tardis tardis_example.yml output_spectrum.dat



Then plot the output_spectrum.dat with your favourite plotting program. Here's an example how to do this with python.
(The only thing you need to install is ipython and matplotlib - in addition to TARDIS's requirements)

.. code-block:: python

    ipython --pylab
    >>> tardis_spec = loadtxt('output_spectrum.dat')
    >>> plot(tardis_spec[:,0], tardis_spec[:,1])


More atomic datasets can be downloaded from :ref:`atomic-data-download`.




Running TARDIS interactively
============================

To get more information from each run of TARDIS one can run it interactively and
have full access to the model properties

.. code-block:: python

    >>> from tardis import run_tardis
    >>> model = run_tardis('myconfig.yml')

This model can then be used for inspecting the run as described
:ref:`physical_quantities`


Graphical User Interface
========================
To get a detailed explanation on gui layout go to :ref:`gui_explanation` .

**To setup and run the GUI(under development) follow these steps:**

The gui can use one of two python bindings for qt, namely PyQt4
and PySide. You can choose which binding is used by setting the
environment variable QT_API in your bash.

**1**. Installing required packages

.. code-block:: none
	
	source activate tardis
	conda install ipython=3.0.0 pyside=1.2.1 shiboken=1.2.1


**2**. Choosing between PySide and PyQt4

.. code-block:: none

	#To choose PySide
	export QT_API=pyside
	
	#To choose PyQt
	export QT_API=pyqt

**3**. An example of creating a model and GUI

To show the gui from the ipython shell use the following commands.

.. code-block:: none

	ipython --pylab=qt4

.. code-block:: python

	>>> from tardis import run_tardis
	>>> mdl = run_tardis('yamlconfigfile.yml', 'atomdatafile.h5')
	>>> from tardis.gui import interface 
	>>> interface.show(mdl)

If you just want to run from a configuration file and show the results, you can 
do that outside the ipython shell.To do this navigate to the folder where you 
installed tardis and go to tardis/tardis/gui, and use the following command.

.. code-block:: none

    python interface.py path-to-yaml-configuration-file path-to-atomic-data-file 
