.. _running:

****************
Quickstart Guide
****************

To run TARDIS requires two files. The atomic database (for more info refer to
:ref:`atomic-data-download`) and a configuration file (more info at
:ref:`config-file`). You can obtain a copy of the atomic database from the
`tardis-refdata <https://github.com/tardis-sn/tardis-refdata>`_ repository
(``atom_data`` subfolder). We recommended to use the
``kurucz_cd23_chianti_H_He.h5`` dataset.


Running TARDIS interactively in a Jupyter notebook (recommended)
================================================================

To get more information from each run of TARDIS one can run it interactively and
have full access to the model properties (as described in :ref:`physical_quantities`)

.. toctree::
    :maxdepth: 1

    examples/run_simple_example.ipynb


Running TARDIS in the commandline
=================================

After installing TARDIS just download the configuration file from the
`tardis-setups <https://github.com/tardis-sn/tardis-setups>`_ and the standard
atomic data set from the `tardis-refdata
<https://github.com/tardis-sn/tardis-refdata>`_ repository and run TARDIS.
Assuming you have ``wget``, you could follow the procedure:


.. code-block:: none

    mkdir tardis_example
    cd tardis_example
    wget https://raw.githubusercontent.com/tardis-sn/tardis-setups/master/2014/2014_kerzendorf_sim/appendix_A1/tardis_example.yml
    wget https://github.com/tardis-sn/tardis-refdata/raw/master/atom_data/kurucz_cd23_chianti_H_He.h5
    tardis tardis_example.yml output_spectrum.dat


Then plot the output_spectrum.dat with your favourite plotting program. Here's an example how to do this with python.
(The only thing you need to install is ipython and matplotlib - in addition to TARDIS's requirements)

.. code-block:: python

    ipython --pylab
    >>> tardis_spec = loadtxt('output_spectrum.dat')
    >>> plot(tardis_spec[:,0], tardis_spec[:,1])

More atomic datasets can be downloaded from :ref:`atomic-data-download`.


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
