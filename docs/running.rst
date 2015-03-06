.. _running:

**************
Running TARDIS
**************

To run TARDIS requires two files. The atomic database (for more info refer to :ref:`atomic-data`) and a
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

The code brings a GUI that let's users explore the model. Currently, this is
in beta version, but already very usable. The current version needs to be
used with ``PySide``. Please in addition to the normal requirements install
anaconda. Then run the code in the following way

.. code-block:: python

    $ipython --pylab=qt4
    >>> from tardis import run_tardis
    >>> mdl = run_tardis('tardis_example.yml', 'kurucz_cd23_chianti_H_He.h5')
    >>> from tardis import gui
    >>> mdlvwr = gui.ModelViewer()
    >>> mdlvwr.show_model(mdl)

