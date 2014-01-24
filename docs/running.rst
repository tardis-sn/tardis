.. _running:

**************
Running TARDIS
**************

To run TARDIS requires two files. The atomic database (for more info refer to :ref:`atomic-data`) and a
configuration file (more info at :ref:`config-file`).

Simple Example
==============

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
    tardis_spec = loadtxt('output_spectrum.dat')
    plot(tardis_spec[:,0], tardis_spec[:,1])


Scripting TARDIS
================

.. code-block:: python

    from tardis import model, simulation
    from tardis.io import config_reader

    tardis_config = config_reader.TARDISConfiguration.from_yaml('myconfig.yml')
    radial1d_mdl = model.Radial1DModel(tardis_config)
    simulation.run_radial1d(radial1d_mdl)

Graphical User Interface
========================

A graphical user interface is being developed for this code using QT. We envision to have this available in the next few minor releases.


