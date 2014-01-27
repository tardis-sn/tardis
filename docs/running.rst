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


**************************
What can I do with TARDIS?
**************************

TARDIS is designed to carry out calculations of synthetic supernova spectra (see Kerzendorf & Sim 2014). A TARDIS calculation requires

1) a model for the supernova eject (density and composition distribution) to be specified;
2) parameters determining the physical and numerical properties of the code are give.

These inputs are specified via the input file (yaml file). As a starting point, examine the input file in the example download (above). The example file sets up a very simple model in which the density distribution is based on the well-known W7 model () and in which uniform abundances are specified. To get you started, here is a short list of settings that can be experimented with to get a feel for using TARDIS. Please contact us if you have questions or woudl like information on more advanced possibilities!

Change the abundances
---------------------

In the example file you can alter the elemental abundances (mass fractions; specified in the "abundances" section). Use this to explore the sensitivity of features to composition.

Change the luminosity
---------------------

You can alter the (emergent) luminosity and see how this affects the synthetic spectrum. Do this by varying the "luminosity requested" in the "supernova" section.

Change the epoch
----------------

In the example file, you can change the epoch for which the synthetic spectrum is calculated (change the value of "time explosion"; specified in the "supernova" section). When doing this you might also change the inner boundary velocity ("start" value in the "velocity" section), and probably the luminosity (see above)!

Experiment with the treatment of line opacity
---------------------------------------------

In the "plasma" section you can change the "line interaction type" between "scatter", "downbranch" and "macroatom" - investigate how important fluorescence is in your synthetic spectrum. (See Kerzendorf & Sim and references therein for the meaning of these settings.)



