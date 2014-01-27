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

1) that a model for the supernova eject (density and composition distribution) is specified;
2) parameters determining the physical and numerical properties of the code are give.

These inputs are specified via the input file (yaml file). As a starting point, examine the input file in the example download (above):

.. code-block:: python

    #New configuration for TARDIS based on YAML
    #IMPORTANT any pure floats need to have a +/- after the e e.g. 2e+5
    #Hopefully pyyaml will fix this soon.
    ---
    #Currently only simple1d is allowed
    tardis_config_version: v1.0
    supernova:
        luminosity_requested: 9.44 log_lsun
        time_explosion: 13 day
    #    distance : 24.2 Mpc

    atom_data: kurucz_atom_pure.h5

    model:
        structure:
            type: specific


            velocity:
                start : 1.1e4 km/s
                stop : 20000 km/s
                num: 20
        
        
            density:
               type : branch85_w7
            
        abundances:
            type: uniform
            O: 0.19
            Mg: 0.03
            Si: 0.52
            S: 0.19
            Ar: 0.04
            Ca: 0.03
                
    plasma:
        disable_electron_scattering: no
        ionization: lte
        excitation: lte
        #radiative_rates_type - currently supported are lte, nebular and detailed
        radiative_rates_type: dilute-blackbody
        #line interaction type - currently supported are scatter, downbranch and macroatom
        line_interaction_type: macroatom

    montecarlo:
        seed: 23111963
        no_of_packets : 1.0e+5
        iterations: 20
    
        black_body_sampling:
            start: 1 angstrom
            stop: 1000000 angstrom
            num: 1.e+6
        last_no_of_packets: 1.e+5
        no_of_virtual_packets: 10

        convergence_criteria:
            type: specific
            damping_constant: 1.0
            threshold: 0.05
            fraction: 0.8
            hold: 3
            t_inner:
                damping_constant: 1.0
    

    spectrum:
        start : 500 angstrom
        stop : 20000 angstrom
        num: 10000
    

This input file sets up a very simple model in which the density distribution is based on the well-known W7 model () and in which uniform abundances are specified. To get you started, here is a short list of settings that can be experimented with to get a feel for using TARDIS. Please contact us if you have questions or woudl like information on more advanced possibilities!

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



