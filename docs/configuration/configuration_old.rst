.. _config-file:

******************
Configuration File
******************

.. :currentmodule:: tardis.config_reader

TARDIS uses the `YAML markup language <https://en.wikipedia.org/wiki/YAML>`_ for its configuration files. There are several sections which allow different
settings for the different aspects of the TARDIS calculation. An example configuration file can be downloaded
:download:`here <example_configs/full_configuration.yml>`.

.. warning::
    One should note that currently floats in YAML need to be specified in a special format:
    any pure floats need to have a +/- after the e e.g. 2e+5

Every configuration file begins with the most basic settings for the model:

.. code-block:: yaml

    ---
    #Currently only simple1d is allowed
    config_type: simple1d

    #luminosity any astropy.unit convertible to erg/s
    #special unit log_lsun(log(luminosity) - log(L_sun)
    luminosity: 9.44 log_lsun

    #time since explosion
    time_explosion: 13 day

    atom_data: ../atom_data/kurucz_atom_chianti_many.h5

The ``config_type`` currently only allows ``simple1d``, but might be expanded in future versions of TARDIS. Many
parameters in the TARDIS configuration files make use of units. One can use any unit that is supported by
`astropy units <http://docs.astropy.org/en/stable/units/index.html>`_ as well as the unit ``log_lsun``
which is :math:`\log(L) - \log(L_\odot)`. Time since explosion just takes a normal time quantity. ``atom_data`` requires
the path to the HDF5 file that contains the atomic data (more information about the HDF5 file can be found
here :ref:`atomic-data`).

Plasma
^^^^^^

The next configuration block describes the plasma parameters:

.. code-block:: yaml

    plasma:
        initial_t_inner: 10000 K
        initial_t_rad: 10000 K
        disable_electron_scattering: no
        plasma_type: nebular
        #radiative_rates_type - currently supported are lte, nebular and detailed
        radiative_rates_type: detailed
        #line interaction type - currently supported are scatter, downbranch and macroatom
        line_interaction_type : macroatom
        w_epsilon : 1.0e-10


``inital_t_inner`` is temperature of the black-body on the inner boundary. ``initial_t_rad`` is the
radiation temperature for all cells. For debugging purposes and to compare to :term:`synapps` calculations one can
disable the electron scattering. TARDIS will issue a warning that this is not physical.
There are currently two ``plasma_type`` options available: ``nebular`` and ``lte`` which tell TARDIS how to run the
ionization equilibrium and level population calculations (see :ref:`plasmas` for more information).
The radiative rates describe how to calculate the :math:`J_\textrm{blue}` needed for the :ref:`nlte` calculations and
:ref:`macroatom` calculations. There are three options for ``radiative_rates_type``: 1) ``lte`` in which
:math:`J_\textrm{blue} = \textrm{Blackbody}(T_\textrm{rad})`, 2) ``nebular`` in which
:math:`J_\textrm{blue} = W \times \textrm{Blackbody}(T_\textrm{rad})`, 3) ``detailed`` in which the :math:`J_\textrm{blue}`
are calculated using an estimator (this is described in ??????).

TARDIS currently supports three different kinds of line interaction: ``scatter`` - a resonance scattering implementation,
``macroatom`` - the most complex form of line interaction described in :ref:`macroatom` and ``downbranch`` a simplified
version of ``macroatom`` in which only downward transitions are allowed.

Finally, ``w_epsilon`` describes the dilution factor to use to calculate :math:`J_\textrm{blue}` that are 0, which
causes problem with the code (so :math:`J_\textrm{blue}` are set to a very small number).


**NLTE**:

.. code-block:: yaml

    nlte:
        coronal_approximation: True
        classical_nebular: False

The NLTE configuration currently allows setting ``coronal_approximation`` which sets all :math:`J_\textrm{blue}` to 0.
This is useful for debugging with :term:`chianti` for example. Furthermore one can enable 'classical_nebular' to set all
:math:`\beta_\textrm{Sobolev}` to 1. Both options are used for checking with other codes and should not be enabled in
normal operations.

Model
^^^^^

The next sections, describing the model, are very hierarchical. The base level is ``model`` and contains two subsections:
``structure`` and ``abundances``. Both sections can either contain a ``file`` subsection which specifies a file and
file type where the information is stored or a number of other sections.


.. code-block:: yaml

    model:
        structure:
            no_of_shells : 20

            velocity:
                type : linear
                v_inner : 1.1e4 km/s
                v_outer : 2e4 km/s


            density:
                #showing different configuration options separated by comments
                #simple uniform:
                #---------------
    #            type : uniform
    #            value : 1e-12 g/cm^3
                #---------------

                #branch85_w7 - fit of seven order polynomial to W7 (like Branch 85):
                #---------------
                type : branch85_w7
                #value : 1e-12
                # default, no need to change!
                #time_0 : 19.9999584 s
                # default, no need to change!
                #density_coefficient : 3e29
                #---------------


    #        file:
    #            type : artis
    #            name : artis_model.dat
    #            v_lowest: 10000.0 km/s
    #            v_highest: 20000.0 km/s


In the ``structure`` section, one can specify a ``file`` section containing a ``type`` parameter
(currently only ``artis`` is supported``) and a ``name`` parameter giving a path top a file. For the ``artis`` type,
one can specify the inner and outermost shell by giving a ``v_lowest`` and ``v_highest`` parameter. This will result in
the selection of certain shells which will be obeyed in the abundance section as well if ``artis`` is selected there as
well.

.. warning::
    If a ``file`` section is given, all other parameters and sections in the ``structure`` section are ignored!

If one doesn't specify a ``file`` section, the code requires two sections (``velocities`` and ``densities``) and a
parameter ``no_of_shells``. ``no_of_shells`` is the requested number of shells for a model. The ``velocity`` section
requires a ``type``. Currently, only ``linear`` is supported and needs two parameters ``v_inner`` and ``v_outer`` with
velocity values for the inner most and outer most shell.

In the ``densities`` section the ``type`` parameter again decides on the parameters. The type ``uniform`` only needs a
 ``value`` parameter with a density compatible quantity. The type ``branch85_w7`` uses a seven order polynomial fit to
 the W7 model and is parametrised by time since explosion. The parameters ``time_0`` and ``density_coefficient`` are set
 to sensible defaults and should not be changed.


.. code-block:: yaml

    #-- continued from model block before --
        abundances:
            #file:
            #    type : artis
            #    name : artis_abundances.dat

            nlte_species : [Si2]
            C: 0.01
            O: 0.01
            Ne: 0.01
            Mg: 0.01
            Si: 0.45
            S: 0.35
            Ar: 0.04
            Ca: 0.03
            Fe: 0.07
            Co: 0.01
            Ni: 0.01


The ``abundance`` section again has a possible ``file`` parameter with ``type`` (currently only ``artis`` is allowed)
and a ``name`` parameter giving a path to a file containing the abundance information.

.. warning::
    In contrast to the ``structure`` section, the ``abundance`` section will not ignore abundances set in the rest of
    the section, but merely will overwrite the abundances given in the file section.

In this section we also specify the species that will be calculated with our :ref:`nlte` formalism using the
``nlte_species`` parameter (they are specified in a list using astrophysical notation, e.g. [Si2, Ca2, Mg2, H1]).
The rest of the section can be used to configure uniform abundances for all shells, by giving the atom name and a
relative abundance fraction. If it does not add up to 1., TARDIS will warn - but normalize the numbers.

MonteCarlo
^^^^^^^^^^

The ``montecarlo`` section describes the parameters for the MonteCarlo radiation transport and convergence criteria:

.. code-block:: yaml

    montecarlo:
        seed: 23111963171620
        no_of_packets : 2.e+4
        iterations: 100

        convergence_criteria:
            type: specific
            damping_constant: 0.5
            threshold: 0.05
            fraction: 0.8
            hold: 3

    #    convergence_criteria:
    #        type: damped
    #        damping_constant: 0.5
    #        t_inner:
    #            damping_constant: 0.7



The ``seed`` parameter seeds the random number generator first for the creation of the packets
(:math:`\nu` and :math:`\mu`) and then the interactions in the actual MonteCarlo process.
The ``no_of_packets`` parameter can take a float number for input convenience and gives the number of packets normally
used in each MonteCarlo loop. The parameters ``last_no_of_packets`` and ``no_of_virtual_packets`` influence the last run
of the MonteCarlo loop when the radiation field should have converged. ``last_no_of_packets`` is normally higher than
``no_of_packets`` to create a less noisy output spectrum. ``no_of_virtual_packets`` can also be set to greater than 0 to
use the Virtual Packet formalism (reference missing ?????). The ``iterations`` parameter describes the maximum number of
MonteCarlo loops executed in a simulation before it ends. Convergence criteria can be used to make the simulation stop
sooner when the convergence threshold has been reached.

The ``convergence_criteria`` section again has a ``type`` keyword. Two types are allowed: ``damped`` and ``specific``.
All convergence criteria can be specified separately for the three variables for which convergence can be checked
(``t_inner``, ``t_rad``, ``ws``) by specifying subsections in the ``convergence_criteria`` of the same name. These
override then the defaults.


#. ``damped`` only has one parameter ``damping-constant`` and does not check for convergence.

#. ``specific`` checks for the convergence threshold specified in ``threshold``. For ``t_rad`` and ``w`` only a given
    fraction (specified in ``fraction``) has to cross the ``threshold``. Once a convergence  threshold is read, the simulation
    needs to hold this state for ``hold`` number of iterations.

Spectrum
^^^^^^^^

The spectrum section defines the

.. code-block:: yaml

    spectrum:
        start : 500 angstrom
        end : 20000 angstrom
        bins : 1000
        sn_distance : lum_density
        #sn_distance : 10 Mpc

Start and end are given as Quantities with units. If they are given in frequency space they are switched around if
necessary. The number of bins is just an integer. Finally the ``sn_distance`` can either be a distance or the special
parameter ``lum_density`` which sets the distance to :math:`\sqrt{\frac{1}{4 \pi}}` to calculate the luminosity density.


Config Reader
^^^^^^^^^^^^^

The YAML file is read by using a classmethod of the :meth:`~tardis.config_reader.Configuration.from_yaml`.



.. automodapi:: tardis.config_reader


