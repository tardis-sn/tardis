******************
Configuration File
******************

TARDIS uses the `YAML markup language <https://en.wikipedia.org/wiki/YAML>`_
for its configuration files. There are several sections which allow different
settings for the different aspects of the TARDIS calculation. An example
configuration file can be downloaded :download:`here
<../examples/tardis_example.yml>`.

.. warning::
    One should note that currently floats in YAML need to be specified in a special format:
    any pure floats need to have a +/- after the e e.g. 2e+5


The TARDIS configuration consists of multiple sections that pertain to certain parts of the code. We will use the
schemas to show what options are available. Our schema mark-up defines names in bold-fat as required.
can be seen here:

.. note::

    The following shows all the options (and their default settings) that are available for TARDIS. No other options
    are allowed or available

Configuration Schema
====================

.. jsonschema:: schemas/base.yml


The base schema links to several other configuration parts


Supernova
---------

.. jsonschema:: schemas/supernova.yml


Model
-----

The next sections, describing the model, are very hierarchical. The base level is ``model`` and contains two subsections:
``structure`` and ``abundances``. Both sections can either contain a ``file`` subsection which specifies a file and
file type where the information is stored or a number of other sections.

.. jsonschema:: schemas/model.yml

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

The ``abundance`` section again has a possible ``file`` parameter with ``type`` (currently only ``artis`` is allowed)
and a ``name`` parameter giving a path to a file containing the abundance information.

.. warning::
    In contrast to the ``structure`` section, the ``abundance`` section will not ignore abundances set in the rest of
    the section, but merely will overwrite the abundances given in the file section.

In this section we also specify the species that will be calculated with our :ref:`nlte` formalism using the
``nlte_species`` parameter (they are specified in a list using astrophysical notation, e.g. [Si2, Ca2, Mg2, H1]).
The rest of the section can be used to configure uniform abundances for all shells, by giving the atom name and a
relative abundance fraction. If it does not add up to 1., TARDIS will warn - but normalize the numbers.



Plasma
------

.. jsonschema:: schemas/plasma.yml

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

NLTE
^^^^

.. code-block:: yaml

    nlte:
        coronal_approximation: True
        classical_nebular: False

The NLTE configuration currently allows setting ``coronal_approximation`` which sets all :math:`J_\textrm{blue}` to 0.
This is useful for debugging with :term:`chianti` for example. Furthermore one can enable 'classical_nebular' to set all
:math:`\beta_\textrm{Sobolev}` to 1. Both options are used for checking with other codes and should not be enabled in
normal operations.



Monte Carlo
-----------

The ``montecarlo`` section describes the parameters for the MonteCarlo radiation transport and convergence criteria:

.. jsonschema:: schemas/montecarlo.yml

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
--------

Start and end are given as Quantities with units. If they are given in
frequency space they are switched around if necessary. The number of bins is
just an integer. Finally, the method option selects the final spectral synthesis mode. Currently, there are three options:

* real: construct spectrum from the real packet population alone
* virtual: use the :doc:`virtual packet scheme <../montecarlo/virtualpackets>` for spectral synthesis
* integrated: use the :doc:`formal integral method <../montecarlo/sourceintegration>` of Lucy 1999

.. warning::

  Currently, the "integrated" mode only works with the downbranching line
  interaction mode. Note also the limitations listed at the bottom of the
  dedicated page.

.. jsonschema:: schemas/spectrum.yml