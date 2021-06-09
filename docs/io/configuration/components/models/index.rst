******
Models
******

CSVY Model
==========

Example CSVY Model
------------------



Model Configuration
===================

Structure
---------

Abundance
---------



Model Converters
================

There are a variety of formats for models from other codes
(both hydro and radiative transfer) that can be converted to TARDIS input files.
Here we aim to provide converters for the most commonly used file formats.

.. toctree::
    :maxdepth: 2

    converters/stella_to_tardis
    converters/cmfgen




    .. jsonschema:: ../schemas/model.yml

Abundances
^^^^^^^^^^
The ``abundance`` section has a possible ``file`` parameter with ``type`` (currently only ``artis`` is allowed)
and a ``name`` parameter giving a path to a file containing the abundance information.

.. warning::
    In contrast to the ``structure`` section, the ``abundance`` section will not ignore abundances set in the rest of
    the section but merely will overwrite the abundances given in the file section.

The rest of the section can be used to configure uniform abundances for all shells, by giving the atom name and a
relative abundance fraction. If it does not add up to 1., TARDIS will warn --- but normalize the numbers.


.. jsonschema:: ../schemas/model_definitions.yml#/definitions/abundances/file

.. jsonschema:: ../schemas/model_definitions.yml#/definitions/abundances/uniform

Structure
^^^^^^^^^
All types of ``structure`` definitions have two keywords in common: ``v_inner_boundary`` and ``v_outer_boundary``.

In the ``structure`` section, one can specify a ``file`` section containing a ``type`` parameter
(currently only ``artis`` is supported``) and a ``name`` parameter giving a path top a file.

If one doesn't specify a ``file`` section, the code requires two sections (``velocities`` and ``densities``) and a
parameter ``no_of_shells``. ``no_of_shells`` is the requested number of shells for a model. The ``velocity`` section
requires a ``type``. Currently, only ``linear`` is supported.

.. jsonschema:: ../schemas/model_definitions.yml#/definitions/structure/file

.. jsonschema:: ../schemas/model_definitions.yml#/definitions/structure/specific

Densities
"""""""""
In the ``densities`` section the ``type`` parameter again decides on the parameters. The type ``uniform`` only needs a
``value`` parameter with a density compatible quantity. The type ``branch85_w7`` uses a seven-order polynomial fit to
the W7 model and is parametrised by time since explosion. The parameters ``time_0`` and ``density_coefficient`` are set
to sensible defaults and should not be changed.

.. jsonschema:: ../schemas/model_definitions.yml#/definitions/density/branch85_w7

.. jsonschema:: ../schemas/model_definitions.yml#/definitions/density/exponential

.. jsonschema:: ../schemas/model_definitions.yml#/definitions/density/power_law

.. jsonschema:: ../schemas/model_definitions.yml#/definitions/density/uniform


TARDIS requires a model of the homologously expanding ejecta in order to run a simulation.
A model will include information like the velocity shell structure, abundances, density, etc.
TARDIS offers two ways of specifying the model: either directly in the configuration YAML file
or separately in a model.csvy file. See `here <https://csvy.org/>`_ for an explanation of
the CSVY file format and :ref:`here <config-file>` for a link to the csvy_model schema.

.. note::
    We highly recommend using the cleaner CSVY format.


The following provides some example setups for a number of typical TARDIS use
cases. You can use these examples as blueprints for your own TARDIS
applications.


Simple Parametrized Models
==========================

These setups specify the ejecta solely via the YAML configuration file:

.. toctree::
    :maxdepth: 2

    examples/tardis_example


API demonstrations
==================

An example on how to use the formal integrator with tardis:

* :ref:`integrator`


Detailed Explosion Models
=========================

Coming soon


CSVY Model
==========

The TARDIS YAML delimiter for CSVY files is ``---``. This means that each CSVY model
file has the following structure: The first line of the file is the YAML delimiter,
followed by the YAML portion of the CSVY model. A line consisting of only the YAML
delimiter separates the YAML portion of the CSVY file from the CSV part. The YAML part
of the CSVY file is for setting model parameters like **v_inner_boundary** or
**model_density_time_0**, while the CSV part of the file is for setting profiles of
physical parameters of the model (e.g. abundances, radiative temperature, dilution factor, etc).
If you use the CSVY model, then you will need to specify the path to the CSVY model file
in the main TARDIS configuration file:

.. jsonschema:: ../schemas/csvy_model.yml

.. code-block::

    csvy_model: /path/to/model.csvy

Example CSVY Model
==================

Below we provide an example model.csvy file.

.. literalinclude:: csvy_full_rad.csvy

Using the Config Model
======================

Although we strongly recommend using the CSVY Model, TARDIS also allows the user
to define custom density and abundance profiles in separate files and reference
these files directly in the main TARDIS configuration file. For further details,
see the following links:

.. toctree::
    :maxdepth: 1

    densityexp/densityexp
    densitypow/densitypow
    densitycust/densitycust
    abundanceuni/abundanceuni
    abundancecust/abundancecust

