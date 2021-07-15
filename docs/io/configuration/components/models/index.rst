.. _model-csvy-and-config:

******
Models
******

TARDIS requires a model of the ejecta in order to run a simulation. A model typically includes information
for the velocity shell structure, density, and abundances. **See :ref:`model` for more information on the
TARDIS model**. TARDIS offers two ways of specifying the model: either directly in the configuration YAML file
or separately in a CSVY file. See `here <https://csvy.org/>`_ for an explanation of the CSVY file format.

TARDIS has several built-in models for the shell structure, density, and abundance. If only these are being used,
we recommend using the YAML configuration method. For creating a custom model (which will be discussed below), we
reccomend using the CSVY method.

.. contents::
    :local:

.. _model-config:

Model Configuration
===================

The following schema shows how the model section of the YAML configuration is set up:

.. jsonschema:: ../schemas/model.yml

For an example of this in use, see :ref:`tardis-example`. This configuration allows for both built-in
structures and abundances, which will be described below, as well as custom models for the shell structure
and abundances. For the latter, while we recommend using the CSVY model, one can specify custom structures
and abundances using other files:

.. jsonschema:: ../schemas/model_definitions.yml#/definitions/structure/file

For more information on a custom structure/density file, see:

.. toctree::
    :maxdepth: 1
    
    densitycust/densitycust

.. jsonschema:: ../schemas/model_definitions.yml#/definitions/abundances/file

For more information on a custom abundence file, see:

.. toctree::
    :maxdepth: 1
        
    abundancecust/abundancecust


Custom Model Tutorial
---------------------

.. toctree::
    :maxdepth: 2

    Custom_TARDIS_Model_Tutorial


.. _converters:

Model Converters
----------------
    
There are a variety of formats for models from other codes
(both hydro and radiative transfer) that can be converted to TARDIS input files.
Here we aim to provide converters for the most commonly used file formats.
    
.. toctree::
    :maxdepth: 2

    converters/stella_to_tardis
    converters/cmfgen




Structure, Density, and Abundance
=================================

TARDIS's built-in models are described in the following sections:

Structure
---------

When using the built-in structure functionality, the code requires two sections (``velocities`` and ``densities``) and a
parameter ``no_of_shells``. ``no_of_shells`` is the requested number of shells for a model. The ``velocity`` section
requires a ``type``. Currently, only ``linear`` is supported.

.. jsonschema:: ../schemas/model_definitions.yml#/definitions/structure/specific


Density
-------

In the ``densities`` section the ``type`` parameter again decides on the parameters. The type ``uniform`` only needs a
``value`` parameter with a density compatible quantity. The type ``branch85_w7`` uses a seven-order polynomial fit to
the W7 model and is parametrised by time since explosion. The parameters ``time_0`` and ``density_coefficient`` are set
to sensible defaults and should not be changed. The physics of these density models are further discussed in
:ref:`model`.

.. jsonschema:: ../schemas/model_definitions.yml#/definitions/density/branch85_w7

.. jsonschema:: ../schemas/model_definitions.yml#/definitions/density/exponential

For more information, see:

.. toctree::
    :maxdepth: 1
    
    densityexp/densityexp

    
.. jsonschema:: ../schemas/model_definitions.yml#/definitions/density/power_law

For more information, see:

.. toctree::
    :maxdepth: 1
    
    densitypow/densitypow

    
.. jsonschema:: ../schemas/model_definitions.yml#/definitions/density/uniform

    
Abundance
---------

This section can be used to configure uniform abundances for all shells, by giving the atom name and a
relative abundance fraction. If it does not add up to 1., TARDIS will warn --- but normalize the numbers.

.. jsonschema:: ../schemas/model_definitions.yml#/definitions/abundances/uniform

For more information, see:

.. toctree::
    :maxdepth: 1
    
    abundanceuni/abundanceuni


.. _csvy-model:

CSVY Model
==========

TARDIS allows users to use a CSVY file to input information about the model. To do this, instead of a
``model`` section, one includes ``csvy_model: <file path to CSVY file>`` in the main TARDIS configuration
file.

The CSVY model has a YAML part as well as a CSV part, separated by the YAML delimiter ``---``. This means
that each CSVY model file has the following structure: The first line of the file is the YAML delimiter,
followed by the YAML portion of the CSVY model, then another line with just the YAML delimiter,
and finally the CSV portion of the CSVY file. This is shown in the example CSVY file later in this section.

The YAML portion of the CSVY file allows the user to use most of the features of the YAML model configuration,
as shown in the schema below:

.. jsonschema:: ../schemas/csvy_model.yml

The CSV part of the CSVY file creates a table that can include information about shell velocities, densities,
and abundances in each cell. The column headers (the first row of the CSV part) may contain ``velocity``,
``density``, ``t_rad``, ``dilution_factor``, or the name of any element or isotope (e.g. ``H``, ``Mg``,
``Ni56``). These columns are explained in the following example:

.. literalinclude:: csvy_full_rad.csvy

Notice that for each column that is used in the CSV section of the file, there is a corresponding field under
``datatype`` in the YAML section of the file. In our example, each of the fields under ``datatype`` has a brief
description to go along with it. While the description is not necessary for any of the fields, the unit section
is required for ``velocity``, ``density``, and ``t_rad``.

Since the ``velocity`` column contains the outer shell velocity, the first entry in the velocity column is the
velocity of the photosphere -- i.e. the inner boundary of the computational domain (see :ref:`model`).
Consequently, **none of the other information in the first row is used**. In our example, there are only two
shells, and the shell will have an inner boundary with a velocity of :math:`9000 km/s`, an outer boundary with a velocity
of :math:`10500 km/s`, a density of :math:`2.0*10^{-10} g/cm^3`, a dilution factor of .8, etc.

.. note::

    None of the CSV columns are required. However, if ``velocity``, ``density``, or the abundances are missing,
    they must be specified in the YAML portion of the file. If ``t_rad`` or ``dilution_factor`` are missing,
    they will be automatically calculated (see :ref:`plasma`).

.. note::

    ``t_rad`` and ``dilution_factor`` are the values of the temperature and dilution factor for the first
    iteration, and will be updated in subsequent iterations (see :ref:`est_and_conv`). To prevent these
    quantities from being changed, you must set the damping constant to zero in the :ref:`Convergence
    Configuration <conv-config>` in the Monte Carlo section of the configuration file.

CSVY Model Tutorial
-------------------

Coming soon!
