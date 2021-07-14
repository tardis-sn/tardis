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

.. _model-config:

Model Configuration
===================

The following schema

Although we strongly recommend using the CSVY Model, TARDIS also allows the user
to define structure and abundance profiles in separate files and reference
these files directly in the main TARDIS configuration file:

.. jsonschema:: ../schemas/model.yml

For an example of this in use, see :ref:`tardis-example`.































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
and abundances in each cell. The column headers (the first row of the CSV part) may contain ??????????????????????????????????????????????????????????????????????????






























.. contents::
    :local:


.. _csvy-model:

CSVY Model
==========

TARDIS allows users to use a CSVY file to input information about the model:

.. jsonschema:: ../schemas/csvy_model.yml

The TARDIS YAML delimiter for CSVY files is ``---``. This means that each CSVY model
file has the following structure: The first line of the file is the YAML delimiter,
followed by the YAML portion of the CSVY model. A line consisting of only the YAML
delimiter separates the YAML portion of the CSVY file from the CSV part. The YAML part
of the CSVY file is for setting model parameters like **v_inner_boundary** or
**model_density_time_0**, while the CSV part of the file is for setting profiles of
physical parameters of the model (e.g. abundances, radiative temperature, dilution factor, etc).
If you use the CSVY model, then you will need to specify the path to the CSVY model file
in the main TARDIS configuration file.

Example CSVY Model
------------------

We provide an example CSVY model file:

.. literalinclude:: csvy_full_rad.csvy


Structure, Density, and Abundance
=================================

Structure
---------

All types of ``structure`` definitions have two keywords in common: ``v_inner_boundary`` and ``v_outer_boundary``.

In the ``structure`` section, one can specify a ``file`` section containing a ``type`` parameter
(currently only ``artis`` is supported``) and a ``name`` parameter giving a path top a file.

If one doesn't specify a ``file`` section, the code requires two sections (``velocities`` and ``densities``) and a
parameter ``no_of_shells``. ``no_of_shells`` is the requested number of shells for a model. The ``velocity`` section
requires a ``type``. Currently, only ``linear`` is supported.

.. jsonschema:: ../schemas/model_definitions.yml#/definitions/structure/file

.. jsonschema:: ../schemas/model_definitions.yml#/definitions/structure/specific


Density
-------

In the ``densities`` section the ``type`` parameter again decides on the parameters. The type ``uniform`` only needs a
``value`` parameter with a density compatible quantity. The type ``branch85_w7`` uses a seven-order polynomial fit to
the W7 model and is parametrised by time since explosion. The parameters ``time_0`` and ``density_coefficient`` are set
to sensible defaults and should not be changed.

.. jsonschema:: ../schemas/model_definitions.yml#/definitions/density/branch85_w7

.. jsonschema:: ../schemas/model_definitions.yml#/definitions/density/exponential

.. toctree::
    :maxdepth: 1
    
    densityexp/densityexp

    
.. jsonschema:: ../schemas/model_definitions.yml#/definitions/density/power_law

.. toctree::
    :maxdepth: 1
    
    densitypow/densitypow

    
.. jsonschema:: ../schemas/model_definitions.yml#/definitions/density/uniform

Custom Density
^^^^^^^^^^^^^^

TARDIS also allows for a custom density profile:

.. toctree::
    :maxdepth: 1
    
    densitycust/densitycust

    
Abundance
---------

The ``abundance`` section has a possible ``file`` parameter with ``type`` (currently only ``artis`` is allowed)
and a ``name`` parameter giving a path to a file containing the abundance information.

.. warning::
    In contrast to the ``structure`` section, the ``abundance`` section will not ignore abundances set in the rest of
    the section but merely will overwrite the abundances given in the file section.

The rest of the section can be used to configure uniform abundances for all shells, by giving the atom name and a
relative abundance fraction. If it does not add up to 1., TARDIS will warn --- but normalize the numbers.


.. jsonschema:: ../schemas/model_definitions.yml#/definitions/abundances/file

.. toctree::
    :maxdepth: 1
    
    abundancecust/abundancecust


.. jsonschema:: ../schemas/model_definitions.yml#/definitions/abundances/uniform

.. toctree::
    :maxdepth: 1
    
    abundanceuni/abundanceuni


Tutorials
=========

.. toctree::
    :maxdepth: 2

    Custom_TARDIS_Model_Tutorial


.. _converters:

Model Converters
================

There are a variety of formats for models from other codes
(both hydro and radiative transfer) that can be converted to TARDIS input files.
Here we aim to provide converters for the most commonly used file formats.

.. toctree::
    :maxdepth: 2

    converters/stella_to_tardis
    converters/cmfgen