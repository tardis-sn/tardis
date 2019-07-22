*****
Model
*****

TARDIS requires a model of the homologously expanding ejecta in order to run a simulation.
A model will include information like the velocity shell structure, abundances, density, etc.
TARDIS offers two ways of specifying the model: either directly in the configuration yaml file
or separately in a model.csvy file. See `here <https://csvy.org/>`_ for an explanation of
the csvy file format and :ref:`here <../running/configuration/config-file>` for a link to the csvy_model schema.

.. note::
    We highly recommend using the cleaner csvy format.

CSVY Model
==========

The TARDIS YAML delimiter for csvy files is ``---``. This means that each csvy model
file has the following structure: The first line of the file is the YAML delimiter,
followed by the YAML portion of the csvy model. A line consisting of only the YAML
delimiter separates the YAML portion of the csvy file from the csv part. The YAML part
of the csvy file is for setting model parameters like **v_inner_boundary** or
**model_density_time_0** while the csv part of the file is for setting profiles of
physical parameters of the model (e.g. abundances, radiative temperature, dilution factor, etc).
If you use the csvy model, then you will need to specify the path to the csvy model file
in the main TARDIS configuration file:

.. code-block::

    csvy_model: /path/to/model.csvy

Example CSVY Model
==================

Below we provide an example model.csvy file.

.. literalinclude:: examples/csvy_full_rad.csvy

Using the Config Model
======================

Although we strongly recommend using the CSVY Model, TARDIS also allows the user
to define custom density and abundance profiles in separate files and reference
these files directly in the main TARDIS configuration file. For further details,
see the following links:

.. toctree::
    :maxdepth: 1

    examples/modelconfig
    examples/densityexp
    examples/densitypow
    examples/densitycust
    examples/abundanceuni
    examples/abundancecust

