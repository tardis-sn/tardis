******
Models
******

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

