.. _config-components:

************************
Configuration Components
************************

TARDIS uses the `YAML markup language <https://en.wikipedia.org/wiki/YAML>`_
for its configuration files. There are several sections which allow different
settings for the different aspects of the TARDIS calculation. An example
configuration file (with a small subset of the options that can be specified) can be downloaded at `example_yml`_.

.. warning::
    One should note that currently floats in YAML need to be specified in a special format:
    any pure floats need to have a +/- after the "e", e.g. 2e+5


The TARDIS configuration consists of multiple sections that pertain to certain parts of the code. We will use the
schemas to show what options are available. Our schema mark-up defines names in bold-fat as required and can be seen here:

.. note::

    The following shows all the options (and their default settings) that are available for TARDIS. No other options
    are allowed or available.

The base schema outlines the main components of the TARDIS configuration. Schemas for these components are linked
in the base schema, and important information about the configuration is linked below.

.. jsonschema:: schemas/base.yml

.. toctree::
    :maxdepth: 3

    supernova/index
    atomic/atomic_data
    plasma/index
    models/index
    montecarlo/index
    spectrum/index