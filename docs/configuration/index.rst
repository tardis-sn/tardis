********************
TARDIS Configuration
********************

Most of the configuration for TARDIS uses the YAML format. If a YAML file is
handed to TARDIS it is first checked using a configuration validator to
see if the required keywords of the required type exist and produces a dictionary
with that has the validated types as well as sensible defaults for missing values
built-in. This validator uses a default validation file that already explains
a lot about the required configuration items.

.. toctree::
    :maxdepth: 1

    configuration
    config_validator