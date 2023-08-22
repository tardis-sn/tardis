.. _config-file:

******************************
Configuration (Required Input)
******************************


The necessary parameters for a TARDIS simulation are provided via a `YAML <https://en.wikipedia.org/wiki/YAML>`_
configuration file. The following sections respectively describes what can or must be included in the
configuration file, shows an example configuration file, describes how TARDIS checks to make sure a configuration
file is valid, and demonstrates how a YAML configuration file is read in.

.. caution:: 
    TARDIS parallelization is not working correctly at the moment and might produce incorrect results.
    Please avoid using it.
    For more information, see issue [#2021](https://github.com/tardis-sn/tardis/issues/2021).
  
.. toctree::
  :maxdepth: 1

  components/index
  example
  config_validator
  read_configuration
