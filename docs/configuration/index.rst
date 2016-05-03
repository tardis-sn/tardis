********************
TARDIS Configuration
********************

Overview
========

The bulk of the parameters for a TARDIS simulation are provided via a YAML
configuration file. In a first step, the integrity of the specified
configuration file is checked by a configuration validator. During this
process, missing entries are replaced by default values. The following sections
and pages provide more details concerning the TARDIS configuration process.

Configuration Details
=====================

* :doc:`In depth description of config options <configuration_old>`
* :doc:`Validation Process <config_validator>`
* :doc:`Example Configuration with Default Values <configuration>`


.. warning::
    We are currently in the process of overhauling the configuration system.
    See `TEP 006 <https://github.com/tardis-sn/tep/blob/master/TEP006_configuration_tags.rst>`_
    and stay tuned!

..
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
