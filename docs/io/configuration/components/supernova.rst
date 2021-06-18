.. _supernova-config:

***********************
Supernova Configuration
***********************

The supernova component of the configuration file contains some key information about the supernova being modeled:

.. jsonschema:: schemas/supernova.yml

During a run of TARDIS, we attempt to converge the output spectrum to match the requested luminosity
(see :ref:`est_and_conv`).