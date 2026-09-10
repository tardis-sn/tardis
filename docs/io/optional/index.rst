.. _optional-input:

***************
Optional Inputs
***************

TARDIS also allows other inputs that are passed as keyword arguments into the ``run_tardis`` function.

.. toctree::
    :maxdepth: 1

    tutorial_callback_example
    tutorial_logging_configuration


Additionally, ``run_tardis`` can take in a filepath for the atomic data and a boolean for virtual packet logging. For example:

.. code-block:: python
    
    from tardis.base import run_tardis

    sim = run_tardis('configuration_file', atom_data='filepath_to_atomic_data', virtual_packet_logging=True)

Both of these are also options in the :ref:`configuration file <config-components>`. The option to pass them inside ``run_tardis``
may be removed in a future release of the code.


File logging
============

Set ``debug.log_file`` in a configuration file to write TARDIS log records to a
plain UTF-8 text file. The file output does not contain terminal color escape
sequences and is produced in addition to the usual stream or widget logging.

.. code-block:: yaml

    debug:
      log_file: tardis.log
