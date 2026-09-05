.. _atomic-data-download:

***********
Atomic Data
***********

Atomic data is provided as a filepath in the configuration file or as a keyword argument in the ``run_tardis``
function. The `Carsus package <https://tardis-sn.github.io/carsus/>`_ can generate atomic data in a format readable
by TARDIS.

The default public atomic dataset can be downloaded with
``tardis.io.atom_data.atom_web_download.download_atom_data``. This stores the
dataset in TARDIS's local data directory for use in a configuration.

.. toctree::

    atomic_data_description
    current_public_table
