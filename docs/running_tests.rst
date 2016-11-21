*************
Running tests
*************

There are two basic categories of tests unit tests in TARDIS 1) the unit
tests 2) integration tests. Unit tests check the outputs of individual functions
while the integration tests check entire runs for different setups of TARDIS.

The Unit tests run very quickly and thus are executed after every suggested change
to TARDIS. The Integration tests are much more costly and thus are only executed
every few days on a dedicated server.

All of them are based on the excellent ``astropy-setup-helpers`` package and
``pytest``.

Running the unit tests
======================

This is very straight forward to run on your own machine. For very simple unit
tests you can run this with:

.. code-block:: shell

    > python setup.py test


Running the more advanced unit tests it requires atomic data that can be
downloaded
(`atom_data <http://opensupernova.org/files/atomic/kurucz_cd23_chianti_H_He.h5.zip>`_).

.. code-block:: shell

    > python setup.py test --args="--atomic-dataset=kurucz_cd23_chianti_H_He.h5"


Running the integration tests
=============================

These tests require reference files against which the results of the various
tardis runs are tested. So you first need to either download the current
reference files (`here <https://github.com/tardis-sn/tardis-refdata>`_)
or generate new ones.

Both of of these require a configuration file for the integration tests:

.. literalinclude:: integration.yml
    :language: yaml

Inside the atomic data directory there needs to be atomic data for each of
the setups that are provided in the ``test_integration`` folder.
If no references are given the first step is to generate them.
The ``--less-packets`` option is useful for debugging purposes and will just
use very few packets to generate the references and thus make the process much
faster - THIS IS ONLY FOR DEBUGGING PURPOSES. The ``-s`` option ensures that
TARDIS prints out the progress:

.. code-block::

    > python setup.py test --args="--integration=integration.yml -m integration
    --generate-reference --less-packets -s"

