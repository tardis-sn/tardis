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
downloaded (XXX Link XXX).

.. code-block:: shell

    > python setup.py test --args="--atomic-dataset=kurucz_cd23_chianti_H_He.h5"


Running the integration tests
=============================

These tests require reference files against which the results of the various
tardis runs are tested. So you first need to either download the current
reference files (XXX link XXX) or generate new ones.

.. code-block::

    > python setup.py



