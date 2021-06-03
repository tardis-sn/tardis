Writing tests
=============

pytest has the following test discovery rules:

 * test_*.py or \*_test.py files

 * Test prefixed classes (without an __init__ method)

 * test\_ prefixed functions and methods

It is recommended to use the `test_*.py` pattern for consistency.

Tests Location
==============

######################
Package specific tests
######################

Each package should include a suite of unit tests, covering as many of the public methods/functions as possible. These tests should be included inside each sub-package, e.g:

.. code-block:: shell

    tardis/montecarlo/tests/

##########
Full tests
##########

Tests that check the full functionality of TARDIS are kept in the `tests` directory

.. code-block:: shell

    tardis/tests/