*************
Running tests
*************

There are two basic categories of tests in TARDIS: 1) the unit
tests, and 2) the integration tests. Unit tests check the outputs of individual functions,
while the integration tests check entire runs for different setups of TARDIS.

The unit tests run very quickly and thus are executed after every suggested change
to TARDIS. The integration tests are much more costly and thus are only executed
every few days on a dedicated server.

All of them are based on the excellent ``astropy-setup-helpers`` package and
`pytest <https://docs.pytest.org/en/latest/>`_.

Running the Unit Tests
======================

This is very straightforward to run on your own machine. For very simple unit
tests, you can run this with:

.. code-block:: shell

    > python setup.py test


Running the more advanced unit tests requires TARDIS Reference data that can be
downloaded
(`tardis_refdata <https://github.com/tardis-sn/tardis-refdata>`_).

.. code-block:: shell

    > python setup.py test --args="--tardis-refdata=/path/to/tardis-refdata/"

Generating Plasma Reference
===========================

You can generate Plasma Reference by the following command:

.. code-block:: shell

    > pytest -rs tardis/plasma/tests/test_complete_plasmas.py 
    --tardis-refdata="/path/to/tardis-refdata/" --generate-reference

Running the Integration Tests
=============================

These tests require reference files against which the results of the various
tardis runs are tested. So you first need to either download the current
reference files (`here <https://github.com/tardis-sn/tardis-refdata>`_)
or generate new ones.

Both of these require a configuration file for the integration tests:

.. literalinclude:: integration.yml
    :language: yaml

Inside the atomic data directory there needs to be atomic data for each of
the setups that are provided in the ``test_integration`` folder.
If no references are given, the first step is to generate them.
The ``--less-packets`` option is useful for debugging purposes and will just
use very few packets to generate the references and thus make the process much
faster --- THIS IS ONLY FOR DEBUGGING PURPOSES. The ``-s`` option ensures that
TARDIS prints out the progress:

.. code-block:: shell

    > python setup.py test --args="--integration=integration.yml -m integration
    --generate-reference --less-packets"

To run the test after having run the ``--generate-references``, all that is
needed is:

.. code-block:: shell

    > python setup.py test --args="--integration=integration.yml -m integration
    --less-packets" --remote-data


Setting up the DokuWiki report
==============================

A normal `DokuWiki <https://www.dokuwiki.org/dokuwiki>`_ installation is performed on the required server. Before the
connection works one is requires to set the option remote access in the
settings. If this is not done the ``dokuwiki`` python plugin will not connect
with the warning ``DokuWikiError: syntax error: line 1, column 0``. One also
has to enable this for users (``remoteuser`` option) otherwise the error:
``ProtocolError for xmlrpc.php?p=xxxxxx&u=tardistester: 403 Forbidden``
will appear.

Another important configuration option is to enable embedded html ``htmlok``;
otherwise, it won't show nice html page reports.

Finally, one has to call the `python setup.py test` with the ``--remote-data``
option to allow posting to an external DokuWiki server.



