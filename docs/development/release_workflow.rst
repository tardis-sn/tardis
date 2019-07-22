*******
Release
*******

This is documentation is mainly intended for core developers of TARDIS. There
are several steps to do a release.


Update the Changelog
====================

The changelog is located in ``CHANGELOG.rst`` in the root directory of TARDIS
and contains a description of changes between versions. One useful idea is to
mention all the pull requests, which can be accomplished with the tool
`github-changes <https://www.npmjs.com/package/github-changes>`_ and invoke it
with

.. code:: shell

  github-changes -o tardis-sn -r tardis --only-pulls --use-commit-body -f tmp_changelog.md

which will write the pull requests into ``tmp_changelog.md``. This can be used
to put into the actual ``CHANGELOG.rst``.

Doing the actual release
========================

This is best done in a clean environment (so a fresh checkout from the master).
We again use the astropy instructions:
`<http://astropy.readthedocs.org/en/latest/development/releasing.html>`_
