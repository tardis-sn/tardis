.. _doc-preview:

*********************
Documentation Preview
*********************

To preview your changes to the documentation please:

#. Enable GitHub Actions to run on your fork (see the *Actions* tab).
#. Make sure your fork is deploying GitHub Pages inside the root of ``gh-pages`` branch (go to *Settings* -> *Pages*).

Then, there are two ways to trigger the build:

#. If the branch you are working on contains the word ``doc`` in it, then every commit pushed to that branch will trigger the build.
#. If your commit message contains the ``[build docs]`` tag, then that commit will trigger the build.

.. note::

    If you forgot to folow any of the two ways described above, you always can trigger push an empty commit
    as follows: ``git commit --allow-empty -m "[build docs]"``


Your preview will be available at ``<username>.github.io/tardis/branch/<branch name>/index.html``.
