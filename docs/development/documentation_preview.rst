.. _doc-preview:

*********************
Documentation Preview
*********************

To preview your changes to the documentation please:

#. Enable GitHub Actions in the *Actions* tab of your fork.
#. Under *Settings -> Pages* in your fork, make sure GitHub Pages is being built from the ``gh-pages`` branch and the ``/ (root)`` folder.

Then, there are two ways to trigger the build:

#. If the branch you are working on contains the word ``doc`` in it, then every commit pushed to that branch will trigger the build.
#. If your commit message contains the ``[build docs]`` tag, then that commit will trigger the build.

.. note::

    You always can trigger a new build by pushing an empty commit: ``git commit --allow-empty -m "[build docs]"``


Your preview will be available at ``<username>.github.io/tardis/branch/<branch name>/index.html``.
