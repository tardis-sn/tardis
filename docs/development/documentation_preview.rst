.. _doc-preview:

*********************
Documentation Preview
*********************

To preview your changes to the documentation please:

#. Enable GitHub Actions to run on your fork (see the *Actions* tab).
#. Make sure your fork is deploying GitHub Pages inside the root of ``gh-pages`` branch (see the *Settings* tab).

Then, there are two ways to trigger the build:

#. If the branch you are working on contains the word ``doc`` in it, then every commit will trigger the build.
#. Add the tag ``[build docs]`` on the commit message before pushing.

Your preview will be available at ``<username>.github.io/tardis/branch/<branch name>/index.html``.
