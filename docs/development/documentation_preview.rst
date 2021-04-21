.. _doc-preview:

*********************
Documentation Preview
*********************

Most of the time it's enough to build the documentation locally and see how things are going. But sometimes 
it's nice to share our changes and get some feedback from other collaborators. 

To preview your changes to the documentation first make sure you have GitHub Actions activated on your fork
(see the Actions tab). Then, there are two ways to trigger the build:

#. If the branch you are working on contains the word ``doc`` in it, then every commit will trigger the build.
#. Add the tag ``[build docs]`` on the commit message before pushing.

Your preview will be available at ``<username>.github.io/tardis/branch/<branch name>/index.html``.
