*********************
Documentation Preview
*********************

Most of the time it's enough to build documentation locally
and see how things are going. But sometimes it's nice to share
our changes to get some feedback from other collaborators.

Unfortunately, GitHub Pages does not provide a simple way to
preview documentation changes from a pull request as Read The Docs
does, but there's a way to achieve this that suits most of our use
cases.


=========
Procedure
=========

Imagine you are working locally on a branch named `new-feature`.

- Make a new branch from the branch you are working with a suitable
    name, for example: `new-feature-docs`.

- Edit `.github/workflows/documentation-build.yml` and add a new
    trigger:

.. code-block: none
    pull_request:
    branches:
      - master

- Push changes and make a new pull request to **your fork** *master* 
    branch.

- Go to fork *Settings* and make sure you is building GitHub Pages from
    *gh-pages* branch.

- The documentation preview should be available at
    `<your-username>.github.io/tardis`.
