*********************
Documentation Preview
*********************

Most of the time it's enough to build the documentation locally
and see how things are going. But sometimes it's nice to share
our changes and get some feedback from other collaborators.

Unfortunately, GitHub Pages does not provide a simple way to
preview documentation changes from a pull request (as ReadTheDocs
does) but there's a way to achieve this that suits for most of our 
use cases.


=========
Procedure
=========

Imagine you are developing a new feature for TARDIS on your local
branch named `new-feature` and follow these steps:

- Checkout to a new branch with a suitable name, 
    like `new-feature-docs`.

- Edit `.github/workflows/documentation-build.yml` and add a new
    trigger below the *push* trigger:

.. code-block: none
    pull_request:
      branches:
        - master

- Push changes and make a new pull request to **your fork's** *master* 
    branch.

- Go to your fork's *Settings* tab make sure GitHub Pages are building
    from the *gh-pages* branch.

- If everything is okm the documentation preview should be available at
    `<your-username>.github.io/tardis`.

- Remember you will need to rebase `new-feature-docs` to `new-feature`
    every time you make changes to `new-feature`.
