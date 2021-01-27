.. _doc-preview:

*********************
Documentation Preview
*********************

Most of the time it's enough to build the documentation locally and see how things are going. But sometimes 
it's nice to share our changes and get some feedback from other collaborators. 

Unfortunately, GitHub Pages does not provide a simple way to preview documentation changes from pull requests
(as ReadTheDocs does) but there's a way to achieve this that suits for most of our use cases.


=========
Procedure
=========

Imagine you are developing a new feature for TARDIS on your local
branch named ``new-feature`` and follow these steps:

1. Go to your fork's *Settings* tab make sure GitHub Pages are building from the *gh-pages* branch.

2. Checkout to a new branch with a suitable name, like ``new-feature-docs``.

3. Edit ``.github/workflows/documentation-build.yml`` and add a new trigger below the *push* trigger::

    pull_request:
      branches:
        - master

4. Push changes and make a new pull request to **your fork's** *master* branch.

5. If everything is ok, the documentation preview should be available at ``<your-username>.github.io/tardis``.

.. note :: Remember you will need to rebase ``new-feature-docs`` to ``new-feature`` every time you push changes to ``new-feature``.


===========
Limitations
===========

This method has a major drawback: you can build just a single preview for your entire fork. This means if
your are working on multiple pull request, your site ``<your-username.github.io/tardis>`` will display just
the latest successful build, overwriting the previous one.
