.. _documentation-guidelines:

************************
Documentation Guidelines
************************

High-quality and consistent documentation is very important at TARDIS. It allows new users to find out how to do something specific using TARDIS, as well as helps developers (like you!) to understand the best practices.

TARDIS uses the popular Python documentation generator `Sphinx <https://www.sphinx-doc.org/>`_. Sphinx translates a set of plain text source files (often written in `reStructuredText <https://www.sphinx-doc.org/en/master/usage/restructuredtext/basics.html>`_) to HTML files, automatically producing cross-references, indices, etc. If you haven't worked with Sphinx before, you should first read their `quickstart <https://www.sphinx-doc.org/en/master/usage/quickstart.html>`_ guide.


Documenting the code you write
==============================

When making or adding changes to the functionality of an aspect of TARDIS, an example notebook or .rst file should be created to demonstrate how it works. A good example of this is the `quickstart notebook <https://tardis-sn.github.io/tardis/quickstart/quickstart.html>`_ which gives an overview of how TARDIS works. You will also need to add the .rst (or notebook) file's path within the relevant ``index.rst`` file - it is usually present within the same (or parent) directory where you've created your file. For example, when the documentation was created for `research papers that used TARDIS <https://tardis-sn.github.io/tardis/research/research_done_using_TARDIS/research_papers.html>`_, the toctree in `this index file <https://github.com/tardis-sn/tardis/blob/master/docs/research/index.rst>`_ was edited to add the path to the new .rst file.

Besides this, the functions and classes in your code must always contain **docstrings**. Read :ref:`this section <docstrings>` of our code quality guidelines to understand their importance and how they should be formatted. Sphinx uses these docstrings to auto-generate the `API documentation <https://tardis-sn.github.io/tardis/api/modules.html>`_ for entire TARDIS package. Please make sure that you have correctly formatted the docstrings by checking how the corresponding module's API looks once you build the documentation.


Building documentation locally
==============================

To build TARDIS documentation locally, use the following commands:

.. code ::

    cd docs
    make html

.. note :: 

    - If you're working on a fresh local copy of the TARDIS repository, you might need to do ``python setup.py develop`` before executing these commands.
    - Use ``DISABLE_NBSPHINX=1 make html`` to disable notebook rendering (fast mode).

After running this command, you can find the built docs (i.e. HTML webpages) in ``docs/_build``. Open the ``index.html`` in your browser to see how the documentation looks like with your edits. Navigate to page where you made changes or file that you added to check whether it looks as intended or not.


Sharing the built documentation in your PR
==========================================

When you make edits in TARDIS documentation and submit a PR, we can only see the changes in source files in GitHub files diff, but not the built documentation (webpages). This is usually fine unless you have made changes in the way documentation pages are structured or anything that affects lot of files. In such cases, you should share the preview of documentation with your changes by building it online. The steps to do this are described :ref:`here <doc-preview>`. This will help us (the reviewers of your PR) to check how the documentation will look once your PR is merged.
