.. _documentation-guidelines:

************************
Documentation Guidelines
************************

High-quality and consistent documentation is very important at TARDIS. It allows new users to find out how to do something specific using TARDIS, as well as helps developers (like you!) to understand the best practices.

TARDIS uses the popular Python documentation generator `Sphinx <https://www.sphinx-doc.org/>`_. Sphinx translates a set of source files (often written in reStructuredText or Jupyter notebooks, see below) to HTML files, automatically producing cross-references, indices, etc. If you haven't worked with Sphinx before, you should first read their `quickstart <https://www.sphinx-doc.org/en/master/usage/quickstart.html>`_ guide.


Documenting the code you write
==============================

When making or adding changes to the functionality of an aspect of TARDIS, an ``.rst`` file or Jupyter notebook (``.ipynb`` file) should be created to demonstrate how it works, and that page must then be included in the documentation. This is described in detail in the following sections.


RST Documentation
-----------------

Documentation not featuring interactive code examples is written in Sphinx's reStructuredText (see `here <https://www.sphinx-doc.org/en/master/usage/restructuredtext/basics.html>`_). Files written in reStructuredText have a ``.rst`` file extension, and are then built as HTML filed by Sphinx during the documentation build. Only the RST file, not the built HTML file, are committed to the repository. Documentation should be clear and concise. See :doc:`../../io/visualization/using_widgets` as a good example of an RST-generated page.


IPYNB Documentation
-------------------

Often, code examples can help explain concepts better. The TARDIS utilizes `Jupyter notebooks <https://jupyter.org/>`_ (``.ipynb`` file extension) to demonstrate features of the code package within our documentation. See :doc:`../../quickstart` or :doc:`../../physics/montecarlo/initialization` for good examples.

TARDIS uses the `nbsphinx <https://nbsphinx.readthedocs.io/>`_ extension to turn these notebooks into HTML pages in the documentation. During a documentation build, nbsphinx runs all notebooks in the documentation with cleared output and places their output in the HTML. **Thus, notebook output must always be cleared before it is submitted** to ensure that the notebooks are run by nbsphinx. Running these notebooks during the documentation build helps ensure that the documentation is kept up-to-date, as notebook output will reflect the current state of the TARDIS code. Additionally, if updates in the code are inconsistent with the documentation, the documentation build will return an error, alerting TARDIS developers to the inconsistency.

An added benefit of IPYNB documentation is the ability to have interactive tutorials. All notebooks in the TARDIS documentation feature a button at the top encouraging users to launch the interactive version of the notebook (see the previously mentioned examples). This directs users to the TARDIS repository on `Binder <https://mybinder.org/>`_, where the notebook can be run using an online Jupyter kernel. Additionally, all notebooks in the Input/Output section of the documentation are automatically linked to on the :doc:`../../tutorials` page.


Including Your Page in the TARDIS Documentation
-----------------------------------------------

Whether your page is written in reStructuredText or as a Jupyter notebook, it must be included in the TARDIS documentation. This has three steps:

1. Determine the appropriate location for the page within the documentation. Feel free to reach out to someone in the TARDIS collaboration for help with this step.
2. Place your file in the corresponding directory in the ``docs/`` directory of the repository. For example, the :doc:`../../io/visualization/using_widgets` is a subpage of "Visualization Tools and Widgets" under the Input/Output section of the documentation, so it is placed in ``docs/io/visualization/``.
3. Include your file in the/a `toctree <https://www.sphinx-doc.org/en/master/usage/restructuredtext/directives.html#directive-toctree>`_ of the corresponding ``index.rst``. For example, :doc:`../../io/visualization/using_widgets` was included in a toctree of ``docs/io/visualization/index.rst``.

.. note::
    
    When new functions or classes are added to the code, in addition to documentation, **docstrings** must always be added. Read :ref:`this section <docstrings>` of our code quality guidelines to understand their importance and how they should be formatted. Sphinx uses these docstrings to auto-generate the `API documentation <https://tardis-sn.github.io/tardis/api/modules.html>`_ for entire TARDIS package. Please make sure that you have correctly formatted the docstrings by checking how the corresponding module's API looks once you build the documentation.


Building documentation locally
==============================

To build TARDIS documentation locally, use the following commands:

.. code::

    cd docs
    make html

.. note:: 

    - If you're working on a fresh local copy of the TARDIS repository, you might need to do ``python setup.py develop`` before executing these commands.
    - Use ``DISABLE_NBSPHINX=1 make html`` to disable notebook rendering (fast mode).
    - Use ``make html CORES=<number of cores>`` to have the documentation build in parallel. Using ``make html CORES=auto`` instructs Sphinx to use all of your device's cores.

After running this command, you can find the built docs (i.e. HTML webpages) in ``docs/_build/html``. Open the ``index.html`` in your browser to see how the documentation looks like with your edits. Navigate to page where you made changes or file that you added to check whether it looks as intended or not.

Additionally, check your terminal for warning messages during the documentation build (often caused by faulty hyperlinks or failing to include the page in the documentation). These should be repaired prior to merging your changes into the documentation.


.. _doc-preview:

Sharing the built documentation in your PR (Documentation Preview)
==================================================================

When you make edits in TARDIS documentation and submit a PR, we can only see the changes in source files in GitHub files diff, but not the built documentation (webpages). This is usually fine unless you have made changes in the documentation itself or changes that could break the Jupyter notebooks used in the documentation. In such cases, you should share the preview of documentation with your changes by building it online via GitHub. This will help us (the reviewers of your PR) to check how the documentation will look once your PR is merged.

To preview your changes to the documentation on GitHub, please:

#. Enable GitHub Actions in the *Actions* tab of your fork.
#. Under *Settings -> Pages* in your fork, make sure GitHub Pages is being built from the ``gh-pages`` branch and the ``/ (root)`` folder.

Then, there are two ways to trigger the build:

#. If the branch you are working on contains the word ``doc`` in it, then every commit pushed to that branch will trigger the build.
#. If your commit message contains the ``[build docs]`` tag, then that commit will trigger the build.

You can check for warning messages via the *Actions* tab of your fork. Your preview will be available at ``<username>.github.io/tardis/branch/<branch name>/index.html``.

.. note::

    You always can trigger a new build by pushing an empty commit: ``git commit --allow-empty -m "[build docs]"``

.. warning::
    
    On GitHub, built documentation files (including ``.ipynb`` files built by Sphinx) can be a maximum of 100 MB. You can check the file sizes after a local documentation build in ``docs/_build/html``, or after a documentation preview on GitHub in the ``gh-pages`` branch. Note that image output in notebooks built by Sphinx are by default in SVG format. For detailed images, these images can be very large. If file size becomes a problem, you will need to change the image format for that notebook by placing ``%config InlineBackend.figure_formats='png2x'`` in a `hidden cell <https://nbsphinx.readthedocs.io/en/0.8.7/hidden-cells.html>`_ at the beginning of the notebook.