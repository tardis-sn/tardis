.. _new-documentation:

Documentation
=============

Guidance for writing, building, previewing, and troubleshooting developer documentation.

.. _new-docstrings:


Docstrings
----------

.. _new-reference-docstring-reference:

Reference: Docstring Reference
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A docstring describes a module, function, class, or method definition. It is
stored as ``object.__doc__`` and is surrounded by triple double quotes.

TARDIS follows the NumPy docstring format:

https://numpydoc.readthedocs.io/en/latest/format.html

Sphinx uses docstrings to auto-generate API documentation:

https://tardis-sn.github.io/tardis/api/modules.html

Example adapted from ``tardis/io/model/parse_density_configuration.py``:

.. code-block:: python

   def parse_density_section_config(
       density_configuration: ConfigurationNameSpace,
       v_middle: u.Quantity,
       time_explosion: u.Quantity,
   ) -> tuple[u.Quantity, u.Quantity]:
       """
       Parse the density section of the configuration file.

       Parameters
       ----------
       density_configuration : ConfigurationNameSpace
           Density configuration namespace.
       v_middle : astropy.units.Quantity
           Middle of the velocity bins.
       time_explosion : astropy.units.Quantity
           Time of the explosion.

       Returns
       -------
       density_0 : astropy.units.Quantity
           Density at time_0.
       time_0 : astropy.units.Quantity
           Time of the density profile.
       """


Docstring conventions:

- Do not include leading or trailing carriage returns.
- Include a carriage return between each segment.
- Start with a summary of the function, class, module, or method.
- The summary should use standard English syntax, start with a capital letter,
  and end with appropriate punctuation.
- The summary should describe purpose, not individual lines or return values.
- Comments on individual lines should be inline comments.
- Variable, module, function, and class names should be written between single
  backticks in prose.
- In the ``Returns`` section, always state the type, even if the variable name is
  omitted.
- Do not include a ``Returns`` section when there is no return value.
- Always list the full path for a variable type if it is not a built-in type,
  such as ``astropy.units.Quantity``.

Returns section format:

.. code-block:: python

   """
   Returns
   -------
   (`optional variable name` : )type
       (optional descriptor)
   """

.. _new-documentation-builds:

Documentation Builds
--------------------

.. _new-how-to-guide-build-documentation-locally:

How-To Guide: Build Documentation Locally
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Build the documentation:

.. code-block:: shell

   cd docs
   make html


Notes:

- On a fresh local copy, you may need to run ``pip install -e .`` first.
- Use ``DISABLE_NBSPHINX=1 make html`` to disable notebook rendering for a faster
  build.
- Use ``make html NCORES=<number of cores>`` to build in parallel.
- Use ``make html NCORES=auto`` to use all available device cores.
- Use ``make html SPHINXOPTS="<insert sphinx options>"`` to pass additional
  Sphinx options.

After the build, open ``docs/_build/html/index.html`` in your browser. Navigate to
the changed or added page and check that it looks as intended. Check the
terminal for warnings, often caused by faulty hyperlinks or missing toctree
entries. Fix warnings before merging.

When checking only RST changes or link structure, use:

.. code-block:: shell

   cd docs
   DISABLE_NBSPHINX=1 make html NCORES=auto


.. _new-reference-documentation-command-reference:

Reference: Documentation Command Reference
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Build documentation:

.. code-block:: shell

   cd docs
   make html


Build without notebook rendering:

.. code-block:: shell

   DISABLE_NBSPHINX=1 make html


Build in parallel:

.. code-block:: shell

   make html NCORES=<number of cores>
   make html NCORES=auto


Pass Sphinx options:

.. code-block:: shell

   make html SPHINXOPTS="<insert sphinx options>"


Built docs location:

.. code-block:: text

   docs/_build/html/index.html


Pull request documentation preview URL:

.. code-block:: text

   https://tardis-sn.github.io/tardis/pull/<pull request number>/index.html


Notebook file-size workaround:

.. code-block:: python

   %config InlineBackend.figure_formats='png2x'


.. _new-how-to-guide-share-built-documentation-in-a-pull-request:

How-To Guide: Share Built Documentation in a Pull Request
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Add the ``build-docs`` label to your pull request. If you cannot add the label,
leave a comment in the pull request or contact a senior member of the
collaboration.

The documentation builds when the label is added. Subsequent commits trigger new
documentation builds while the label remains present.

The built documentation is available at:

.. code-block:: text

   https://tardis-sn.github.io/tardis/pull/<pull request number>/index.html


It is also linked automatically in pull request comments.

To view build logs, go to the Actions tab in the TARDIS repository and select
the ``docs`` workflow. Search documentation builds by branch to find the relevant
log.

.. _new-how-to-guide-troubleshoot-documentation-builds:

How-To Guide: Troubleshoot Documentation Builds
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Documentation should be free of warnings and errors.

Common causes:

- Errors often mean notebooks are incompatible with new code. Update notebooks
  to reflect the code changes.
- Warnings often come from incorrect RST syntax in links, section headers,
  tables of contents, or similar structures.
- Warnings can also come from docstrings that do not follow the NumPy docstring
  format.
- GitHub built documentation files, including ``.ipynb`` files built by Sphinx,
  can be at most 100 MB.
- Check file sizes after a local documentation build in ``docs/_build/html``.
- Notebook image output built by Sphinx defaults to SVG. Detailed SVG images can
  be very large.
- If file size becomes a problem, add
  ``%config InlineBackend.figure_formats='png2x'`` in a hidden cell at the
  beginning of the notebook.

The Sublime and Sphinx Guide is a useful resource for RST syntax:

https://sublime-and-sphinx-guide.readthedocs.io/en/latest/index.html

Reach out for help if documentation issues are difficult to resolve.

.. _new-documentation-pages:

Documentation Pages
-------------------

.. _new-reference-developer-documentation-map:

Reference: Developer Documentation Map
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

General developer workflow pages:

- Reporting issues
- Git workflow
- Documentation guidelines
- Running tests
- Benchmarks
- Code quality
- Developer FAQ

Advanced core team development pages:

- Continuous integration
- Updating regression data
- Matterbridge
- Debugging ``numba_montecarlo``

.. _new-how-to-guide-document-code-changes:

How-To Guide: Document Code Changes
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

When you make or add functionality in TARDIS, create or update an ``.rst`` file or
Jupyter notebook (``.ipynb``) to demonstrate how the feature works. Include the
page in the documentation.

For RST documentation:

1. Use reStructuredText for pages that do not feature interactive code examples.
2. Commit the ``.rst`` source file, not built HTML.
3. Keep the documentation clear and concise.

For notebook documentation:

1. Use Jupyter notebooks when code examples help explain concepts.
2. TARDIS uses ``nbsphinx`` to turn notebooks into HTML pages.
3. During documentation builds, ``nbsphinx`` runs notebooks with cleared output and
   places generated output in the HTML.
4. Always clear notebook output before submitting notebooks.
5. In VS Code, use the "Clear All Outputs" command.
6. In JupyterLab, use ``Edit > Clear Outputs of All Cells``.

Running notebooks during the documentation build helps keep documentation
up to date. If code updates are inconsistent with documentation, the
documentation build returns an error.

TARDIS notebook documentation can provide interactive tutorials. Documentation
notebooks include a launch button that directs readers to Binder, where the
notebook can run with an online Jupyter kernel. Notebooks in the Input/Output
section are automatically linked from the tutorials page.

When new functions or classes are added, add docstrings as well as page-level
documentation. Sphinx uses docstrings to auto-generate the API documentation for
the TARDIS package. Build the documentation and check how the corresponding
module API looks.

.. _new-explanation-documentation-system:

Explanation: Documentation System
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

High-quality and consistent documentation helps users find specific tasks and
helps developers understand best practices.

TARDIS uses Sphinx to generate documentation. Sphinx translates source files,
often reStructuredText or Jupyter notebooks, into HTML files and automatically
produces cross-references and indices. Developers new to Sphinx should read the
Sphinx quickstart guide.

RST documentation is used for pages without interactive code examples. Notebook
documentation is used when code examples help explain concepts. TARDIS uses
``nbsphinx`` to render notebooks into HTML.

Only source documentation files are committed. Built HTML is not committed.

Notebook output must always be cleared before submission because ``nbsphinx`` runs
notebooks during the build and inserts fresh output into the generated HTML.
This helps keep documentation synchronized with the current codebase.

.. _new-how-to-guide-include-a-new-documentation-page:

How-To Guide: Include a New Documentation Page
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Whether the page is reStructuredText or a Jupyter notebook:

1. Determine the appropriate location in the documentation. Ask the TARDIS
   collaboration for help if needed.
2. Place the file in the corresponding directory under ``docs/``.
3. Include the file in a ``toctree`` in the corresponding ``index.rst``.

Example: a page under "Visualization Tools and Widgets" in the Input/Output
section belongs in a corresponding directory under ``docs/`` and must be included
in that section's ``index.rst``.

For example, a new developer page at
``docs/contributing/development/new_tooling.rst`` would be linked from
``docs/contributing/development/index.rst`` like this:

.. code-block:: rst

   .. toctree::
       :maxdepth: 2

       issues
       git_workflow
       documentation_guidelines
       new_tooling
