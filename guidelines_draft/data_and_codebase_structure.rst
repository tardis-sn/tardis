.. _data_and_codebase_structure:

Data And Codebase Structure
===========================

Guidance for local development setup, codebase-specific design notes, and debugging support.

.. _developer-faq:

Developer FAQ
-------------

.. _reference-developer-faq:

Reference: Developer FAQ
~~~~~~~~~~~~~~~~~~~~~~~~

Constants in TARDIS are taken from Astropy. The ``tardis.constants`` module
imports all constants from ``astropy.constants.astropy13constants``.

Class design and inheritance:

- If only the constructor changed, use a classmethod.
- If overriding other methods, use a subclass.

TARDIS uses Ruff to check PEP 8 compliance.

.. raw:: html

   <span style="color:red">Added: examples for classmethod versus subclass guidance.</span>


Use a classmethod when construction changes but behavior stays the same, such as
adding ``from_config(...)`` or ``from_csvy(...)`` constructors around an existing
class. Use a subclass when method behavior changes after construction, such as a
specialized solver or reader that overrides parsing or calculation methods.

.. _development-install:

Development Install
-------------------

.. _how-to-guide-install-tardis-in-development-mode:

How-To Guide: Install TARDIS in Development Mode
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

From the repository root, run:

.. code-block:: shell

   pip install -e .


TARDIS is designed to be usable from the source tree. An editable install makes
``tardis`` import from your clone, so local edits are available immediately in new
Python sessions.

.. raw:: html

   <span style="color:red">Added: import-check example after editable install.</span>


After installation, confirm that Python imports from your checkout:

.. code-block:: shell

   python -c "import tardis; print(tardis.__file__)"

.. _numba-debugging:

Numba Debugging
---------------

.. _how-to-guide-debug-numbamontecarlo:

How-To Guide: Debug ``numba_montecarlo``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

TARDIS provides experimental debugging configurations for deeper debugging when
interfacing with the ``montecarlo_numba`` module.

PyCharm debugging configurations, related scripts, and ``.yml`` files are in
``tardis.scripts.debug``. They currently include single-packet mode.

Relevant files:

- ``tardis_example_single.yml``: configuration file for the single-packet TARDIS
  run.
- ``run_numba_single.py``: Python script that runs the ``.yml`` file.
- ``run_numba_single.xml``: PyCharm debug configuration file used with the files
  above.

This method is experimental.
