.. _new-data_and_codebase_structure:

Data And Codebase Structure
===========================

Guidance for local development setup, codebase-specific design notes, and debugging support.

.. _new-developer-faq:

Developer FAQ
-------------

.. _new-reference-developer-faq:

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

.. _new-numba-debugging:

Numba Debugging
---------------

.. _new-how-to-guide-debug-numbamontecarlo:

How-To Guide: Debug ``numba_montecarlo``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This section summarizes practical debugging workflows for TARDIS's
Numba-powered Monte Carlo path. It is based on Numba's current debugging
guidance:

- https://numba.readthedocs.io/en/stable/developer/debugging.html
- https://numba.readthedocs.io/en/stable/user/troubleshoot.html#debugging-jit-compiled-code-with-gdb

Notebook Debugging Workflow
^^^^^^^^^^^^^^^^^^^^^^^^^^^

For most debugging tasks, use a Jupyter notebook so you can inspect state
interactively between steps.

Suggested structure:

1. Cell 1: imports and environment flags.
2. Cell 2: load config and inputs.
3. Cell 3 and later: run short Monte Carlo steps and inspect intermediate
   objects.

Start with small, reproducible runs with fewer packets or iterations so each
cell executes quickly.

Recommended Debugging Order
^^^^^^^^^^^^^^^^^^^^^^^^^^^

When investigating a bug in the Numba Monte Carlo path, start with the fastest
loop and escalate only as needed:

1. Disable JIT to separate Python-level logic errors from compiled-code issues.
2. Re-enable JIT and use Numba debug information with low optimization for
   step-through debugging.
3. Use GDB or Memcheck for hard crashes or suspected native memory corruption.

Set Environment Variables Before Imports
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Set debug-related environment variables before importing Numba or TARDIS modules
that trigger Numba compilation:

.. code-block:: python

   import os

   # Toggle this to "1" for Python-level debugging, "0" for normal JIT behavior.
   os.environ["NUMBA_DISABLE_JIT"] = "1"

   # Useful when JIT is enabled and you need better debugger visibility.
   os.environ["NUMBA_DEBUGINFO"] = "1"
   os.environ["NUMBA_OPT"] = "0"
   os.environ["NUMBA_EXTEND_VARIABLE_LIFETIMES"] = "1"

After changing these values, restart the notebook kernel and rerun from the
first cell.

Debug Python-Level Logic With JIT Disabled
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Numba documents ``NUMBA_DISABLE_JIT=1`` as the first-line option for debugging
logic problems. This makes ``@jit`` and ``@njit`` functions execute as regular
Python functions.

In a notebook cell:

.. code-block:: python

   import os

   os.environ["NUMBA_DISABLE_JIT"] = "1"

   from tardis import run_tardis

   sim = run_tardis("path/to/your_config.yml")

This is usually the best way to use standard Python debugger breakpoints, get
clearer Python exceptions, and determine whether a failure is specific to JIT
compilation.

Debug Numba Type Errors
^^^^^^^^^^^^^^^^^^^^^^^

Type errors are common while developing Numba code. When a Numba typing error
appears:

1. Read the first typing-error block that mentions the failing function.
2. Reduce the failing call to the smallest input that still reproduces the
   error.
3. Check array dtypes, shapes, optional values, units, and object-mode
   operations that Numba cannot compile.
4. Run with ``NUMBA_DISABLE_JIT=1`` to separate ordinary Python logic problems
   from Numba typing or compilation constraints.

If the code works with JIT disabled but fails when JIT is enabled, focus on
types and operations at the compiled function boundary.

Debug Compiled CPU Code With GDB
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For issues that only appear with JIT enabled, run with debug-oriented Numba
settings:

.. code-block:: shell

   NUMBA_DEBUGINFO=1 NUMBA_OPT=0 NUMBA_EXTEND_VARIABLE_LIFETIMES=1 \
       gdb -q --args python -c "from tardis import run_tardis; run_tardis('path/to/your_config.yml')"

Then in ``gdb``:

.. code-block:: text

   run
   bt
   info args
   info locals

Notes from Numba guidance that are especially relevant:

- ``NUMBA_DEBUGINFO=1`` enables debug symbols for jitted functions.
- ``NUMBA_OPT=0`` improves stepping and reduces optimized-out variables.
- ``NUMBA_EXTEND_VARIABLE_LIFETIMES=1`` makes local variable inspection more
  predictable.
- Debug information increases memory use and can slow compilation and
  execution.

For notebook users, once you suspect a native-level issue, reproduce it in a
small script or one-line command and attach GDB there. This is usually more
reliable than trying to drive ``gdb`` from inside an active notebook session.

Use Memcheck For Suspected Native Memory Errors
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If you suspect out-of-bounds access or use-after-free behavior in compiled code,
Numba recommends Valgrind Memcheck with suppression files from matching Python
and Numba versions.

.. code-block:: shell

   valgrind --tool=memcheck \
       --suppressions=${CPYTHON_SRC_DIR}/Misc/valgrind-python.supp \
       --suppressions=${NUMBA_SRC_DIR}/numba/misc/valgrind-numba.supp \
       python -c "from tardis import run_tardis; run_tardis('path/to/your_config.yml')"

Use this workflow only after reducing the case as much as possible; it is slower
and noisier than Python-level debugging.
