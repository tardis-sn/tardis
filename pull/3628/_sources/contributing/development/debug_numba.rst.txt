****************************
Debugging Numba-powered code
****************************

This page summarizes practical debugging workflows for TARDIS's Numba-powered
code. It is based on Numba's current debugging guidance:

- https://numba.readthedocs.io/en/stable/developer/debugging.html
- https://numba.readthedocs.io/en/stable/user/troubleshoot.html#debugging-jit-compiled-code-with-gdb

Notebook debugging workflow
======================================

For most debugging tasks, use a Jupyter notebook so you can inspect state
interactively between steps.

Suggested structure:

1. Cell 1: imports + environment flags.
2. Cell 2: load config / inputs.
3. Cell 3+: run short Monte Carlo steps and inspect intermediate objects.

Start with small, reproducible runs (fewer packets / iterations) so each cell
executes quickly and you can iterate rapidly.

Recommended debugging workflow
==============================

When investigating a bug in the Numba Monte Carlo path, start with the fastest
loop and only escalate if needed.

1. Disable JIT to separate Python-level logic errors from compiled-code issues.
2. Re-enable JIT and use Numba debug info + low optimization for step-through.
3. Use GDB or Memcheck only for hard crashes / suspected memory corruption.

Set environment variables before importing Numba/TARDIS
=======================================================

In notebooks, set debug-related environment variables in the first code cell,
before importing modules that trigger Numba compilation:

.. code-block:: python

	import os

	# Toggle this to "1" for Python-level debugging, "0" for normal JIT behavior.
	os.environ["NUMBA_DISABLE_JIT"] = "1"

	# Useful when JIT is enabled and you need better debugger visibility.
	os.environ["NUMBA_DEBUGINFO"] = "1"
	os.environ["NUMBA_OPT"] = "0"
	os.environ["NUMBA_EXTEND_VARIABLE_LIFETIMES"] = "1"

After changing these values, restart the notebook kernel and re-run from the
first cell.

1) Python-level debugging (disable JIT)
=======================================

Numba documents ``NUMBA_DISABLE_JIT=1`` as the first-line option for debugging
logic problems. This makes ``@jit``/``@njit`` functions execute as regular
Python functions.

In a notebook cell:

.. code-block:: python

	import os
	os.environ["NUMBA_DISABLE_JIT"] = "1"

	from tardis import run_tardis
	sim = run_tardis("path/to/your_config.yml")

This is usually the best way to:

- use standard Python debugger breakpoints,
- get clearer Python exceptions,
- verify whether a failure is specific to JIT compilation.

Notebook tips for live inspection:

- use short cells to isolate one step at a time,
- print/plot intermediate estimators instead of only final spectra,
- if state becomes unclear, restart kernel and re-run sequentially.

2) Diagnose and fix type errors
===============================

Many Numba failures are typing failures (for example ``TypingError`` and
"cannot unify" messages).

They can often be resolved by simply disabled JIT and allowing pure Python
errors to reveal the underlying type issues.

Detailed workflow for harder problems:

- isolate the failing call in a cell that executes only the problematic
	function where possible,
- read the failing operation and operand types directly from the exception,
- force stable dtypes at entry and exit points (for example with ``numpy.asarray(...,
	dtype=...)``),
- avoid mixed return types across branched code paths,
- replace empty or mixed Python lists with typed containers or fixed-dtype
	arrays.

When needed, temporarily split a large jitted function into smaller helpers so
the failing expression is easier to localize.

1) Debug compiled CPU code with GDB
===================================

For issues that only appear with JIT enabled, run with debug-oriented Numba
settings:

.. code-block:: bash

	NUMBA_DEBUGINFO=1 NUMBA_OPT=0 NUMBA_EXTEND_VARIABLE_LIFETIMES=1 \
	gdb -q --args python -c "from tardis import run_tardis; run_tardis('path/to/your_config.yml')"

Then in ``gdb``:

.. code-block:: gdb

	run
	bt
	info args
	info locals

Notes from Numba guidance that are especially relevant:

- ``NUMBA_DEBUGINFO=1`` enables debug symbols for jitted functions.
- ``NUMBA_OPT=0`` improves stepping and reduces "optimized out" variables.
- ``NUMBA_EXTEND_VARIABLE_LIFETIMES=1`` makes local variable inspection more
  predictable.
- Debug info increases memory usage and can slow compilation/execution.

For notebook users: once you suspect a native-level issue, reproduce it in a
small script or one-liner command and attach GDB there. This is typically more
reliable than trying to drive ``gdb`` from inside an active notebook session.

4) Memcheck for suspected native memory errors
==============================================

If you suspect out-of-bounds access or use-after-free behavior in compiled
code, Numba recommends running under Valgrind Memcheck with suppression files
from the matching Python and Numba versions.

.. code-block:: bash

	valgrind --tool=memcheck \
	  --suppressions=${CPYTHON_SRC_DIR}/Misc/valgrind-python.supp \
	  --suppressions=${NUMBA_SRC_DIR}/contrib/valgrind-numba.supp \
	  python -c "from tardis import run_tardis; run_tardis('path/to/your_config.yml')"

This is significantly slower than notebook iteration, so use it only after you
have a minimal reproducer.

Practical notebook caveats
==========================

- Environment variables must be set before Numba compilation; kernel restart is
  often required.
- Re-running cells out of order can hide state-dependent bugs.
- Keep one "known-good" top-to-bottom execution path for reproducibility.
