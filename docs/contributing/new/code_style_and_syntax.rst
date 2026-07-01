.. _new-code_style_and_syntax:

Code Style And Syntax
=====================

Guidance for writing readable, maintainable TARDIS code.

.. _new-code-quality:

Code Quality
------------

.. _new-reference-code-quality-reference:

Quick Reference
~~~~~~~~~~~~~~~

TARDIS follows PEP 8.

Rule locations:

- Permanent Ruff rules: ``pyproject.toml``
- Non-permanent Ruff rules: ``.ruff.toml``

Naming conventions:

- Functions: lowercase with underscores.
- Variables: lowercase with underscores.
- Classes: CapWords.

Pre-commit hooks are no longer used in the TARDIS developer workflow. Run Ruff
and tests directly from the active TARDIS development environment.

.. _new-explanation-code-quality:

Explanation: Code Quality
~~~~~~~~~~~~~~~~~~~~~~~~~

Code quality helps new developers understand previous work. In open-source
software such as TARDIS, high-quality code is essential.

High-quality code:

- Does what it is supposed to do.
- Does not contain defects or problems.
- If anything unexpected happens, crash informatively.
- Raises or handles exceptions deliberately.
- Is easy to read, maintain, and extend.

TARDIS follows PEP 8. Ruff handles many style concerns such as whitespace,
string quotes, and code layout, but developers are still responsible for naming
conventions:

- Function names should be lowercase with words separated by underscores.
- Variable names use the same convention as function names.
- Class names should use CapWords.

.. raw:: html

   <span style="color:red">Added: concrete TARDIS naming example.</span>


For example, ``tardis/io/model/parse_density_configuration.py`` uses descriptive
function names such as ``parse_density_section_config``,
``parse_density_from_csvy``, and ``calculate_power_law_density``. These names are
preferable to short generic names because they state the model component and
operation being performed.

Edge Cases and Exceptions
^^^^^^^^^^^^^^^^^^^^^^^^^

Code should anticipate errors that may occur during execution. If an exception
is likely and can be handled, handle it. If an edge case would cause incorrect
behavior, raise an appropriate exception with a useful message.

Example:

.. code-block:: python

   def _calculate_plotting_data(self, packets_mode, packet_wvl_range, distance):
       if packets_mode not in ["virtual", "real"]:
           raise ValueError(
               "Invalid value passed to packets_mode. Only "
               "allowed values are 'virtual' or 'real'"
           )
       # Rest of the code ...


Here, ``packets_mode`` must be either ``"virtual"`` or ``"real"``. A specific
``ValueError`` tells the user what went wrong and prevents invalid code from
continuing. This kind of guard is useful when an invalid value is likely enough
to occur and easy enough for the user to fix. TARDIS should not add exception
branches for every theoretical edge case; when a failure is unlikely, cannot be
handled locally, or would make the code misleadingly verbose, it can be better
to let the underlying error surface naturally.

.. _new-naming-conventions:

Naming Conventions
------------------

.. _new-explanation-variable-and-function-names:

Explanation: Variable and Function Names
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Use names that describe the scientific or codebase concept, not only the local
implementation detail. Prefer names such as ``density_configuration``,
``time_explosion``, ``velocity_field_index``, and ``density_field_index`` when those
objects refer to specific model quantities or schema locations.

Avoid names that are too generic for review, such as ``data``, ``arr``, ``x``, or
``idx``, unless the scope is very small and the meaning is obvious from the
surrounding line.

In TARDIS code today, ``idx`` or an ``_idx`` suffix is used for short-lived
integer positions in a local loop, NumPy array, or packet collection. The
indexed object should be named nearby, and the value should not escape the small
local scope. Use ``index`` or an ``_index`` suffix when the value has durable
meaning: a named schema field, a pandas/DataFrame index-like object, a stable
lookup location, or a variable that will be read by reviewers outside the
immediate line. If a variable stores the position of a named field, use the
``_index`` suffix and include the concept being indexed:

.. code-block:: python

   velocity_field_index = [
       field.name for field in csvy_model_config.datatype.fields
   ].index("velocity")

   density_field_index = [
       field.name for field in csvy_model_config.datatype.fields
   ].index("density")


Do not use bare ``idx`` for schema fields, table columns, pandas indexes, or
other values whose meaning matters beyond a tiny loop.

Good examples:

.. code-block:: python

   velocity_field_index = fields.index("velocity")
   density_field_index = fields.index("density")
   packet_idx = packet_indices[i]


Poor examples:

.. code-block:: python

   idx = fields.index("density")
   i = packet_indices[i]
   data = parse_density_section_config(...)


.. raw:: html

   <span style="color:red">Added: naming convention guidance for ``_idx``, ``_index``, and TARDIS-style variable/function names.</span>

.. _new-typed-numpy-arrays:

Typed NumPy Arrays
------------------

.. _new-reference-numpy-array-typing:

Reference: NumPy Array Typing
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Use ``numpy.typing.NDArray`` for typed NumPy array annotations in new or touched
code when the function accepts or returns concrete NumPy arrays. NumPy provides
``numpy.typing.NDArray`` as a runtime-available alias for annotating arrays with a
dtype and unspecified shape.

.. code-block:: python

   from typing import Any

   import numpy as np
   import numpy.typing as npt


   def normalize(values: npt.NDArray[np.float64]) -> npt.NDArray[np.float64]:
       return values / np.sum(values)


   def as_array(values: npt.ArrayLike) -> npt.NDArray[Any]:
       return np.asarray(values)


Prefer typed arrays over bare ``np.ndarray`` when the dtype is meaningful to the
reader or reviewer:

.. code-block:: python

   def calculate_density_after_time(
       density_0: npt.NDArray[np.float64],
       time_0: float,
       time_explosion: float,
   ) -> npt.NDArray[np.float64]:
       ...


The typed NumPy API is stricter than runtime NumPy. For example, type checkers
will discourage patterns that create object arrays accidentally or mutate array
dtypes directly. If a function accepts broad array-like input, use
``npt.ArrayLike`` for the input and return a typed ``npt.NDArray[...]`` after
conversion.

See the official NumPy typing reference:
https://numpy.org/doc/stable/reference/typing.html

.. raw:: html

   <span style="color:red">Added: typed NumPy ndarray guidance and examples.</span>
