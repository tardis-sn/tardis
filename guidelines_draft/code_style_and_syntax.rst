.. _code_style_and_syntax:

Code Style And Syntax
=====================

Guidance for writing readable, maintainable TARDIS code.

.. _code-quality:

Code Quality
------------

.. _explanation-code-quality:

Explanation: Code Quality
~~~~~~~~~~~~~~~~~~~~~~~~~

Code quality helps new developers understand previous work. In open-source
software such as TARDIS, high-quality code is essential.

High-quality code:

- Does what it is supposed to do.
- Does not contain defects or problems.
- Handles edge cases safely.
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
continuing.

.. _reference-code-quality-reference:

Reference: Code Quality Reference
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

TARDIS follows PEP 8.

Rule locations:

- Permanent Ruff rules: ``pyproject.toml``
- Non-permanent Ruff rules: ``.ruff.toml``

Naming conventions:

- Functions: lowercase with underscores.
- Variables: lowercase with underscores.
- Classes: CapWords.

Pre-commit install:

.. code-block:: shell

   pip install pre-commit
   pre-commit install


.. raw:: html

   <span style="color:red">Deleted: removed duplicate Ruff command examples from the reference because the Ruff how-to gives runnable commands.</span>

.. _naming-conventions:

Naming Conventions
------------------

.. _explanation-variable-and-function-names:

Explanation: Variable and Function Names
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Use names that describe the scientific or codebase concept, not only the local
implementation detail. Prefer names such as ``density_configuration``,
``time_explosion``, ``velocity_field_index``, and ``density_field_index`` when those
objects refer to specific model quantities or schema locations.

Avoid names that are too generic for review, such as ``data``, ``arr``, ``x``, or
``idx``, unless the scope is very small and the meaning is obvious from the
surrounding line. If a variable stores the position of a named field, use the
``_index`` suffix and include the concept being indexed:

.. code-block:: python

   velocity_field_index = [
       field.name for field in csvy_model_config.datatype.fields
   ].index("velocity")

   density_field_index = [
       field.name for field in csvy_model_config.datatype.fields
   ].index("density")


Use ``_idx`` only for short-lived loop or array positions where the indexed object
is already named nearby. Use ``_index`` for persistent variables, pandas indexes,
schema indexes, or public-facing names where readability matters more than
brevity.

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

.. _pre-commit:

Pre-Commit
----------

.. _how-to-guide-use-pre-commit:

How-To Guide: Use Pre-Commit
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Pre-commit hooks are optional local tools that run checks before commits.

Install pre-commit:

.. code-block:: shell

   pip install pre-commit


Set up hooks for the repository:

.. code-block:: shell

   pre-commit install


This only needs to be done once per repository. The hooks then run automatically
on each commit.

.. raw:: html

   <span style="color:red">Added: example command for checking a TARDIS file before committing.</span>


To run the hooks manually on a specific file, use:

.. code-block:: shell

   pre-commit run --files tardis/io/model/parse_density_configuration.py

.. _ruff:


Ruff
----

.. _how-to-guide-run-ruff:

How-To Guide: Run Ruff
~~~~~~~~~~~~~~~~~~~~~~

TARDIS follows PEP 8 and uses Ruff for linting and formatting.

Install Ruff:

.. code-block:: shell

   conda install -c conda-forge ruff


Lint code:

.. code-block:: shell

   ruff check <source_file_or_directory>


Lint and fix automatically fixable issues:

.. code-block:: shell

   ruff check <source_file_or_directory> --fix


.. raw:: html

   <span style="color:red">Added: TARDIS-specific Ruff examples.</span>


For example, lint one parser module:

.. code-block:: shell

   ruff check tardis/io/model/parse_density_configuration.py


Or lint and fix a package area:

.. code-block:: shell

   ruff check tardis/io/model --fix


TARDIS adopts linting rules used by Astropy. Permanent rules are defined in
``pyproject.toml``. Non-permanent rules are defined in ``.ruff.toml``. Add new rules
to ``.ruff.toml``. To add a permanent rule, open a pull request against
``pyproject.toml``.

.. _typed-numpy-arrays:

Typed NumPy Arrays
------------------

.. _reference-numpy-array-typing:

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
