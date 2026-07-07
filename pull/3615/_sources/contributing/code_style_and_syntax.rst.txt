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

TARDIS follows PEP 8 using Ruff. Run Ruff directly from the active TARDIS development environment.

Rule locations:

- Permanent Ruff rules: ``pyproject.toml``
- Non-permanent Ruff rules: ``.ruff.toml``

Naming conventions: https://peps.python.org/pep-0008/#naming-conventions

.. _new-naming-conventions:

Naming Conventions
------------------

We prefer names to be descriptive. For example, ``tardis/io/model/parse_density_configuration.py`` uses descriptive
function names such as ``parse_density_section_config``,
``parse_density_from_csvy``, and ``calculate_power_law_density``. These names are
preferable to short generic names because they state the operation being performed 
and the model component being operated on.

.. _new-explanation-variable-and-function-names:

Explanation: Variable and Function Names
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Use names that describe the scientific or codebase concept, not only the local
implementation detail. Prefer names such as ``density_configuration``,
``time_explosion``, ``velocity_field_index``, and ``density_field_index`` when those
objects refer to specific model quantities or schema locations. Ensure that if in
other parts of the codebase similar concepts exist, use those names.

Avoid names that are too generic for review, such as ``data``, ``arr``, ``x``, or
``idx``, unless the scope is very small and the meaning is obvious from the
surrounding line.

In TARDIS, the ``idx`` or an ``_idx`` suffix should be used for short-lived
integer positions in a local loop, NumPy array, or packet collection. The
indexed object should be named nearby, and the value should not escape the small
local scope. Use ``index`` or an ``_index`` suffix when the value is used outside the scope: 
a named schema field, a pandas/DataFrame index-like object, a
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

.. _new-typed-numpy-arrays:

Type Hinting
------------

Use type hints for all function definitions in new or touched code. Type hints 
improve readability and maintainability, and they enable static type checking.

.. code-block:: python
    from tardis.opacities.opacity_state_numba import (
        OpacityStateNumba,
    )

    def macro_atom_interaction(
        activation_level_id: int,
        current_shell_id: int,
        opacity_state: OpacityStateNumba,
    ):
        pass

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


Prefer typed arrays over bare ``np.ndarray``.

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
