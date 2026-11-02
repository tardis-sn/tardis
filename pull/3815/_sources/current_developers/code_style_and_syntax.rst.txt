.. _code_style_and_syntax:

Code Style And Syntax
=====================

Guidance for writing readable, maintainable TARDIS code.

.. _code-quality:

Code Quality
------------

.. _reference-code-quality-reference:

Quick Reference
~~~~~~~~~~~~~~~

TARDIS follows PEP 8 using Ruff. Run Ruff directly from the active TARDIS development environment.

Rule locations:

- Permanent Ruff rules: ``pyproject.toml``
- Non-permanent Ruff rules: ``.ruff.toml``

Naming conventions: https://peps.python.org/pep-0008/#naming-conventions

.. _naming-conventions:

Naming Conventions
------------------

We prefer names to be descriptive. For example, ``tardis/io/model/parse_density_configuration.py`` uses descriptive
function names such as ``parse_density_section_config``,
``parse_density_from_csvy``, and ``calculate_power_law_density``. These names are
preferable to short generic names because they state the operation being performed 
and the model component being operated on.

While functions and methods should use a verb-noun style, class properties 
and variables should use a noun style. For example, 
``calculate_power_law_density`` is a good function name, 
while ``power_law_density`` is a good variable or property name.

.. _explanation-variable-and-function-names:

Explanation: Variable and Function Names
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Use names that describe the scientific or codebase concept, not only the local
implementation detail. Prefer names such as ``density_configuration``,
``time_explosion``, ``velocity_field_index``, and ``density_field_index`` when those
objects refer to specific model quantities or schema locations. Ensure that if in
other parts of the codebase similar concepts exist, use those names.

If you use an object property e.g. ``packet.velocity``, it is reasonable to
assign it to ``velocity`` in a local scope if it is used more than once or twice.

Avoid names that are too generic for review, such as ``data``, ``arr``, ``x``, or
``idx``, unless the scope is very small and the meaning is obvious from the
surrounding line.

.. _index_vs_idx:

Index vs. IDX in TARDIS
~~~~~~~~~~~~~

In the TARDIS codebase, we attempt to use index and idx as separate designations 
with different meanings. We use {var_name}_idx to specify that a variable refers to the 
integer position in some array-like object, with pythonic zero-indexed location 
(i.e., any pandas object accessed via .iloc would be done so with an idx 
variable). Alternatively, we use {var_name}_index as a variable that is explicitly not 
an integer location, but is still used for lookup via a key to a dictionary-like 
object (e.g., a pandas object accessed with .loc). 

Consequently, an object named index in the tardis codebase will often store idx
values. The dictionary-like is the index, and the values are idxs. 

Good examples:

.. code-block:: python

    block_start_idx = opacity_state.macro_block_edge_index[
        absorbing_activation_level_idx
    ]
    block_end_idx = opacity_state.macro_block_edge_index[
        absorbing_activation_level_idx + 1
    ]
    emission_transition_probability = 0.0
    probability_emission_event = np.random.random()

    for deactivation_channel_idx in range(block_start_idx, block_end_idx):
        deactivation_probability = opacity_state.transition_probabilities[
            deactivation_channel_idx, current_shell_id
        ]
        emission_transition_probability += deactivation_probability

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
