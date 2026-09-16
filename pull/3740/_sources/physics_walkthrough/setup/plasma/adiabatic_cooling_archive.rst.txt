.. _adiabatic_cooling_archive:

Archived adiabatic-cooling calculation
--------------------------------------

This page preserves the adiabatic electron-gas cooling calculation removed
with the legacy continuum plasma properties. It is documentation only: the
classic plasma workflow no longer assembles or evaluates this property.

Cooling rate
~~~~~~~~~~~~

For each model shell, the legacy property used the electron number density
:math:`n_e`, electron temperature :math:`T_e`, Boltzmann constant :math:`k_B`,
and time since explosion :math:`t_\mathrm{explosion}` to calculate

.. math::

    C_\mathrm{adiabatic} =
        \frac{3 n_e k_B T_e}{t_\mathrm{explosion}}.

The original implementation, with ``K_B`` defined as the CGS value of the
Boltzmann constant, was:

.. code-block:: python

    cool_rate_adiabatic = (
        3.0 * electron_densities * K_B * t_electrons
    ) / time_explosion

Transition-probability representation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The resulting per-shell cooling rates were represented as a single
macro-atom deactivation row. Its MultiIndex was
``("k", "adiabatic", -1)`` with levels named ``source_level_idx``,
``destination_level_idx``, and ``transition_type``. The generic helper used
by the deleted property was:

.. code-block:: python

    def cooling_rate_series2dataframe(cooling_rate_series, destination_level_idx):
        index_names = [
            "source_level_idx",
            "destination_level_idx",
            "transition_type",
        ]
        index = pd.MultiIndex.from_tuples(
            [("k", destination_level_idx, -1)], names=index_names
        )
        return pd.DataFrame(cooling_rate_series.values[np.newaxis], index=index)

The property passed ``destination_level_idx="adiabatic"`` to this helper:

.. code-block:: python

    cool_rate_adiabatic = cooling_rate_series2dataframe(
        cool_rate_adiabatic, destination_level_idx="adiabatic"
    )

This representation and calculation are retained for scientific traceability
and are not an available classic-plasma API.
