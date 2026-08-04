from dataclasses import dataclass
from itertools import pairwise

import numpy as np
import pandas as pd
from astropy import units as u

from tardis import constants as const
from tardis.plasma.electron_energy_distribution import (
    ThermalElectronEnergyDistribution,
)
from tardis.plasma.equilibrium.rates import (
    AnalyticCorrectedPhotoionizationCoeffSolver,
    CollisionalIonizationSeaton,
    EstimatedPhotoionizationCoeffSolver,
    SpontaneousRecombinationCoeffSolver,
    ThermalCollisionalRateSolver,
)
from tardis.plasma.equilibrium.rates.heating_cooling_rates import (
    BoundFreeThermalRates,
    FreeFreeThermalRates,
)


def _cumulative_integrate_by_blocks(
    values: np.ndarray, frequencies: np.ndarray, block_references: np.ndarray
) -> np.ndarray:
    """Cumulatively integrate frequency-dependent values within each block."""
    integrated = np.zeros_like(values)
    for start, stop in pairwise(block_references):
        if stop - start < 2:
            continue
        increments = (
            np.diff(frequencies[start:stop])[:, np.newaxis]
            * (values[start + 1 : stop] + values[start : stop - 1])
            / 2.0
        )
        integrated[start + 1 : stop] = np.cumsum(increments, axis=0)
        totals = integrated[stop - 1]
        nonzero = totals != 0.0
        integrated[start + 1 : stop, nonzero] /= totals[nonzero]
    return integrated


def _normalize_rows(rates: pd.DataFrame) -> tuple[pd.Series, np.ndarray]:
    """Return per-cell channel probabilities and contiguous probabilities."""
    total = rates.sum(axis=0)
    probabilities = rates.divide(total.replace(0.0, 1.0), axis=1)
    # Macro-atom cooling kernels consume arrays as (cell, transition), while
    # rate tables are conventionally stored as (transition, cell).
    return total.to_numpy(), probabilities.to_numpy().T


def _macro_atom_collisional_index(
    index: pd.Index, *, reverse: bool = False
) -> pd.Index:
    """Use the transition-level names expected by macro-atom kernels."""
    if reverse:
        index = index.swaplevel(
            "level_number_source", "level_number_destination"
        )
    redundant_levels = [
        level
        for level in ("ion_number_source", "ion_number_destination")
        if level in index.names
    ]
    if redundant_levels:
        index = index.droplevel(redundant_levels)
    level_names = (
        {
            "level_number_source": "level_number_upper",
            "level_number_destination": "level_number_lower",
        }
        if reverse
        else {
            "level_number_source": "level_number_lower",
            "level_number_destination": "level_number_upper",
        }
    )
    return index.rename(level_names)


@dataclass
class EquilibriumContinuumState:
    """Continuum coefficients consumed by the continuum macro-atom solver."""

    radiative_ionization_rate: pd.DataFrame
    radiative_recombination_rate: pd.DataFrame
    collisional_excitation_rate: pd.DataFrame
    collisional_deexcitation_rate: pd.DataFrame
    collisional_ionization_rate: pd.DataFrame
    collisional_recombination_rate: pd.DataFrame
    delta_E_yg: pd.Series  # noqa: N815 - retained continuum-state field name
    collisional_excitation_cooling_probability: pd.Series
    collisional_excitation_cooling_array: np.ndarray
    collisional_excitation_references: pd.Index
    collisional_ionization_cooling_probability: pd.Series
    collisional_ionization_cooling_array: np.ndarray
    radiative_recombination_cooling_probability: pd.Series
    radiative_recombination_cooling_array: np.ndarray
    free_free_cooling_probability: pd.Series

    @classmethod
    def from_plasma(
        cls,
        plasma: object,
        estimators: dict[str, object] | None = None,
    ) -> EquilibriumContinuumState:
        """Build continuum inputs from standard equilibrium-plasma outputs."""
        electron_distribution = ThermalElectronEnergyDistribution(
            0 * u.erg,
            plasma.t_electrons * u.K,
            plasma.electron_densities.values * u.cm**-3,
        )
        photo_data = plasma.atomic_data.photoionization_data
        if (
            estimators is not None
            and estimators["photoionization_rate_estimator"] is not None
        ):
            corrected_solver = EstimatedPhotoionizationCoeffSolver(
                plasma.atomic_data.level2continuum_edge_idx
            )
            radiative_ionization_rate = corrected_solver.solve(
                {
                    "photoionization_rate_estimator": estimators[
                        "photoionization_rate_estimator"
                    ],
                    "stimulated_recombination_rate_estimator": estimators[
                        "stimulated_recombination_rate_estimator"
                    ],
                },
                plasma.level_number_density,
                plasma.lte_level_number_density,
                plasma.ion_number_density,
                plasma.lte_ion_number_density,
            )
        else:
            corrected_solver = AnalyticCorrectedPhotoionizationCoeffSolver(
                photo_data
            )
            radiative_ionization_rate = corrected_solver.solve(
                plasma.dilute_planckian_radiation_field,
                electron_distribution.temperature,
                plasma.lte_level_number_density,
                plasma.level_number_density,
                plasma.lte_ion_number_density,
                plasma.ion_number_density,
            )
        spontaneous = SpontaneousRecombinationCoeffSolver(photo_data).solve(
            electron_distribution.temperature
        )
        lower_ions = spontaneous.index.droplevel("level_number")
        upper_ions = pd.MultiIndex.from_arrays(
            [
                lower_ions.get_level_values("atomic_number"),
                lower_ions.get_level_values("ion_number") + 1,
            ],
            names=["atomic_number", "ion_number"],
        )
        level_factor = plasma.lte_level_number_density.loc[
            spontaneous.index
        ].divide(
            plasma.lte_ion_number_density.loc[upper_ions].to_numpy()
            * electron_distribution.number_density.value
        )
        radiative_recombination_rate = spontaneous * level_factor

        lines = plasma.atomic_data.lines
        collisional_solver = ThermalCollisionalRateSolver(
            plasma.atomic_data.levels,
            lines,
            plasma.atomic_data.collision_data_temperatures,
            plasma.atomic_data.yg_data,
            "cmfgen",
        )
        collisional_rates = collisional_solver.solve(
            electron_distribution.temperature
        )
        source = collisional_rates.index.get_level_values("level_number_source")
        destination = collisional_rates.index.get_level_values(
            "level_number_destination"
        )
        excitation_mask = source < destination
        collisional_excitation_rate = collisional_rates[excitation_mask]
        collisional_deexcitation_rate = collisional_rates[~excitation_mask]

        collisional_ionization_rate = CollisionalIonizationSeaton(
            photo_data
        ).solve(electron_distribution.temperature)
        collisional_recombination_rate = collisional_ionization_rate.multiply(
            level_factor.loc[collisional_ionization_rate.index]
        )

        level_energy = plasma.atomic_data.levels.energy
        transition_index = collisional_excitation_rate.index
        lower_index = pd.MultiIndex.from_arrays(
            [
                transition_index.get_level_values("atomic_number"),
                transition_index.get_level_values("ion_number"),
                transition_index.get_level_values("level_number_source"),
            ],
            names=["atomic_number", "ion_number", "level_number"],
        )
        upper_index = pd.MultiIndex.from_arrays(
            [
                transition_index.get_level_values("atomic_number"),
                transition_index.get_level_values("ion_number"),
                transition_index.get_level_values("level_number_destination"),
            ],
            names=["atomic_number", "ion_number", "level_number"],
        )
        delta_E_yg = pd.Series(
            level_energy.loc[upper_index].values
            - level_energy.loc[lower_index].values,
            index=collisional_excitation_rate.index,
        )
        excitation_cooling = (
            collisional_excitation_rate
            * plasma.level_number_density.loc[lower_index].values
            * electron_distribution.number_density.value
            * delta_E_yg.values[:, np.newaxis]
        )
        ionization_energy = (
            photo_data.groupby(level=[0, 1, 2]).first().nu * const.h.cgs.value
        )
        ionization_cooling = (
            collisional_ionization_rate
            * plasma.level_number_density.loc[collisional_ionization_rate.index]
            * electron_distribution.number_density.value
            * ionization_energy.loc[collisional_ionization_rate.index].values[
                :,
                np.newaxis,
            ]
        )
        stimulated_cooling_estimator = (
            None
            if estimators is None
            else estimators.get("stimulated_recombination_cooling_estimator")
        )
        recombination_cooling = BoundFreeThermalRates(
            photo_data
        ).calculate_cooling_rate(
            plasma.ion_number_density,
            electron_distribution,
            level_factor,
            stimulated_cooling_estimator,
        )

        exc_total, exc_array = _normalize_rows(excitation_cooling)
        ion_total, ion_array = _normalize_rows(ionization_cooling)
        fb_total, fb_array = _normalize_rows(recombination_cooling)
        free_free_solver = FreeFreeThermalRates()
        free_free_heating_factor = free_free_solver.heating_factor(
            plasma.ion_number_density,
            electron_distribution.number_density.cgs.value,
        )
        free_free_cooling = free_free_solver.cooling_rate(
            electron_distribution.temperature.cgs.value,
            free_free_heating_factor,
        )
        total_cooling = exc_total + ion_total + fb_total + free_free_cooling
        free_free = np.divide(
            free_free_cooling,
            np.where(total_cooling == 0.0, 1.0, total_cooling),
        )
        collisional_excitation_rate.index = _macro_atom_collisional_index(
            collisional_excitation_rate.index
        )
        collisional_deexcitation_rate.index = _macro_atom_collisional_index(
            collisional_deexcitation_rate.index,
            reverse=True,
        )
        deexc_index = collisional_deexcitation_rate.index
        deexc_lower = pd.MultiIndex.from_arrays(
            [
                deexc_index.get_level_values("atomic_number"),
                deexc_index.get_level_values("ion_number"),
                deexc_index.get_level_values("level_number_lower"),
            ],
            names=["atomic_number", "ion_number", "level_number"],
        )
        deexc_upper = pd.MultiIndex.from_arrays(
            [
                deexc_index.get_level_values("atomic_number"),
                deexc_index.get_level_values("ion_number"),
                deexc_index.get_level_values("level_number_upper"),
            ],
            names=["atomic_number", "ion_number", "level_number"],
        )
        delta_E_yg = pd.Series(
            level_energy.loc[deexc_upper].values
            - level_energy.loc[deexc_lower].values,
            index=deexc_index,
        )
        excitation_references = pd.MultiIndex.from_arrays(
            [
                transition_index.get_level_values("atomic_number"),
                transition_index.get_level_values("ion_number"),
                transition_index.get_level_values("level_number_destination"),
            ],
            names=["atomic_number", "ion_number", "level_number"],
        )

        return cls(
            radiative_ionization_rate,
            radiative_recombination_rate,
            collisional_excitation_rate,
            collisional_deexcitation_rate,
            collisional_ionization_rate,
            collisional_recombination_rate,
            delta_E_yg,
            np.divide(
                exc_total, np.where(total_cooling == 0.0, 1.0, total_cooling)
            ),
            exc_array,
            excitation_references,
            np.divide(
                ion_total, np.where(total_cooling == 0.0, 1.0, total_cooling)
            ),
            ion_array,
            np.divide(
                fb_total, np.where(total_cooling == 0.0, 1.0, total_cooling)
            ),
            fb_array,
            free_free,
        )
