"""Structured continuum inputs consumed by the macro-atom solver."""

from dataclasses import dataclass

import numpy as np
import numpy.typing as npt
import pandas as pd
from astropy import units as u

from tardis import constants as const
from tardis.plasma.electron_energy_distribution import (
    ThermalElectronEnergyDistribution,
)
from tardis.plasma.equilibrium.inputs import ContinuumCoefficientState
from tardis.plasma.equilibrium.rates import (
    ThermalCollisionalRateSolver,
)
from tardis.plasma.equilibrium.rates.heating_cooling_rates import (
    BoundFreeThermalRates,
    FreeFreeThermalRates,
)


@dataclass(frozen=True)
class ContinuumMacroAtomState:
    """Continuum rates and cooling branches consumed by macro atoms."""

    radiative_ionization_rate: pd.DataFrame
    radiative_recombination_rate: pd.DataFrame
    collisional_excitation_rate: pd.DataFrame
    collisional_deexcitation_rate: pd.DataFrame
    collisional_ionization_rate: pd.DataFrame
    collisional_recombination_rate: pd.DataFrame
    delta_E_yg: pd.Series  # noqa: N815 - established macro-atom field name
    collisional_excitation_cooling_probability: npt.NDArray[np.float64]
    collisional_excitation_cooling_array: npt.NDArray[np.float64]
    collisional_excitation_references: pd.MultiIndex
    collisional_ionization_cooling_probability: npt.NDArray[np.float64]
    collisional_ionization_cooling_array: npt.NDArray[np.float64]
    radiative_recombination_cooling_probability: npt.NDArray[np.float64]
    radiative_recombination_cooling_array: npt.NDArray[np.float64]
    free_free_cooling_probability: npt.NDArray[np.float64]

    @classmethod
    def from_equilibrium(
        cls,
        atomic_data: object,
        lines: pd.DataFrame,
        photo_ion_cross_sections: pd.DataFrame,
        continuum_coefficients: ContinuumCoefficientState,
        level_to_continuum_saha_factor: pd.DataFrame,
        ion_number_density: pd.DataFrame,
        level_number_density: pd.DataFrame,
        electron_densities: pd.Series,
        t_electrons: npt.ArrayLike,
        stimulated_recombination_cooling_estimator: pd.DataFrame | None = None,
    ) -> ContinuumMacroAtomState:
        """Format accepted equilibrium rates for continuum macro atoms."""
        columns = level_number_density.columns
        photo_ion_index = photo_ion_cross_sections.index.unique()
        upper_ion_index = pd.MultiIndex.from_arrays(
            [
                photo_ion_index.get_level_values("atomic_number"),
                photo_ion_index.get_level_values("ion_number") + 1,
            ],
            names=["atomic_number", "ion_number"],
        )
        level_population = level_number_density.loc[photo_ion_index].to_numpy()
        stimulated_recombination_population = (
            level_to_continuum_saha_factor.loc[photo_ion_index].to_numpy()
            * ion_number_density.loc[upper_ion_index].to_numpy()
            * electron_densities.to_numpy()
        )
        population_correction = np.divide(
            stimulated_recombination_population,
            level_population,
            out=np.zeros((len(photo_ion_index), len(columns))),
            where=level_population != 0.0,
        )
        radiative_ionization_rate = pd.DataFrame(
            continuum_coefficients.photoionization.loc[
                photo_ion_index
            ].to_numpy()
            - population_correction
            * continuum_coefficients.stimulated_recombination.loc[
                photo_ion_index
            ].to_numpy(),
            index=photo_ion_index,
            columns=columns,
        )
        radiative_recombination_rate = (
            continuum_coefficients.spontaneous_recombination.loc[
                photo_ion_index
            ].set_axis(columns, axis="columns")
            * level_to_continuum_saha_factor.loc[photo_ion_index]
        )
        collisional_ionization_rate = (
            continuum_coefficients.collisional_ionization.loc[
                photo_ion_index
            ].set_axis(columns, axis="columns")
        )
        collisional_recombination_rate = (
            collisional_ionization_rate
            * level_to_continuum_saha_factor.loc[photo_ion_index]
        )

        electron_distribution = ThermalElectronEnergyDistribution(
            0.0 * u.erg,
            np.asarray(t_electrons) * u.K,
            electron_densities.to_numpy() / u.cm**3,
        )
        continuum_species = photo_ion_index.droplevel("level_number").unique()
        line_species = lines.index.droplevel(
            ["level_number_lower", "level_number_upper"]
        )
        continuum_lines = lines.loc[line_species.isin(continuum_species)]
        collisional_rates = ThermalCollisionalRateSolver(
            atomic_data.levels,
            continuum_lines,
            atomic_data.collision_data_temperatures,
            atomic_data.yg_data.loc[continuum_lines.index],
            collision_strengths_type="cmfgen",
        ).solve(electron_distribution.temperature)
        collisional_rates.columns = columns
        source_levels = collisional_rates.index.get_level_values(
            "level_number_source"
        )
        destination_levels = collisional_rates.index.get_level_values(
            "level_number_destination"
        )
        excitation_mask = source_levels < destination_levels
        collisional_excitation_rate = collisional_rates.loc[
            excitation_mask
        ].copy()
        collisional_deexcitation_rate = collisional_rates.loc[
            ~excitation_mask
        ].copy()

        transition_index = collisional_excitation_rate.index
        lower_level_index = pd.MultiIndex.from_arrays(
            [
                transition_index.get_level_values("atomic_number"),
                transition_index.get_level_values("ion_number"),
                transition_index.get_level_values("level_number_source"),
            ],
            names=["atomic_number", "ion_number", "level_number"],
        )
        upper_level_index = pd.MultiIndex.from_arrays(
            [
                transition_index.get_level_values("atomic_number"),
                transition_index.get_level_values("ion_number"),
                transition_index.get_level_values("level_number_destination"),
            ],
            names=["atomic_number", "ion_number", "level_number"],
        )
        excitation_energy = pd.Series(
            atomic_data.levels.energy.loc[upper_level_index].to_numpy()
            - atomic_data.levels.energy.loc[lower_level_index].to_numpy(),
            index=transition_index,
        )
        excitation_cooling = (
            collisional_excitation_rate
            * level_number_density.loc[lower_level_index].to_numpy()
            * electron_densities.to_numpy()
            * excitation_energy.to_numpy()[:, None]
        )
        threshold_energy = (
            photo_ion_cross_sections.groupby(level=[0, 1, 2]).first().nu
            * const.h.cgs.value
        )
        ionization_cooling = (
            collisional_ionization_rate
            * level_number_density.loc[photo_ion_index].to_numpy()
            * electron_densities.to_numpy()
            * threshold_energy.loc[photo_ion_index].to_numpy()[:, None]
        )
        recombination_cooling = BoundFreeThermalRates(
            photo_ion_cross_sections
        ).calculate_cooling_rate(
            ion_number_density,
            electron_distribution,
            level_to_continuum_saha_factor,
            stimulated_recombination_cooling_estimator,
        )

        excitation_total, excitation_array = _normalize_cooling_rates(
            excitation_cooling
        )
        ionization_total, ionization_array = _normalize_cooling_rates(
            ionization_cooling
        )
        recombination_total, recombination_array = _normalize_cooling_rates(
            recombination_cooling
        )
        free_free_solver = FreeFreeThermalRates()
        free_free_factor = free_free_solver.heating_factor(
            ion_number_density, electron_densities.to_numpy()
        )
        free_free_cooling = free_free_solver.cooling_rate(
            np.asarray(t_electrons), free_free_factor
        ).to_numpy()
        total_cooling = (
            excitation_total
            + ionization_total
            + recombination_total
            + free_free_cooling
        )
        safe_total_cooling = np.where(total_cooling == 0.0, 1.0, total_cooling)

        collisional_excitation_rate.index = _macro_atom_collisional_index(
            collisional_excitation_rate.index
        )
        collisional_deexcitation_rate.index = _macro_atom_collisional_index(
            collisional_deexcitation_rate.index, reverse=True
        )
        deexcitation_index = collisional_deexcitation_rate.index
        deexcitation_lower_index = pd.MultiIndex.from_arrays(
            [
                deexcitation_index.get_level_values("atomic_number"),
                deexcitation_index.get_level_values("ion_number"),
                deexcitation_index.get_level_values("level_number_lower"),
            ],
            names=["atomic_number", "ion_number", "level_number"],
        )
        deexcitation_upper_index = pd.MultiIndex.from_arrays(
            [
                deexcitation_index.get_level_values("atomic_number"),
                deexcitation_index.get_level_values("ion_number"),
                deexcitation_index.get_level_values("level_number_upper"),
            ],
            names=["atomic_number", "ion_number", "level_number"],
        )
        delta_E_yg = pd.Series(
            atomic_data.levels.energy.loc[deexcitation_upper_index].to_numpy()
            - atomic_data.levels.energy.loc[
                deexcitation_lower_index
            ].to_numpy(),
            index=deexcitation_index,
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
            excitation_total / safe_total_cooling,
            excitation_array,
            excitation_references,
            ionization_total / safe_total_cooling,
            ionization_array,
            recombination_total / safe_total_cooling,
            recombination_array,
            free_free_cooling / safe_total_cooling,
        )


def _normalize_cooling_rates(
    rates: pd.DataFrame,
) -> tuple[npt.NDArray[np.float64], npt.NDArray[np.float64]]:
    """Return shell totals and shell-by-transition branch probabilities."""
    totals = rates.sum(axis=0).to_numpy(dtype=np.float64)
    safe_totals = np.where(totals == 0.0, 1.0, totals)
    probabilities = rates.to_numpy(dtype=np.float64).T / safe_totals[:, None]
    return totals, np.ascontiguousarray(probabilities)


def _macro_atom_collisional_index(
    index: pd.MultiIndex, *, reverse: bool = False
) -> pd.MultiIndex:
    """Use the transition index names expected by macro-atom kernels."""
    if reverse:
        index = index.swaplevel(
            "level_number_source", "level_number_destination"
        )
    redundant_levels = [
        level_name
        for level_name in (
            "ion_number_source",
            "ion_number_destination",
        )
        if level_name in index.names
    ]
    if redundant_levels:
        index = index.droplevel(redundant_levels)
    return index.rename(
        {
            "level_number_source": (
                "level_number_upper" if reverse else "level_number_lower"
            ),
            "level_number_destination": (
                "level_number_lower" if reverse else "level_number_upper"
            ),
        }
    )
