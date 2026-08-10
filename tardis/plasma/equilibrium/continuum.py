from collections.abc import Collection
from dataclasses import dataclass
from itertools import pairwise

import numpy as np
import numpy.typing as npt
import pandas as pd
from astropy import units as u

from tardis import constants as const
from tardis.plasma.electron_energy_distribution import (
    ThermalElectronEnergyDistribution,
)
from tardis.plasma.equilibrium.ion_populations import (
    AnalyticIonPopulationSolver,
)
from tardis.plasma.equilibrium.population_state import (
    PopulationState,
    SingleElementPopulationState,
)
from tardis.plasma.equilibrium.rate_matrix import (
    IonRateMatrix,
    LevelRateMatrix,
)
from tardis.plasma.equilibrium.rates import (
    AnalyticCorrectedPhotoionizationCoeffSolver,
    CollisionalIonizationSeaton,
    EstimatedPhotoionizationCoeffSolver,
    RadiativeRatesSolver,
    SpontaneousRecombinationCoeffSolver,
    ThermalCollisionalRateSolver,
)
from tardis.plasma.equilibrium.rates.heating_cooling_rates import (
    BoundFreeThermalRates,
    FreeFreeThermalRates,
)
from tardis.plasma.equilibrium.rates.photoionization_strengths import (
    apply_lte_level_to_ion_factor,
)
from tardis.plasma.equilibrium.rates.util import (
    reindex_ionization_rate_dataframe,
    reindex_level_resolved_recombination_rate_dataframe,
)
from tardis.plasma.properties.ion_population import (
    SahaFactor,
    ThermalPhiSahaLTE,
)
from tardis.plasma.properties.level_population import LevelNumberDensity


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


def select_species_rows(
    frame: pd.DataFrame,
    species: Collection[tuple[int, int]],
) -> pd.DataFrame:
    """Select rows whose atomic number and ion number are configured."""
    selected_species = set(species)
    if frame.empty or not selected_species:
        return frame.iloc[0:0]
    atomic_numbers = frame.index.get_level_values("atomic_number")
    ion_numbers = frame.index.get_level_values("ion_number")
    mask = np.fromiter(
        (
            (int(atomic_number), int(ion_number)) in selected_species
            for atomic_number, ion_number in zip(
                atomic_numbers, ion_numbers, strict=True
            )
        ),
        dtype=bool,
        count=len(frame),
    )
    return frame.loc[mask]


@dataclass(frozen=True)
class PreparedContinuumRates:
    """Density-independent continuum coefficients for charge solves."""

    radiative_rates: pd.DataFrame
    collisional_rates: pd.DataFrame
    photoionization_rates: pd.DataFrame
    recombination_coefficients: pd.DataFrame
    collisional_ionization_coefficients: pd.DataFrame
    collisional_recombination_coefficients: pd.DataFrame


def prepare_continuum_rates(
    *,
    atomic_data: object,
    lines: pd.DataFrame,
    continuum_species: Collection[tuple[int, int]],
    radiation_field: object,
    electron_temperatures: pd.Series,
    thermal_g_electron: object,
    beta_electron: object,
    thermal_lte_level_boltzmann_factor: pd.DataFrame,
    thermal_lte_partition_function: pd.DataFrame,
    ionization_data: pd.DataFrame,
    photoionization_rate_estimator: object | None,
    stimulated_recombination_rate_estimator: object | None,
    beta_sobolev: pd.DataFrame,
) -> PreparedContinuumRates:
    """Build continuum rates that do not depend on the trial density."""
    columns = electron_temperatures.index
    selected_species = frozenset(continuum_species)
    line_atomic_numbers = lines.index.get_level_values("atomic_number")
    line_ion_numbers = lines.index.get_level_values("ion_number")
    line_mask = np.fromiter(
        (
            (int(atomic_number), int(ion_number)) in selected_species
            for atomic_number, ion_number in zip(
                line_atomic_numbers, line_ion_numbers, strict=True
            )
        ),
        dtype=bool,
        count=len(lines),
    )
    selected_lines = lines.loc[line_mask]
    prepared_beta_sobolev = beta_sobolev.reindex(
        selected_lines.index, columns=columns
    ).copy(deep=True)
    radiative_rates = RadiativeRatesSolver(selected_lines).solve(
        radiation_field,
        beta_sobolevs=prepared_beta_sobolev,
    )
    radiative_rates = radiative_rates.reindex(columns=columns)
    collisional_rates = ThermalCollisionalRateSolver(
        atomic_data.levels,
        selected_lines,
        atomic_data.collision_data_temperatures,
        atomic_data.yg_data.loc[selected_lines.index],
        "cmfgen",
    ).solve(electron_temperatures.to_numpy() * u.K)
    collisional_rates.columns = columns

    estimators = {
        "photoionization_rate_estimator": photoionization_rate_estimator,
        "stimulated_recombination_rate_estimator": (
            stimulated_recombination_rate_estimator
        ),
    }
    photoionization_coefficients, stimulated_coefficients = (
        EstimatedPhotoionizationCoeffSolver(
            atomic_data.level2continuum_edge_idx
        ).solve(estimators)
    )
    photoionization_coefficients = select_species_rows(
        photoionization_coefficients, selected_species
    ).reindex(columns=columns)
    stimulated_coefficients = select_species_rows(
        stimulated_coefficients, selected_species
    ).reindex(columns=columns)
    spontaneous_coefficients = SpontaneousRecombinationCoeffSolver(
        atomic_data.photoionization_data
    ).solve(electron_temperatures.to_numpy() * u.K)
    spontaneous_coefficients.columns = columns
    spontaneous_coefficients = select_species_rows(
        spontaneous_coefficients, selected_species
    )
    thermal_phi_lte = ThermalPhiSahaLTE.calculate(
        thermal_g_electron,
        beta_electron,
        thermal_lte_partition_function,
        ionization_data,
    )
    level_to_ion_factor = SahaFactor(None).calculate(
        thermal_phi_lte.loc[:, columns],
        thermal_lte_level_boltzmann_factor.loc[:, columns],
        thermal_lte_partition_function.loc[:, columns],
    )
    recombination_coefficients = apply_lte_level_to_ion_factor(
        spontaneous_coefficients.add(stimulated_coefficients, fill_value=0.0),
        level_to_ion_factor.reindex(spontaneous_coefficients.index),
    )
    collisional_strengths = CollisionalIonizationSeaton(
        atomic_data.photoionization_data
    ).solve(electron_temperatures.to_numpy() * u.K)
    collisional_strengths.columns = columns
    collisional_strengths = select_species_rows(
        collisional_strengths, selected_species
    )
    collisional_recombination_coefficients = collisional_strengths.multiply(
        level_to_ion_factor.reindex(collisional_strengths.index)
    )
    return PreparedContinuumRates(
        radiative_rates=radiative_rates,
        collisional_rates=collisional_rates,
        photoionization_rates=reindex_ionization_rate_dataframe(
            photoionization_coefficients,
            recombination=False,
        ),
        recombination_coefficients=(
            reindex_level_resolved_recombination_rate_dataframe(
                recombination_coefficients
            )
        ),
        collisional_ionization_coefficients=(
            reindex_ionization_rate_dataframe(
                collisional_strengths,
                recombination=False,
            )
        ),
        collisional_recombination_coefficients=(
            reindex_level_resolved_recombination_rate_dataframe(
                collisional_recombination_coefficients
            )
        ),
    )


class ContinuumPopulationSolver(AnalyticIonPopulationSolver):
    """Solve continuum populations at trial electron densities."""

    _TRANSITION_INDEX_NAMES = [
        "atomic_number",
        "ion_number",
        "ion_number_source",
        "ion_number_destination",
        "level_number_source",
        "level_number_destination",
    ]

    def __init__(
        self,
        *,
        atomic_data: object,
        lines: pd.DataFrame,
        continuum_species: Collection[tuple[int, int]],
        radiation_field: object,
        electron_temperatures: npt.ArrayLike,
        elemental_number_density: pd.DataFrame,
        general_level_boltzmann_factor: pd.DataFrame,
        general_partition_function: pd.DataFrame,
        thermal_g_electron: npt.ArrayLike,
        beta_electron: npt.ArrayLike,
        thermal_lte_level_boltzmann_factor: pd.DataFrame,
        thermal_lte_partition_function: pd.DataFrame,
        ionization_data: pd.DataFrame,
        nebular_phi: pd.DataFrame,
        photoionization_rate_estimator: object | None,
        stimulated_recombination_rate_estimator: object | None,
        max_solver_iterations: int = 100,
    ) -> None:
        super().__init__(
            IonRateMatrix(),
            None,
            None,
            elemental_number_density,
            max_solver_iterations=max_solver_iterations,
        )
        columns = elemental_number_density.columns
        electron_temperature_values = np.asarray(electron_temperatures)
        if electron_temperature_values.shape != (len(columns),):
            raise ValueError(
                "electron_temperatures must contain one value per shell"
            )
        self._continuum_atomic_data = atomic_data
        self._continuum_lines = lines
        self._continuum_species = frozenset(continuum_species)
        self._continuum_radiation_field = radiation_field
        self._continuum_electron_temperatures = pd.Series(
            electron_temperature_values, index=columns
        )
        self._continuum_general_level_boltzmann_factor = (
            general_level_boltzmann_factor.copy(deep=True)
        )
        self._continuum_general_partition_function = (
            general_partition_function.copy(deep=True)
        )
        self._continuum_thermal_g_electron = np.asarray(thermal_g_electron)
        self._continuum_beta_electron = np.asarray(beta_electron)
        self._continuum_thermal_lte_level_boltzmann_factor = (
            thermal_lte_level_boltzmann_factor.copy(deep=True)
        )
        self._continuum_thermal_lte_partition_function = (
            thermal_lte_partition_function.copy(deep=True)
        )
        self._continuum_ionization_data = ionization_data.copy(deep=True)
        self._continuum_nebular_phi = nebular_phi.copy(deep=True)
        self._continuum_photoionization_rate_estimator = (
            photoionization_rate_estimator
        )
        self._continuum_stimulated_recombination_rate_estimator = (
            stimulated_recombination_rate_estimator
        )

    def prepare_charge_conservation_solve(
        self, beta_sobolev: pd.DataFrame | None
    ) -> None:
        """Prepare continuum coefficients for one charge-conservation solve."""
        if not self._continuum_species:
            return
        if beta_sobolev is None:
            raise ValueError(
                "Configured continuum rates require Sobolev escape probabilities"
            )
        self._continuum_charge_conservation_rates = prepare_continuum_rates(
            atomic_data=self._continuum_atomic_data,
            lines=self._continuum_lines,
            continuum_species=self._continuum_species,
            radiation_field=self._continuum_radiation_field,
            electron_temperatures=self._continuum_electron_temperatures,
            thermal_g_electron=self._continuum_thermal_g_electron,
            beta_electron=self._continuum_beta_electron,
            thermal_lte_level_boltzmann_factor=(
                self._continuum_thermal_lte_level_boltzmann_factor
            ),
            thermal_lte_partition_function=(
                self._continuum_thermal_lte_partition_function
            ),
            ionization_data=self._continuum_ionization_data,
            photoionization_rate_estimator=(
                self._continuum_photoionization_rate_estimator
            ),
            stimulated_recombination_rate_estimator=(
                self._continuum_stimulated_recombination_rate_estimator
            ),
            beta_sobolev=beta_sobolev,
        )

    @classmethod
    def _empty_transition_rates(cls, columns: pd.Index) -> pd.DataFrame:
        """Return an empty rate frame with the charge-rate transition index."""
        return pd.DataFrame(
            index=pd.MultiIndex.from_tuples(
                [], names=cls._TRANSITION_INDEX_NAMES
            ),
            columns=columns,
            dtype=float,
        )

    @staticmethod
    def _calculate_level_factors_from_solution(
        general_level_boltzmann_factor: pd.DataFrame,
        solution: SingleElementPopulationState,
        species: tuple[int, int],
    ) -> pd.DataFrame:
        """Convert solved within-ion level fractions to level factors."""
        level_fraction = solution.normalized_level_populations.loc[species]
        ion_fraction = solution.normalized_ion_populations.loc[species]
        level_fraction = level_fraction.divide(ion_fraction, axis="columns")
        reference = general_level_boltzmann_factor.loc[species]
        ground_fraction = level_fraction.iloc[0].replace(0.0, np.nan)
        factors = level_fraction.multiply(
            reference.iloc[0] / ground_fraction, axis="columns"
        )
        return factors.fillna(reference)

    def _build_nebular_rates(
        self,
        electron_number_density: pd.Series,
    ) -> pd.DataFrame:
        """Build bidirectional rates for every uncovered ionization edge."""
        transitions: list[tuple[int, int, int, int, int, int]] = []
        rates: list[npt.NDArray[np.float64]] = []
        columns = electron_number_density.index
        nebular_phi = self._continuum_nebular_phi.loc[:, columns].copy(
            deep=True
        )
        zero_phi = nebular_phi == 0.0
        if zero_phi.any().any():
            positive_phi = (
                self._continuum_nebular_phi.where(
                    self._continuum_nebular_phi > 0.0
                )
                .min()
                .min()
            )
            if not np.isfinite(positive_phi):
                raise ValueError(
                    "Nebular phi must contain a positive value when zero "
                    "ionization factors are present"
                )
            nebular_phi = nebular_phi.mask(zero_phi, 1e-10 * positive_phi)
        for atomic_number in self.elemental_number_density.index:
            configured_ions = {
                ion_number
                for species_atomic_number, ion_number in self._continuum_species
                if species_atomic_number == atomic_number
            }
            for ion_number in range(int(atomic_number)):
                if ion_number in configured_ions:
                    continue
                phi_index = (atomic_number, ion_number + 1)
                try:
                    ionization_rate = nebular_phi.loc[phi_index, columns]
                except KeyError as exc:
                    raise ValueError(
                        "Nebular phi does not cover ionization edge "
                        f"{phi_index}"
                    ) from exc
                transitions.extend(
                    [
                        (
                            atomic_number,
                            ion_number,
                            ion_number,
                            ion_number + 1,
                            0,
                            0,
                        ),
                        (
                            atomic_number,
                            ion_number,
                            ion_number + 1,
                            ion_number,
                            0,
                            0,
                        ),
                    ]
                )
                rates.extend(
                    [
                        np.asarray(ionization_rate, dtype=float),
                        electron_number_density.to_numpy(dtype=float),
                    ]
                )
        return pd.DataFrame(
            rates,
            index=pd.MultiIndex.from_tuples(
                transitions, names=self._TRANSITION_INDEX_NAMES
            ),
            columns=columns,
            dtype=float,
        )

    def _calculate_continuum_rate_frames(
        self,
        electron_number_density: pd.Series,
    ) -> tuple[pd.DataFrame | None, tuple[pd.DataFrame, ...]]:
        """Apply trial electron-density scaling to continuum coefficients."""
        columns = electron_number_density.index
        if not self._continuum_species:
            empty_rates = self._empty_transition_rates(columns)
            return None, (empty_rates,) * 4
        charge_conservation_rates = self._continuum_charge_conservation_rates
        radiative_rates = charge_conservation_rates.radiative_rates.loc[
            :, columns
        ]
        collisional_rates = charge_conservation_rates.collisional_rates.loc[
            :, columns
        ]
        raw_level_rate_matrices = LevelRateMatrix(
            self._continuum_atomic_data.levels
        ).build_raw_rate_matrices(
            radiative_rates,
            collisional_rates,
            electron_number_density.to_numpy(),
        )
        photoionization_rates = (
            charge_conservation_rates.photoionization_rates.loc[:, columns]
        )
        recombination_coefficients = (
            charge_conservation_rates.recombination_coefficients.loc[:, columns]
        )
        recombination_rates = recombination_coefficients.multiply(
            electron_number_density, axis="columns"
        )
        collisional_ionization_rates = (
            charge_conservation_rates.collisional_ionization_coefficients.loc[
                :, columns
            ].multiply(electron_number_density, axis="columns")
        )
        collisional_recombination_rates = charge_conservation_rates.collisional_recombination_coefficients.loc[
            :, columns
        ].multiply(electron_number_density**2, axis="columns")
        return raw_level_rate_matrices, (
            photoionization_rates,
            recombination_rates,
            collisional_ionization_rates,
            collisional_recombination_rates,
        )

    def solve_charge_state_at_electron_density(
        self,
        electron_number_density: pd.Series,
    ) -> PopulationState:
        """Solve every element at a trial electron density."""
        columns = electron_number_density.index
        column_positions = self.elemental_number_density.columns.get_indexer(
            columns
        )
        if (column_positions < 0).any():
            raise ValueError("Unknown shell in electron_number_density")
        columns = self.elemental_number_density.columns.take(column_positions)
        electron_number_density = pd.Series(
            electron_number_density.to_numpy(dtype=float), index=columns
        )
        selected_elemental_density = self.elemental_number_density.loc[
            :, columns
        ]
        if (selected_elemental_density <= 0.0).any().any():
            raise ValueError(
                "Continuum charge-state solves require a positive abundance "
                "for every selected element and shell"
            )
        (
            raw_level_rate_matrices,
            configured_rate_frames,
        ) = self._calculate_continuum_rate_frames(electron_number_density)
        nebular_rates = self._build_nebular_rates(electron_number_density)
        elemental_populations: dict[int, SingleElementPopulationState] = {}
        for atomic_number in self.elemental_number_density.index:
            atomic_number_mask = (
                nebular_rates.index.get_level_values("atomic_number")
                == atomic_number
            )
            elemental_rate_frames = tuple(
                rate_frame.loc[
                    rate_frame.index.get_level_values("atomic_number")
                    == atomic_number
                ]
                for rate_frame in configured_rate_frames
            )
            elemental_nebular_rates = nebular_rates.loc[atomic_number_mask]
            selected_level_matrices = None
            if raw_level_rate_matrices is not None:
                level_mask = (
                    raw_level_rate_matrices.index.get_level_values(
                        "atomic_number"
                    )
                    == atomic_number
                )
                if level_mask.any():
                    selected_level_matrices = raw_level_rate_matrices.loc[
                        level_mask
                    ].copy(deep=True)
            matrix_set = self.rate_matrix_solver.solve_ion_and_level(
                int(atomic_number),
                raw_level_rate_matrices=selected_level_matrices,
                photoion_rates_df=elemental_rate_frames[0],
                recomb_rates_df=elemental_rate_frames[1],
                collisional_ionization_rates_df=elemental_rate_frames[2],
                collision_recombination_rates_df=elemental_rate_frames[3],
                nebular_rates_df=elemental_nebular_rates,
            )
            elemental_populations[int(atomic_number)] = (
                self._build_elemental_population_solution(
                    matrix_set,
                    self.elemental_number_density.loc[atomic_number, columns],
                )
            )

        ion_number_density = pd.concat(
            [
                solution.ion_populations
                for solution in elemental_populations.values()
            ]
        ).sort_index()
        ion_number_density.columns = columns
        level_boltzmann_factor = (
            self._continuum_general_level_boltzmann_factor.loc[:, columns].copy(
                deep=True
            )
        )
        level_number_density = LevelNumberDensity(None).calculate(
            level_boltzmann_factor,
            ion_number_density,
            level_boltzmann_factor.index,
            self._continuum_general_partition_function.loc[:, columns],
        )
        level_number_density.columns = columns
        for atomic_number, solution in elemental_populations.items():
            for species in self._continuum_species:
                if species[0] == atomic_number:
                    level_boltzmann_factor.loc[species] = (
                        self._calculate_level_factors_from_solution(
                            level_boltzmann_factor,
                            solution,
                            species,
                        )
                    )
            level_number_density.loc[solution.level_populations.index] = (
                solution.level_populations
            )
        return PopulationState(
            electron_densities=electron_number_density,
            elemental_populations=elemental_populations,
            ion_number_density=ion_number_density,
            level_number_density=level_number_density,
            level_boltzmann_factor=level_boltzmann_factor,
        )


@dataclass
class ContinuumRateState:
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


def build_continuum_rate_state(
    plasma: object,
    estimators: dict[str, object] | None = None,
) -> ContinuumRateState:
    """Build continuum inputs from equilibrium-plasma outputs."""
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
        radiative_ionization_rate = corrected_solver.solve_corrected(
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
    radiative_recombination_rate = apply_lte_level_to_ion_factor(
        spontaneous, level_factor
    )

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

    collisional_ionization_rate = CollisionalIonizationSeaton(photo_data).solve(
        electron_distribution.temperature
    )
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

    return ContinuumRateState(
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
        np.divide(fb_total, np.where(total_cooling == 0.0, 1.0, total_cooling)),
        fb_array,
        free_free,
    )
