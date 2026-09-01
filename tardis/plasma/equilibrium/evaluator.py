from collections.abc import Mapping
from copy import deepcopy
from dataclasses import dataclass

import astropy.units as u
import numpy as np
import numpy.typing as npt
import pandas as pd
from scipy.optimize import root

from tardis.opacities.sobolevs_from_levels import (
    calculate_sobolev_opacities_from_level_densities,
)
from tardis.plasma.electron_energy_distribution import (
    ThermalElectronEnergyDistribution,
)
from tardis.plasma.equilibrium.inputs import (
    ContinuumRateCoefficients,
    LevelEquationRates,
    ShellNumberDensity,
    SobolevInputs,
)
from tardis.plasma.equilibrium.rate_matrix import RateMatrix
from tardis.plasma.equilibrium.rates.collisional_ionization_rates import (
    CollisionalIonizationRateSolver,
)
from tardis.plasma.equilibrium.rates.photoionization_rates import (
    EstimatedPhotoionizationRateSolver,
)
from tardis.plasma.properties.general import BetaElectron, ThermalGElectron
from tardis.plasma.properties.ion_population import (
    SahaFactor,
    ThermalPhiSahaLTE,
)
from tardis.plasma.properties.partition_function import (
    ThermalLevelBoltzmannFactorLTE,
    ThermalLTEPartitionFunction,
)
from tardis.transport.montecarlo.estimators.estimators_continuum import (
    EstimatorsContinuum,
)


@dataclass(frozen=True)
class PlasmaEquilibriumEvaluation:
    """Fixed-density and terminal equilibrium outputs."""

    trial_electron_density: pd.Series
    normalized_population: pd.DataFrame
    diagnostic_ion_ratio: pd.Series
    trial_beta_sobolev: pd.DataFrame
    trial_level_residual: pd.DataFrame
    charge_solved_electron_density: pd.Series | None
    absolute_level_population: pd.DataFrame
    ion_population: pd.DataFrame | None
    tau_sobolev: pd.DataFrame
    beta_sobolev: pd.DataFrame
    level_residual: pd.DataFrame
    charge_residual: pd.Series | None
    electron_residual: pd.Series | None
    total_heating: pd.Series | None
    fractional_heating: pd.Series | None


def calculate_nlte_level_population_residual(
    level_fractions: npt.NDArray[np.float64],
    level_rates: LevelEquationRates,
    rate_matrix_solver: RateMatrix,
    j_blues: pd.DataFrame,
    thermal_electron_energy_distribution: ThermalElectronEnergyDistribution,
    species: tuple[int, int],
    shell_number_density: ShellNumberDensity,
    sobolev: SobolevInputs,
    level_density: npt.NDArray[np.float64] | None = None,
) -> tuple[npt.NDArray[np.float64], npt.NDArray[np.float64], float]:
    """Calculate the reduced fixed-density NLTE level residual.

    The auxiliary continuum ratio is eliminated algebraically as

    ``ionized_to_neutral_ratio = (total_ionization_rates @ b)
    / total_recombination_rate``.

    Parameters
    ----------
    level_fractions : numpy.ndarray
        Candidate normalized bound-level fractions.
    level_rates : LevelEquationRates
        Density-specific ionization, recombination, and ionization-loss rates.
    rate_matrix_solver : RateMatrix
        Shared bound-bound matrix owner.
    j_blues : pandas.DataFrame
        Fixed mean intensities for this residual shell.
    thermal_electron_energy_distribution : ThermalElectronEnergyDistribution
        Candidate shell electron distribution.
    species : tuple[int, int]
        Atomic and ion number of the reduced NLTE species.
    number_density_per_shell : NumberDensityPerShell
        Absolute level-density state and selected-level positions.
    sobolev : SobolevInputs
        Line geometry used to calculate candidate beta.
    level_density : numpy.ndarray, optional
        Authoritative absolute level-density state used for final-state
        closure. When omitted, the temporary state implied by the
        ionized-to-neutral ratio is reconstructed from
        reconstructed from ``hydrogen_number_density``.

    Returns
    -------
    tuple[numpy.ndarray, numpy.ndarray, float]
        Component level residual, candidate Sobolev escape probabilities, and
        the reduced ionized-to-neutral ratio.
    """
    level_fractions = np.asarray(level_fractions, dtype=np.float64)
    total_ionization_rates = level_rates.ionization

    total_recombination_rates = level_rates.recombination.copy()
    total_recombination_rates[0] = 0.0
    total_recombination_rate = total_recombination_rates.sum()

    ionized_to_neutral_ratio = float(
        np.dot(total_ionization_rates, level_fractions)
        / total_recombination_rate
    )

    if level_density is None:
        level_density = shell_number_density.level_number_density.copy()
        level_density[shell_number_density.species_level_positions] = (
            shell_number_density.hydrogen_number_density
            / (1.0 + ionized_to_neutral_ratio)
            * level_fractions
        )
    else:
        level_density = np.asarray(level_density, dtype=np.float64).copy()

    _, beta_sobolev = calculate_sobolev_opacities_from_level_densities(
        level_density,
        sobolev.line_indices,
        sobolev.lines_lower_level_index,
        sobolev.lines_upper_level_index,
        sobolev.g_lower,
        sobolev.g_upper,
        sobolev.metastable_upper,
        sobolev.nlte_lines_mask,
        sobolev.tau_coefficient,
    )

    sobolev_escape_probability_frame = pd.DataFrame(
        beta_sobolev[:, None],
        index=sobolev.line_index,
        columns=j_blues.columns,
    )
    level_rate_matrix = (
        rate_matrix_solver.assemble_matrices(
            j_blues,
            thermal_electron_energy_distribution,
            sobolev_escape_probability_frame,
        )
        .loc[species]
        .iloc[0]
        .copy()
    )
    level_rate_matrix += level_rates.ionization_loss_matrix
    level_rate_matrix[0, :] = 1.0

    recombination_rates_column = total_recombination_rates.copy()
    normalization_vector = np.zeros(level_fractions.size, dtype=np.float64)
    normalization_vector[0] = 1.0
    level_residual = (
        level_rate_matrix @ level_fractions
        + recombination_rates_column * ionized_to_neutral_ratio
        - normalization_vector
    )
    return level_residual, beta_sobolev, ionized_to_neutral_ratio


@dataclass(frozen=True)
class _LevelSolveState:
    """Intermediate level results shared by the optional final stages."""

    continuum_rate_coefficients: tuple[ContinuumRateCoefficients, ...]
    level_to_continuum_saha_factor: pd.DataFrame
    collisional_ionization_rate_coefficient: pd.DataFrame
    lte_ionization_factor: pd.DataFrame
    fractions: tuple[npt.NDArray[np.float64], ...]
    ionized_to_neutral_ratios: tuple[float, ...]
    full_population: pd.DataFrame
    normalized_population: pd.DataFrame
    trial_escape_probability: pd.DataFrame
    trial_level_residual: pd.DataFrame
    hydrogen_index: pd.MultiIndex


class PlasmaEquilibriumEvaluator:
    """Evaluate the reduced fixed-density equilibrium composition.

    The four per-shell input tuples are kept separate so each level-rate
    dependency remains visible at construction and evaluation call sites.
    """

    def __init__(
        self,
        photoionization_cross_sections: pd.DataFrame,
        level2continuum_edge_idx: pd.Series,
        estimators_continuum: EstimatorsContinuum,
        time_simulation: u.Quantity,
        volume: u.Quantity,
        levels: pd.DataFrame,
        ionization_data: pd.Series,
        rate_matrix_solver: RateMatrix,
        j_blues: pd.DataFrame,
        shell_number_densities: tuple[ShellNumberDensity, ...],
        sobolev_inputs: tuple[SobolevInputs, ...],
        level_population_index: pd.MultiIndex,
        hydrogen_species: tuple[int, int],
        elemental_number_density: pd.DataFrame,
        maximum_electron_density: npt.ArrayLike,
        ion_population_solver: object | None = None,
        ion_population_arguments: Mapping[str, object] | None = None,
        thermal_balance_solver: object | None = None,
        thermal_balance_arguments: Mapping[str, object] | None = None,
        reference_electron_temperature: u.Quantity | None = None,
    ) -> None:
        """Construct an evaluator from fixed scientific inputs.

        Parameters
        ----------
        photoionization_cross_sections : pandas.DataFrame
            Photoionization cross sections for the selected continuum species.
        level2continuum_edge_idx : pandas.Series
            Continuum-estimator positions indexed by level.
        estimators_continuum : EstimatorsContinuum
            Fixed post-Monte-Carlo continuum estimators.
        time_simulation : astropy.units.Quantity
            Monte Carlo simulation time used to normalize estimators.
        volume : astropy.units.Quantity
            Shell volumes used to normalize estimators.
        levels : pandas.DataFrame
            Atomic level energies and statistical weights.
        ionization_data : pandas.Series
            Ionization energies indexed by atomic and ion number.
        rate_matrix_solver : RateMatrix
            Shared bound-bound matrix owner.
        j_blues : pandas.DataFrame
            Fixed post-Monte-Carlo line mean intensities.
        shell_number_densities : tuple[ShellNumberDensity, ...]
            Per-shell absolute number densities, level number densities and species level positions.
        sobolev_inputs : tuple[SobolevInputs, ...]
            Per-shell Sobolev line inputs.
        level_population_index : pandas.MultiIndex
            Complete level-population index used for returned absolute levels.
        hydrogen_species : tuple[int, int]
            Atomic and ion number of the reduced NLTE species.
        elemental_number_density : pandas.DataFrame
            Elemental number density indexed by atomic number and shell.
        maximum_electron_density : array-like
            Maximum electron density used by the charge residual.
        ion_population_solver : IonPopulationSolver, optional
            Existing authoritative charge solver.
        ion_population_arguments : mapping, optional
            Fixed keyword arguments for ``ion_population_solver.solve``.
        thermal_balance_solver : ThermalBalanceSolver, optional
            Existing thermal-rate owner.
        thermal_balance_arguments : mapping, optional
            Fixed keyword arguments for ``thermal_balance_solver.solve``.
        reference_electron_temperature : astropy.units.Quantity, optional
            Temperature at which fixed thermal rate coefficients were built.
        """
        self.photoionization_cross_sections = photoionization_cross_sections
        self.level2continuum_edge_idx = level2continuum_edge_idx
        self.estimators_continuum = estimators_continuum
        self.time_simulation = time_simulation
        self.volume = volume
        self.estimated_photoionization_rate_solver = (
            EstimatedPhotoionizationRateSolver(
                photoionization_cross_sections,
                level2continuum_edge_idx,
                estimators_continuum,
                time_simulation,
                volume,
            )
        )
        self.collisional_ionization_rate_solver = (
            CollisionalIonizationRateSolver(photoionization_cross_sections)
        )
        self.levels = levels
        self.ionization_data = ionization_data
        self.rate_matrix_solver = rate_matrix_solver
        self.j_blues = j_blues.copy(deep=True)
        self.shell_number_densities = shell_number_densities
        self.sobolev_inputs = sobolev_inputs
        self.level_population_index = level_population_index
        self.hydrogen_species = hydrogen_species
        self.elemental_number_density = elemental_number_density.copy(deep=True)
        self.maximum_electron_density = np.asarray(
            maximum_electron_density, dtype=np.float64
        )
        self.ion_population_solver = ion_population_solver
        self.ion_population_arguments = deepcopy(ion_population_arguments or {})
        self.thermal_balance_solver = thermal_balance_solver
        self.thermal_balance_arguments = deepcopy(
            thermal_balance_arguments or {}
        )
        self.reference_electron_temperature = reference_electron_temperature

        shell_count = len(self.elemental_number_density.columns)
        if not all(
            len(inputs) == shell_count
            for inputs in (shell_number_densities, sobolev_inputs)
        ):
            raise ValueError(
                "Each level input tuple must contain one entry per shell."
            )
        if len(self.j_blues.columns) != shell_count:
            raise ValueError("j_blues must contain one column per shell.")

    def _calculate_continuum_rate_coefficients(
        self, electron_temperature: npt.NDArray[np.float64]
    ) -> tuple[
        tuple[ContinuumRateCoefficients, ...],
        pd.DataFrame,
        pd.DataFrame,
        pd.DataFrame,
    ]:
        """Calculate candidate-temperature continuum coefficients."""
        (
            photoionization,
            stimulated_recombination,
            spontaneous_recombination,
        ) = self.estimated_photoionization_rate_solver.solve_coefficients(
            electron_temperature * u.K
        )
        collisional_ionization = (
            self.collisional_ionization_rate_solver.solve_coefficients(
                electron_temperature * u.K
            )
        )
        beta_electron = BetaElectron(None).calculate(electron_temperature)
        level_boltzmann_factor = ThermalLevelBoltzmannFactorLTE(None).calculate(
            self.levels.energy,
            self.levels.g,
            beta_electron,
            self.levels.index,
        )
        partition_function = ThermalLTEPartitionFunction(None).calculate(
            level_boltzmann_factor
        )
        thermal_saha_factor = ThermalPhiSahaLTE(None).calculate(
            ThermalGElectron(None).calculate(beta_electron),
            beta_electron,
            partition_function,
            self.ionization_data,
        )
        level_to_continuum_saha_factor = SahaFactor(None).calculate(
            thermal_saha_factor,
            level_boltzmann_factor,
            partition_function,
        )
        species_saha_factor = level_to_continuum_saha_factor.loc[
            self.hydrogen_species
        ]
        shell_index = self.elemental_number_density.columns
        species_index = species_saha_factor.index
        if len(species_saha_factor.columns) != len(shell_index):
            raise ValueError(
                "Level-to-continuum Saha factors must contain one column per shell."
            )
        species_saha_factor.columns = shell_index
        level_to_continuum_saha_factor.columns = shell_index
        collisional_ionization_rate_coefficient = collisional_ionization.copy()
        collisional_ionization_rate_coefficient.columns = shell_index
        thermal_saha_factor.columns = shell_index
        rate_frames = []
        for frame in (
            photoionization,
            stimulated_recombination,
            spontaneous_recombination,
            collisional_ionization,
        ):
            species_frame = frame.loc[self.hydrogen_species].reindex(
                index=species_index
            )
            if len(species_frame.columns) != len(shell_index):
                raise ValueError(
                    "Continuum rate coefficients must contain one column per shell."
                )
            species_frame.columns = shell_index
            rate_frames.append(species_frame)
        (
            photoionization,
            stimulated_recombination,
            spontaneous_recombination,
            collisional_ionization,
        ) = rate_frames
        spontaneous_recombination *= species_saha_factor
        stimulated_recombination *= species_saha_factor
        collisional_recombination = collisional_ionization * species_saha_factor

        continuum_rate_coefficients_per_shell = tuple(
            ContinuumRateCoefficients(
                photoionization.iloc[:, shell_idx].to_numpy(dtype=np.float64),
                collisional_ionization.iloc[:, shell_idx].to_numpy(
                    dtype=np.float64
                ),
                spontaneous_recombination.iloc[:, shell_idx].to_numpy(
                    dtype=np.float64
                ),
                stimulated_recombination.iloc[:, shell_idx].to_numpy(
                    dtype=np.float64
                ),
                collisional_recombination.iloc[:, shell_idx].to_numpy(
                    dtype=np.float64
                ),
            )
            for shell_idx in range(len(shell_index))
        )

        return (
            continuum_rate_coefficients_per_shell,
            level_to_continuum_saha_factor,
            collisional_ionization_rate_coefficient,
            thermal_saha_factor,
        )

    def _calculate_level_state(
        self,
        shell_idx: int,
        electron_density: float,
        electron_temperature: float,
        level_fractions: npt.NDArray[np.float64],
        continuum_rates: ContinuumRateCoefficients,
        population_geometry: ShellNumberDensity,
        sobolev_inputs: SobolevInputs,
        level_density: npt.NDArray[np.float64] | None = None,
    ) -> tuple[npt.NDArray[np.float64], npt.NDArray[np.float64], float]:
        """Evaluate one level residual at an explicit electron density."""
        total_ionization_rates = continuum_rates.photoionization + (
            continuum_rates.collisional_ionization * electron_density
        )
        total_recombination_rates = (
            continuum_rates.spontaneous_recombination
            + continuum_rates.stimulated_recombination
            + continuum_rates.collisional_recombination * electron_density
        ) * electron_density
        level_rates = LevelEquationRates(
            total_ionization_rates,
            total_recombination_rates,
            -np.diag(total_ionization_rates),
        )
        shell_j_blues = self.j_blues.iloc[:, [shell_idx]].copy()
        shell_j_blues.columns = pd.RangeIndex(1)
        electron_distribution = ThermalElectronEnergyDistribution(
            0.0 * u.erg,
            np.array([electron_temperature]) * u.K,
            np.array([electron_density]) / u.cm**3,
        )
        residual, beta_sobolev, ionized_to_neutral_ratio = (
            calculate_nlte_level_population_residual(
                level_fractions,
                level_rates,
                self.rate_matrix_solver,
                shell_j_blues,
                electron_distribution,
                self.hydrogen_species,
                population_geometry,
                sobolev_inputs,
                level_density,
            )
        )
        return residual, beta_sobolev, ionized_to_neutral_ratio

    def _calculate_level_solution(
        self,
        shell_idx: int,
        electron_density: float,
        electron_temperature: float,
        level_seed: npt.NDArray[np.float64],
        continuum_rates: ContinuumRateCoefficients,
        population_geometry: ShellNumberDensity,
        sobolev_inputs: SobolevInputs,
    ) -> tuple[
        npt.NDArray[np.float64],
        float,
        npt.NDArray[np.float64],
        npt.NDArray[np.float64],
    ]:
        """Solve one shell's normalized reduced level equations."""

        def residual(
            level_fractions: npt.NDArray[np.float64],
        ) -> npt.NDArray[np.float64]:
            return self._calculate_level_state(
                shell_idx,
                electron_density,
                electron_temperature,
                level_fractions,
                continuum_rates,
                population_geometry,
                sobolev_inputs,
            )[0]

        solution = root(residual, level_seed, options={"xtol": 1e-12})
        if not solution.success or not np.isfinite(solution.x).all():
            raise ValueError("Reduced level-population solve did not converge.")
        if np.any(solution.x < 0.0):
            raise ValueError(
                "Reduced level-population solve returned negative fractions."
            )
        if solution.x.sum() <= 0.0:
            raise ValueError(
                "Reduced level-population solve returned zero fractions."
            )
        fractions = solution.x / solution.x.sum()
        level_residual, beta_sobolev, ionized_to_neutral_ratio = (
            self._calculate_level_state(
                shell_idx,
                electron_density,
                electron_temperature,
                fractions,
                continuum_rates,
                population_geometry,
                sobolev_inputs,
            )
        )
        if ionized_to_neutral_ratio < 0.0:
            raise ValueError(
                "Reduced level-population solve returned a negative "
                "ionized-to-neutral ratio."
            )
        return (
            fractions,
            ionized_to_neutral_ratio,
            beta_sobolev,
            level_residual,
        )

    def _hydrogen_index(self) -> pd.MultiIndex:
        """Return the level index for the reduced NLTE species."""
        level_index = self.level_population_index
        return level_index[
            (
                level_index.get_level_values("atomic_number")
                == self.hydrogen_species[0]
            )
            & (
                level_index.get_level_values("ion_number")
                == self.hydrogen_species[1]
            )
        ]

    def _solve_levels(
        self,
        trial_density: pd.Series,
        electron_temperature: npt.NDArray[np.float64],
        level_seed: pd.DataFrame,
        continuum_rate_coefficients: tuple[ContinuumRateCoefficients, ...],
        level_to_continuum_saha_factor: pd.DataFrame,
        collisional_ionization_rate_coefficient: pd.DataFrame,
        lte_ionization_factor: pd.DataFrame,
    ) -> _LevelSolveState:
        """Solve all shell-local reduced level equations."""
        fractions: list[npt.NDArray[np.float64]] = []
        ionized_to_neutral_ratios: list[float] = []
        trial_escape_probabilities: list[npt.NDArray[np.float64]] = []
        trial_residuals: list[npt.NDArray[np.float64]] = []
        full_population = pd.DataFrame(
            np.column_stack(
                [
                    geometry.level_number_density
                    for geometry in self.shell_number_densities
                ]
            ),
            index=self.level_population_index,
            columns=self.elemental_number_density.columns,
        )
        for shell_idx, (
            continuum_rates,
            population_geometry,
            sobolev_inputs,
        ) in enumerate(
            zip(
                continuum_rate_coefficients,
                self.shell_number_densities,
                self.sobolev_inputs,
                strict=True,
            )
        ):
            if len(level_seed) == len(self.level_population_index):
                seed = level_seed.iloc[
                    population_geometry.species_level_positions, shell_idx
                ].to_numpy(dtype=float)
            else:
                seed = level_seed.iloc[:, shell_idx].to_numpy(dtype=float)
            seed = seed / seed.sum()
            shell_result = self._calculate_level_solution(
                shell_idx,
                trial_density.iloc[shell_idx],
                electron_temperature[shell_idx],
                seed,
                continuum_rates,
                population_geometry,
                sobolev_inputs,
            )
            fractions.append(shell_result[0])
            ionized_to_neutral_ratios.append(shell_result[1])
            trial_escape_probabilities.append(shell_result[2])
            trial_residuals.append(shell_result[3])

        hydrogen_index = self._hydrogen_index()
        columns = self.elemental_number_density.columns
        return _LevelSolveState(
            continuum_rate_coefficients,
            level_to_continuum_saha_factor,
            collisional_ionization_rate_coefficient,
            lte_ionization_factor,
            tuple(fractions),
            tuple(ionized_to_neutral_ratios),
            full_population,
            pd.DataFrame(
                np.asarray(fractions).T, index=hydrogen_index, columns=columns
            ),
            pd.DataFrame(
                np.asarray(trial_escape_probabilities).T,
                index=self.sobolev_inputs[0].line_index,
                columns=columns,
            ),
            pd.DataFrame(
                np.asarray(trial_residuals).T,
                index=hydrogen_index,
                columns=columns,
            ),
            hydrogen_index,
        )

    def _solve_charge_population(
        self,
        state: _LevelSolveState,
        trial_electron_distribution: ThermalElectronEnergyDistribution,
        ion_population_arguments: dict[str, object],
    ) -> tuple[pd.DataFrame | None, pd.Series | None]:
        """Solve optional charge conservation from the trial level state."""
        if self.ion_population_solver is None:
            return None, None

        estimated_levels = state.full_population.copy(deep=True)
        for shell_idx, shell_number_density in enumerate(
            self.shell_number_densities
        ):
            temporary_ion = shell_number_density.hydrogen_number_density / (
                1.0 + state.ionized_to_neutral_ratios[shell_idx]
            )
            positions = shell_number_density.species_level_positions
            estimated_levels.iloc[positions, shell_idx] = (
                temporary_ion * state.fractions[shell_idx]
            )
        solver_arguments = deepcopy(ion_population_arguments)
        solver_arguments.update(
            thermal_electron_energy_distribution=trial_electron_distribution,
            estimated_level_population=estimated_levels,
            level_to_continuum_saha_factor=(
                state.level_to_continuum_saha_factor
            ),
        )
        return self.ion_population_solver.solve(
            lte_ionization_factor=state.lte_ionization_factor,
            **solver_arguments,
        )

    def _build_absolute_levels(
        self,
        state: _LevelSolveState,
        ion_population: pd.DataFrame | None,
    ) -> pd.DataFrame:
        """Build absolute level populations from normalized shell solutions."""
        absolute_levels = state.full_population.copy(deep=True)
        if ion_population is None:
            for shell_idx, population_geometry in enumerate(
                self.shell_number_densities
            ):
                absolute_levels.iloc[
                    population_geometry.species_level_positions, shell_idx
                ] = (
                    population_geometry.hydrogen_number_density
                    / (1.0 + state.ionized_to_neutral_ratios[shell_idx])
                    * state.fractions[shell_idx]
                )
            return absolute_levels

        level_atomic_numbers = self.level_population_index.get_level_values(
            "atomic_number"
        )
        level_ion_numbers = self.level_population_index.get_level_values(
            "ion_number"
        )
        for (
            atomic_number,
            ion_number,
        ), ion_density in ion_population.iterrows():
            positions = np.flatnonzero(
                (level_atomic_numbers == atomic_number)
                & (level_ion_numbers == ion_number)
            )
            if positions.size == 0:
                continue
            base_ion_density = state.full_population.iloc[positions].sum()
            nonzero_density = base_ion_density != 0.0
            if not nonzero_density.any():
                continue
            base_values = base_ion_density.to_numpy()
            safe_base_values = np.where(
                nonzero_density.to_numpy(), base_values, 1.0
            )
            scaled_values = (
                state.full_population.iloc[positions].to_numpy()
                / safe_base_values
                * ion_density.to_numpy()
            )
            scaled_values[:, ~nonzero_density.to_numpy()] = (
                state.full_population.iloc[positions].to_numpy()[
                    :, ~nonzero_density.to_numpy()
                ]
            )
            absolute_levels.iloc[positions] = scaled_values

        hydrogen_ion_population = ion_population.loc[self.hydrogen_species]
        absolute_levels.iloc[
            self.level_population_index.get_indexer(state.hydrogen_index)
        ] = state.normalized_population.multiply(
            hydrogen_ion_population, axis=1
        ).to_numpy()
        return absolute_levels

    def _calculate_final_opacities(
        self,
        absolute_levels: pd.DataFrame,
        state: _LevelSolveState,
    ) -> tuple[pd.DataFrame, pd.DataFrame]:
        """Calculate terminal Sobolev optical depths and escape probabilities."""
        final_tau = []
        final_escape_probabilities = []
        for shell_idx, sobolev_inputs in enumerate(self.sobolev_inputs):
            tau, beta = calculate_sobolev_opacities_from_level_densities(
                absolute_levels.iloc[:, shell_idx].to_numpy(dtype=float),
                sobolev_inputs.line_indices,
                sobolev_inputs.lines_lower_level_index,
                sobolev_inputs.lines_upper_level_index,
                sobolev_inputs.g_lower,
                sobolev_inputs.g_upper,
                sobolev_inputs.metastable_upper,
                sobolev_inputs.nlte_lines_mask,
                sobolev_inputs.tau_coefficient,
            )
            final_tau.append(tau)
            final_escape_probabilities.append(beta)
        columns = self.elemental_number_density.columns
        line_index = self.sobolev_inputs[0].line_index
        return (
            pd.DataFrame(
                np.asarray(final_tau).T, index=line_index, columns=columns
            ),
            pd.DataFrame(
                np.asarray(final_escape_probabilities).T,
                index=line_index,
                columns=columns,
            ),
        )

    def _calculate_final_residual(
        self,
        final_density: pd.Series,
        electron_temperature: npt.NDArray[np.float64],
        absolute_levels: pd.DataFrame,
        state: _LevelSolveState,
    ) -> pd.DataFrame:
        """Re-evaluate reduced equations at the terminal density."""
        residuals = []
        for shell_idx, (
            continuum_rates,
            population_geometry,
            sobolev_inputs,
        ) in enumerate(
            zip(
                state.continuum_rate_coefficients,
                self.shell_number_densities,
                self.sobolev_inputs,
                strict=True,
            )
        ):
            residuals.append(
                self._calculate_level_state(
                    shell_idx,
                    final_density.iloc[shell_idx],
                    electron_temperature[shell_idx],
                    state.fractions[shell_idx],
                    continuum_rates,
                    population_geometry,
                    sobolev_inputs,
                    absolute_levels.iloc[:, shell_idx].to_numpy(dtype=float),
                )[0]
            )
        return pd.DataFrame(
            np.asarray(residuals).T,
            index=state.hydrogen_index,
            columns=self.elemental_number_density.columns,
        )

    def _calculate_heating(
        self,
        final_electron_distribution: ThermalElectronEnergyDistribution,
        absolute_levels: pd.DataFrame,
        ion_population: pd.DataFrame | None,
        state: _LevelSolveState,
    ) -> tuple[pd.Series | None, pd.Series | None]:
        """Calculate optional thermal-balance output."""
        if self.thermal_balance_solver is None:
            return None, None
        solver_arguments = deepcopy(self.thermal_balance_arguments)
        electron_rate_solver = self.rate_matrix_solver.electron_rate_solver
        collisional_rates = electron_rate_solver.solve(
            final_electron_distribution.temperature
        )
        transition_count = len(collisional_rates) // 2
        if transition_count:
            transition_index = (
                electron_rate_solver.all_collisional_strengths_index
            )
            collisional_excitation = collisional_rates.iloc[
                :transition_count
            ].copy()
            collisional_deexcitation = collisional_rates.iloc[
                transition_count:
            ].copy()
            collisional_excitation.index = transition_index
            collisional_deexcitation.index = transition_index
            collisional_excitation.columns = absolute_levels.columns
            collisional_deexcitation.columns = absolute_levels.columns
            if self.reference_electron_temperature is not None:
                reference_rates = electron_rate_solver.solve(
                    self.reference_electron_temperature
                )
                reference_excitation = reference_rates.iloc[
                    :transition_count
                ].copy()
                reference_deexcitation = reference_rates.iloc[
                    transition_count:
                ].copy()
                reference_excitation.index = transition_index
                reference_deexcitation.index = transition_index
                reference_excitation.columns = absolute_levels.columns
                reference_deexcitation.columns = absolute_levels.columns
                fixed_excitation = self.thermal_balance_arguments[
                    "collisional_excitation_rate_coefficient"
                ]
                fixed_deexcitation = self.thermal_balance_arguments[
                    "collisional_deexcitation_rate_coefficient"
                ]
                collisional_excitation = fixed_excitation.multiply(
                    collisional_excitation / reference_excitation
                ).reindex(fixed_excitation.index)
                collisional_deexcitation = fixed_deexcitation.multiply(
                    collisional_deexcitation / reference_deexcitation
                ).reindex(fixed_deexcitation.index)
            solver_arguments.update(
                collisional_excitation_rate_coefficient=(
                    collisional_excitation
                ),
                collisional_deexcitation_rate_coefficient=(
                    collisional_deexcitation
                ),
            )
        solver_arguments.update(
            thermal_electron_distribution=final_electron_distribution,
            level_population=absolute_levels,
            ion_population=ion_population,
            collisional_ionization_rate_coefficient=(
                state.collisional_ionization_rate_coefficient
            ),
            level_population_ratio=state.level_to_continuum_saha_factor,
        )
        solver = deepcopy(self.thermal_balance_solver)
        return solver.solve(**solver_arguments)

    def evaluate(
        self,
        trial_electron_density: npt.ArrayLike,
        electron_temperature: u.Quantity | npt.ArrayLike,
        level_seed: pd.DataFrame,
    ) -> PlasmaEquilibriumEvaluation:
        """Evaluate levels and the existing optional charge/thermal stages."""
        trial_density = pd.Series(
            np.asarray(trial_electron_density, dtype=np.float64),
            index=self.elemental_number_density.columns,
        )
        if np.any(trial_density <= 0.0):
            raise ValueError("Trial electron densities must be positive.")
        if isinstance(electron_temperature, u.Quantity):
            temperatures = np.asarray(
                electron_temperature.to_value("K"), dtype=float
            )
        else:
            temperatures = np.asarray(electron_temperature, dtype=float)
        if temperatures.size != len(trial_density):
            raise ValueError("One electron temperature is required per shell.")
        trial_electron_distribution = ThermalElectronEnergyDistribution(
            0.0 * u.erg,
            temperatures * u.K,
            trial_density.to_numpy() / u.cm**3,
        )
        (
            continuum_rate_coefficients,
            level_to_continuum_saha_factor,
            collisional_ionization_rate_coefficient,
            lte_ionization_factor,
        ) = self._calculate_continuum_rate_coefficients(temperatures)
        level_state = self._solve_levels(
            trial_density,
            temperatures,
            level_seed,
            continuum_rate_coefficients,
            level_to_continuum_saha_factor,
            collisional_ionization_rate_coefficient,
            lte_ionization_factor,
        )
        ion_population, solved_density = self._solve_charge_population(
            level_state,
            trial_electron_distribution,
            self.ion_population_arguments,
        )

        final_density = (
            trial_density.copy(deep=True)
            if solved_density is None
            else pd.Series(
                np.asarray(solved_density, dtype=np.float64),
                index=trial_density.index,
            )
        )
        final_electron_distribution = ThermalElectronEnergyDistribution(
            0.0 * u.erg,
            temperatures * u.K,
            final_density.to_numpy() / u.cm**3,
        )
        absolute_levels = self._build_absolute_levels(
            level_state, ion_population
        )

        charge_residual = None
        electron_residual = None
        if ion_population is not None and solved_density is not None:
            charges = ion_population.index.get_level_values(
                "ion_number"
            ).to_numpy()
            charge_density = ion_population.multiply(charges, axis=0).sum()
            charge_residual = (charge_density - solved_density) / np.where(
                self.maximum_electron_density == 0.0,
                1.0,
                self.maximum_electron_density,
            )
            electron_residual = (solved_density - trial_density) / trial_density
        final_tau, final_escape_probabilities = self._calculate_final_opacities(
            absolute_levels, level_state
        )
        total_heating, fractional_heating = self._calculate_heating(
            final_electron_distribution,
            absolute_levels,
            ion_population,
            level_state,
        )
        final_residual = self._calculate_final_residual(
            final_density, temperatures, absolute_levels, level_state
        )

        return PlasmaEquilibriumEvaluation(
            trial_density,
            level_state.normalized_population,
            pd.Series(
                level_state.ionized_to_neutral_ratios,
                index=trial_density.index,
            ),
            level_state.trial_escape_probability,
            level_state.trial_level_residual,
            solved_density,
            absolute_levels,
            ion_population,
            final_tau,
            final_escape_probabilities,
            final_residual,
            charge_residual,
            electron_residual,
            total_heating,
            fractional_heating,
        )
