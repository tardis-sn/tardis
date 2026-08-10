import logging

import astropy.units as u
import numpy.typing as npt
import pandas as pd

from tardis.plasma.electron_energy_distribution import (
    ThermalElectronEnergyDistribution,
)
from tardis.plasma.equilibrium.charge_conservation import (
    ChargeConservationSolver,
)
from tardis.plasma.equilibrium.continuum import ContinuumPopulationSolver
from tardis.plasma.equilibrium.ion_populations import (
    calculate_level_to_ion_population_factor,
    get_lower_ion_level_index,
    get_upper_ion_population_index,
)
from tardis.plasma.equilibrium.rates import (
    CollisionalIonizationSeaton,
    ThermalCollisionalRateSolver,
)
from tardis.plasma.equilibrium.rates.heating_cooling_rates import (
    BoundFreeThermalRates,
    CollisionalBoundThermalRates,
    CollisionalIonizationThermalRates,
    FreeFreeThermalRates,
)
from tardis.plasma.equilibrium.sobolev_solver import (
    SobolevPopulationSolver,
)
from tardis.plasma.equilibrium.thermal_balance import ThermalBalanceSolver
from tardis.plasma.properties.base import (
    Input,
    PreviousIterationProperty,
    ProcessingPlasmaProperty,
)
from tardis.plasma.properties.ion_population import (
    IonNumberDensity,
    ThermalPhiSahaLTE,
    calculate_block_ids_from_dataframe,
)
from tardis.plasma.properties.level_population import (
    LevelNumberDensity,
)
from tardis.plasma.radiation_field import EstimatedRadiationField

logger = logging.getLogger(__name__)

__all__ = [
    "BoundFreeHeatingEstimator",
    "CollisionalDeexcitationRateCoefficient",
    "CollisionalExcitationRateCoefficient",
    "CollisionalIonizationRateCoefficient",
    "FreeFreeHeatingEstimator",
    "HydrogenContinuumFractionalHeating",
    "HydrogenContinuumIonRatio",
    "HydrogenContinuumLTEPopulations",
    "IIPContinuumPopulations",
    "Iteration",
    "PhotoionizationRateEstimator",
    "PreviousIonNumberDensity",
    "PreviousLevelNumberDensity",
    "StimulatedRecombinationCoolingEstimator",
    "StimulatedRecombinationRateEstimator",
]


class Iteration(Input):
    """Store the current continuum iteration number.

    Attributes
    ----------
    iteration : int
        Current iteration number.
    """

    outputs = ("iteration",)


class CollisionalIonizationRateCoefficient(Input):
    """Store collisional-ionization rate coefficients.

    Attributes
    ----------
    collisional_ionization_rate_coefficient : pd.DataFrame
        Collisional ionization rate coefficients from estimators.
    """

    outputs = ("collisional_ionization_rate_coefficient",)


class CollisionalExcitationRateCoefficient(Input):
    """Store collisional-excitation rate coefficients.

    Attributes
    ----------
    collisional_excitation_rate_coefficient : pd.DataFrame
        Collisional excitation rate coefficients from estimators.
    """

    outputs = ("collisional_excitation_rate_coefficient",)


class CollisionalDeexcitationRateCoefficient(Input):
    """Store collisional-deexcitation rate coefficients.

    Attributes
    ----------
    collisional_deexcitation_rate_coefficient : pd.DataFrame
        Collisional de-excitation rate coefficients from estimators.
    """

    outputs = ("collisional_deexcitation_rate_coefficient",)


class PhotoionizationRateEstimator(Input):
    """Store Monte Carlo photoionization-rate estimators.

    Attributes
    ----------
    photoionization_rate_estimator : pd.DataFrame
        Photoionization rate estimator.
    """

    outputs = ("photoionization_rate_estimator",)


class StimulatedRecombinationRateEstimator(Input):
    """Store Monte Carlo stimulated-recombination estimators.

    Attributes
    ----------
    stimulated_recombination_rate_estimator : pd.DataFrame
        Stimulated recombination rate estimator.
    """

    outputs = ("stimulated_recombination_rate_estimator",)


class FreeFreeHeatingEstimator(Input):
    """Store Monte Carlo free-free heating estimators.

    Attributes
    ----------
    free_free_heating_estimator : pandas.DataFrame
        Free-free heating estimator from the Monte Carlo transport.
    """

    outputs = ("free_free_heating_estimator",)


class BoundFreeHeatingEstimator(Input):
    """Store Monte Carlo bound-free heating estimators.

    Attributes
    ----------
    bound_free_heating_estimator : pandas.DataFrame
        Bound-free heating estimator from the Monte Carlo transport.
    """

    outputs = ("bound_free_heating_estimator",)


class StimulatedRecombinationCoolingEstimator(Input):
    """Store Monte Carlo stimulated-recombination cooling estimators.

    Attributes
    ----------
    stimulated_recombination_cooling_estimator : pandas.DataFrame
        Stimulated recombination cooling estimator from Monte Carlo transport.
    """

    outputs = ("stimulated_recombination_cooling_estimator",)


class PreviousIonNumberDensity(PreviousIterationProperty):
    """Store ion populations from the previous plasma iteration.

    Attributes
    ----------
    previous_ion_number_density : pd.DataFrame
        Previous iteration ion number density.
    """

    outputs = ("previous_ion_number_density",)

    def set_initial_value(self, kwargs: dict[str, object]) -> None:
        """Initialize the previous ion population as unavailable."""
        self._set_initial_value(None)


class PreviousLevelNumberDensity(PreviousIterationProperty):
    """Store level populations from the previous plasma iteration.

    Attributes
    ----------
    previous_level_number_density : pd.DataFrame
        Previous iteration level number density.
    """

    outputs = ("previous_level_number_density",)

    def set_initial_value(self, kwargs: dict[str, object]) -> None:
        """Initialize the previous level population as unavailable."""
        self._set_initial_value(None)


class LTEIonNumberDensity(IonNumberDensity):
    outputs = ("lte_ion_number_density",)
    latex_name = ("N_{i,j}^*",)

    def calculate(
        self,
        thermal_phi_lte,
        thermal_lte_partition_function,
        number_density,
        electron_densities,
        block_ids,
    ):
        return self.calculate_with_n_electron(
            thermal_phi_lte,
            thermal_lte_partition_function,
            number_density,
            electron_densities,
            block_ids,
            1e-20,  # ION_ZERO_THRESHOLD from ion_population.py
        )


class LTELevelNumberDensity(LevelNumberDensity):
    outputs = ("lte_level_number_density",)
    latex_name = ("N_{i,j,k}^*",)

    def _calculate_dilute_lte(
        self,
        thermal_lte_level_boltzmann_factor,
        lte_ion_number_density,
        levels,
        thermal_lte_partition_function,
    ):
        return super()._calculate_dilute_lte(
            thermal_lte_level_boltzmann_factor,
            lte_ion_number_density,
            levels,
            thermal_lte_partition_function,
        )


class IIPContinuumPopulations(ProcessingPlasmaProperty):
    """Calculate the configured-continuum population state in one step."""

    outputs = (
        "hydrogen_continuum_level_boltzmann_factor",
        "hydrogen_continuum_partition_function",
        "ion_number_density",
        "electron_densities",
        "level_number_density",
    )

    @staticmethod
    def _partition_function(
        level_boltzmann_factor: pd.DataFrame,
    ) -> pd.DataFrame:
        return level_boltzmann_factor.groupby(
            level=["atomic_number", "ion_number"]
        ).sum()

    def calculate(
        self,
        atomic_data: object,
        nlte_data: object,
        continuum_interaction_species: pd.MultiIndex,
        t_electrons: npt.ArrayLike,
        dilute_planckian_radiation_field: object,
        j_blues: pd.DataFrame | None,
        previous_beta_sobolev: pd.DataFrame | None,
        previous_electron_densities: pd.Series | None,
        previous_level_number_density: pd.DataFrame | None,
        previous_ion_number_density: pd.DataFrame | None,
        number_density: pd.DataFrame,
        general_level_boltzmann_factor: pd.DataFrame,
        thermal_g_electron: npt.ArrayLike,
        beta_electron: npt.ArrayLike,
        thermal_lte_level_boltzmann_factor: pd.DataFrame,
        thermal_lte_partition_function: pd.DataFrame,
        ionization_data: pd.DataFrame,
        phi: pd.DataFrame,
        g: pd.Series,
        lines_lower_level_index: npt.ArrayLike,
        lines_upper_level_index: npt.ArrayLike,
        metastability: pd.Series,
        lines: pd.DataFrame,
        time_explosion: u.Quantity,
        photoionization_rate_estimator: object | None,
        stimulated_recombination_rate_estimator: object | None,
    ) -> tuple[
        pd.DataFrame,
        pd.DataFrame,
        pd.DataFrame,
        pd.Series,
        pd.DataFrame,
    ]:
        """Solve and publish the coherent continuum population state.

        Parameters
        ----------
        atomic_data : object
            Atomic data used to assemble bound-bound and bound-free rates.
        nlte_data : object
            Configured NLTE species and line-population policy.
        continuum_interaction_species : pandas.MultiIndex
            Ion stages with explicit continuum level blocks.
        t_electrons : array-like
            Electron temperatures for every shell.
        dilute_planckian_radiation_field : object
            Radiation field used before estimator-backed updates.
        j_blues : pandas.DataFrame, optional
            Line mean intensities for estimator-backed updates.
        previous_beta_sobolev : pandas.DataFrame
            Escape-probability warm start.
        previous_electron_densities : pandas.Series
            Previous electron-density state retained by the property graph.
        previous_level_number_density : pandas.DataFrame
            Level-population warm start for the Sobolev iteration.
        previous_ion_number_density : pandas.DataFrame
            Previous ion state retained by the property graph; it is not a
            frozen contribution to charge conservation.
        number_density : pandas.DataFrame
            Elemental number densities indexed by atomic number.
        general_level_boltzmann_factor : pandas.DataFrame
            Baseline level factors for all ions.
        thermal_g_electron, beta_electron : array-like
            Thermal electron statistical weight and inverse temperature.
        thermal_lte_level_boltzmann_factor : pandas.DataFrame
            Thermal LTE level factors used for detailed balance.
        thermal_lte_partition_function : pandas.DataFrame
            Thermal LTE partition functions.
        ionization_data : pandas.DataFrame
            Ionization energies indexed by ion stage.
        phi : pandas.DataFrame
            Nebular ionization factors for uncovered ionization edges.
        g : pandas.Series
            Statistical weights for the atomic levels.
        lines_lower_level_index, lines_upper_level_index : array-like
            Positional lower- and upper-level indexes for each line.
        metastability : pandas.Series
            Metastability flags for the atomic levels.
        lines : pandas.DataFrame
            Atomic line data.
        time_explosion : astropy.units.Quantity
            Time since explosion.
        photoionization_rate_estimator, stimulated_recombination_rate_estimator : object, optional
            Paired Monte Carlo bound-free estimators.

        Returns
        -------
        tuple[pandas.DataFrame, pandas.DataFrame, pandas.DataFrame, pandas.Series, pandas.DataFrame]
            Level factors, partition functions, ion populations, electron
            densities, and level populations from one coherent state.
        """
        level_boltzmann_factor = general_level_boltzmann_factor.copy(deep=True)
        if (photoionization_rate_estimator is None) != (
            stimulated_recombination_rate_estimator is None
        ):
            raise ValueError(
                "Photoionization and stimulated-recombination estimators must be supplied together"
            )
        general_partition_function = self._partition_function(
            level_boltzmann_factor
        )
        if photoionization_rate_estimator is None:
            population_solver = ContinuumPopulationSolver(
                atomic_data=atomic_data,
                lines=lines,
                continuum_species=(),
                radiation_field=dilute_planckian_radiation_field,
                electron_temperatures=t_electrons,
                elemental_number_density=number_density,
                general_level_boltzmann_factor=level_boltzmann_factor,
                general_partition_function=general_partition_function,
                thermal_g_electron=thermal_g_electron,
                beta_electron=beta_electron,
                thermal_lte_level_boltzmann_factor=(
                    thermal_lte_level_boltzmann_factor
                ),
                thermal_lte_partition_function=(thermal_lte_partition_function),
                ionization_data=ionization_data,
                nebular_phi=phi,
                photoionization_rate_estimator=None,
                stimulated_recombination_rate_estimator=None,
            )
            population_state = ChargeConservationSolver(
                number_density, population_solver
            ).solve(
                electron_number_density_seed=(
                    None
                    if previous_electron_densities is None
                    else previous_electron_densities.copy(deep=True)
                )
            )
            updated_partition_function = self._partition_function(
                population_state.level_boltzmann_factor
            )
            return (
                population_state.level_boltzmann_factor,
                updated_partition_function,
                population_state.ion_number_density,
                population_state.electron_densities,
                population_state.level_number_density,
            )
        if j_blues is None:
            raise ValueError(
                "Estimator-backed continuum updates require j_blues"
            )
        if any(
            value is None
            for value in (
                previous_beta_sobolev,
                previous_electron_densities,
                previous_level_number_density,
                previous_ion_number_density,
            )
        ):
            raise ValueError(
                "Estimator-backed continuum updates require a complete previous population state"
            )

        radiation_field = EstimatedRadiationField(
            pd.DataFrame(
                j_blues,
                index=lines.index,
                columns=number_density.columns,
            )
        )

        population_solver = ContinuumPopulationSolver(
            atomic_data=atomic_data,
            lines=lines,
            continuum_species=continuum_interaction_species,
            radiation_field=radiation_field,
            electron_temperatures=t_electrons,
            elemental_number_density=number_density,
            general_level_boltzmann_factor=level_boltzmann_factor,
            general_partition_function=general_partition_function,
            thermal_g_electron=thermal_g_electron,
            beta_electron=beta_electron,
            thermal_lte_level_boltzmann_factor=(
                thermal_lte_level_boltzmann_factor
            ),
            thermal_lte_partition_function=(thermal_lte_partition_function),
            ionization_data=ionization_data,
            nebular_phi=phi,
            photoionization_rate_estimator=(photoionization_rate_estimator),
            stimulated_recombination_rate_estimator=(
                stimulated_recombination_rate_estimator
            ),
        )
        charge_conservation_solver = ChargeConservationSolver(
            number_density, population_solver
        )
        sobolev_solver = SobolevPopulationSolver(
            charge_conservation_solver,
            lines,
            time_explosion,
            g,
            lines_lower_level_index,
            lines_upper_level_index,
            metastability,
            nlte_species=set(nlte_data.nlte_species),
        )
        population_state, _, _, _ = sobolev_solver.solve(
            previous_level_number_density.copy(deep=True),
            previous_beta_sobolev,
            previous_electron_densities.copy(deep=True),
        )
        partition_function = self._partition_function(
            population_state.level_boltzmann_factor
        )
        return (
            population_state.level_boltzmann_factor,
            partition_function,
            population_state.ion_number_density,
            population_state.electron_densities,
            population_state.level_number_density,
        )


class HydrogenContinuumLTEPopulations(ProcessingPlasmaProperty):
    """
    Calculate LTE reference populations for the hydrogen continuum plasma.

    Attributes
    ----------
    lte_ion_number_density : pandas.DataFrame
        LTE ion number density indexed by atomic number and ion number with
        columns corresponding to zones.
    lte_level_number_density : pandas.DataFrame
        LTE level number density indexed by atomic number, ion number, and
        level number with columns corresponding to zones.
    """

    outputs = ("lte_ion_number_density", "lte_level_number_density")

    def calculate(
        self,
        number_density: pd.DataFrame,
        electron_densities: pd.Series,
        levels: pd.MultiIndex,
        hydrogen_continuum_partition_function: pd.DataFrame,
        thermal_g_electron: npt.ArrayLike,
        beta_electron: npt.ArrayLike,
        ionization_data: pd.DataFrame,
        thermal_lte_level_boltzmann_factor: pd.DataFrame,
        thermal_lte_partition_function: pd.DataFrame,
    ) -> tuple[pd.DataFrame, pd.DataFrame]:
        """Calculate LTE reference populations at the current electron density.

        Parameters
        ----------
        number_density : pandas.DataFrame
            Elemental number densities.
        electron_densities : pandas.Series
            Charge-conserving electron densities from the coupled state.
        levels : pandas.MultiIndex
            Atomic level index.
        hydrogen_continuum_partition_function : pandas.DataFrame
            Partition functions from the coupled level factors.
        thermal_g_electron, beta_electron : array-like
            Thermal electron degeneracy and inverse temperature factors.
        ionization_data : pandas.DataFrame
            Ionization energies indexed by ion stage.
        thermal_lte_level_boltzmann_factor : pandas.DataFrame
            Thermal LTE level factors.
        thermal_lte_partition_function : pandas.DataFrame
            Thermal LTE partition functions.

        Returns
        -------
        tuple[pandas.DataFrame, pandas.DataFrame]
            LTE ion and level number densities.
        """
        thermal_phi_lte = ThermalPhiSahaLTE.calculate(
            thermal_g_electron,
            beta_electron,
            thermal_lte_partition_function,
            ionization_data,
        )

        block_ids = calculate_block_ids_from_dataframe(
            hydrogen_continuum_partition_function
        )

        lte_ion_number_density, block_ids = LTEIonNumberDensity(
            self.plasma_parent, electron_densities=electron_densities
        ).calculate(
            thermal_phi_lte,
            thermal_lte_partition_function,
            number_density,
            electron_densities,
            block_ids,
        )

        lte_level_number_density = LTELevelNumberDensity(
            self.plasma_parent
        ).calculate(
            thermal_lte_level_boltzmann_factor,
            lte_ion_number_density,
            levels,
            thermal_lte_partition_function,
        )

        return lte_ion_number_density, lte_level_number_density


class HydrogenContinuumFractionalHeating(ProcessingPlasmaProperty):
    """
    Calculate fractional heating rates for hydrogen continuum iterations.

    Attributes
    ----------
    fractional_heating : pandas.DataFrame or None
        Fractional contribution of thermal processes to hydrogen continuum
        heating in each zone. The value is ``None`` during the first iteration.
    """

    outputs = ("fractional_heating",)

    def __init__(self, plasma_parent, photo_ion_cross_sections):
        super().__init__(plasma_parent)
        self._update_inputs()
        self.photoionization_data = photo_ion_cross_sections

    def calculate(
        self,
        atomic_data,
        lines,
        t_electrons,
        dilute_planckian_radiation_field,
        lte_level_number_density,
        lte_ion_number_density,
        level_number_density,
        ion_number_density,
        electron_densities,
        iteration,
        collisional_ionization_rate_coefficient,
        collisional_deexcitation_rate_coefficient,
        collisional_excitation_rate_coefficient,
        free_free_heating_estimator,
        bound_free_heating_estimator,
        stimulated_recombination_cooling_estimator,
    ):
        """Calculate the fractional heating for one continuum update."""
        if iteration == 0:
            return None

        if collisional_ionization_rate_coefficient is None:
            collisional_ionization_rate_coefficient = (
                CollisionalIonizationSeaton(self.photoionization_data).solve(
                    t_electrons * u.K
                )
            )
        if (
            collisional_excitation_rate_coefficient is None
            or collisional_deexcitation_rate_coefficient is None
        ):
            thermal_collisional_rates = ThermalCollisionalRateSolver(
                atomic_data.levels,
                lines,
                atomic_data.collision_data_temperatures,
                atomic_data.yg_data,
                "cmfgen",
            ).solve(t_electrons * u.K)
            source = thermal_collisional_rates.index.get_level_values(
                "level_number_source"
            )
            destination = thermal_collisional_rates.index.get_level_values(
                "level_number_destination"
            )
            collisional_excitation_rate_coefficient = thermal_collisional_rates[
                source < destination
            ]
            collisional_deexcitation_rate_coefficient = (
                thermal_collisional_rates[source >= destination]
            )

        collisional_excitation_rate_coefficient.index = (
            collisional_excitation_rate_coefficient.index.rename(
                {
                    "level_number_source": "level_number_lower",
                    "level_number_destination": "level_number_upper",
                }
            )
        )
        collisional_deexcitation_rate_coefficient.index = (
            collisional_deexcitation_rate_coefficient.index.rename(
                {
                    "level_number_source": "level_number_lower",
                    "level_number_destination": "level_number_upper",
                }
            )
        )

        bound_free_thermal_rates = BoundFreeThermalRates(
            self.photoionization_data
        )
        free_free_thermal_rates = FreeFreeThermalRates()
        collisional_ionization_thermal_rates = (
            CollisionalIonizationThermalRates(self.photoionization_data)
        )
        collisional_bound_thermal_rates = CollisionalBoundThermalRates(lines)

        electron_distribution = ThermalElectronEnergyDistribution(
            0,
            t_electrons * u.K,
            electron_densities.values * u.cm**-3,
        )

        # TODO: make more general indices that work for non-Hydrogen species.
        lower_ion_level_index = get_lower_ion_level_index(
            lte_level_number_density
        )
        upper_ion_population_index = get_upper_ion_population_index(
            lte_ion_number_density
        )

        level_population_ratio = calculate_level_to_ion_population_factor(
            lte_level_number_density.loc[lower_ion_level_index],
            lte_ion_number_density.loc[upper_ion_population_index],
            electron_distribution.number_density,
        )

        _total_heating_rate, fractional_heating_rate = ThermalBalanceSolver(
            bound_free_solver=bound_free_thermal_rates,
            free_free_solver=free_free_thermal_rates,
            collisional_ionization_solver=collisional_ionization_thermal_rates,
            collisional_bound_solver=collisional_bound_thermal_rates,
        ).solve(
            electron_distribution,
            level_number_density,
            ion_number_density,
            collisional_ionization_rate_coefficient,
            collisional_deexcitation_rate_coefficient,
            collisional_excitation_rate_coefficient,
            free_free_heating_estimator,
            level_population_ratio,
            dilute_planckian_radiation_field,
            bound_free_heating_estimator,
            stimulated_recombination_cooling_estimator,
        )

        return fractional_heating_rate


class HydrogenContinuumIonRatio(ProcessingPlasmaProperty):
    """
    Calculate the hydrogen ionization ratio.

    Attributes
    ----------
    ion_ratio : pandas.Series
        Ratio of H II to H I number density in each zone.
    """

    outputs = ("ion_ratio",)

    @staticmethod
    def calculate(ion_number_density):
        """Calculate the H II to H I population ratio."""
        ion_ratio = (
            ion_number_density.loc[(1, 1)] / ion_number_density.loc[(1, 0)]
        )
        return ion_ratio
