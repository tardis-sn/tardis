import logging

import astropy.units as u

from tardis.plasma.electron_energy_distribution import (
    ThermalElectronEnergyDistribution,
)
from tardis.plasma.equilibrium.ion_populations import (
    AnalyticEquilibriumIonPopulationSolver,
    EstimatedEquilibriumIonPopulationSolver,
    LTEIonPopulationSolver,
    calculate_level_to_ion_population_factor,
    get_lower_ion_level_index,
    get_upper_ion_population_index,
)
from tardis.plasma.equilibrium.level_populations import LevelPopulationSolver
from tardis.plasma.equilibrium.rate_matrix import (
    EquilibriumIonRateMatrix,
    LevelRateMatrix,
    LTEIonRateMatrix,
)
from tardis.plasma.equilibrium.rates import (
    AnalyticPhotoionizationRateSolver,
    CollisionalIonizationRateSolver,
    CollisionalIonizationSeaton,
    EstimatedPhotoionizationRateSolver,
    RadiativeRatesSolver,
    ThermalCollisionalRateSolver,
)
from tardis.plasma.equilibrium.rates.heating_cooling_rates import (
    BoundFreeThermalRates,
    CollisionalBoundThermalRates,
    CollisionalIonizationThermalRates,
    FreeFreeThermalRates,
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
    "HydrogenContinuumIonPopulations",
    "HydrogenContinuumIonRatio",
    "HydrogenContinuumLTEPopulations",
    "HydrogenContinuumLevelBoltzmannFactor",
    "HydrogenContinuumLevelNumberDensity",
    "HydrogenContinuumPartitionFunction",
    "Iteration",
    "PhotoionizationRateEstimator",
    "PreviousIonNumberDensity",
    "PreviousLevelNumberDensity",
    "StimulatedRecombinationCoolingEstimator",
    "StimulatedRecombinationRateEstimator",
]


class Iteration(Input):
    """
    Attributes
    ----------
    iteration : int
        Current iteration number.
    """

    outputs = ("iteration",)


class CollisionalIonizationRateCoefficient(Input):
    """
    Attributes
    ----------
    collisional_ionization_rate_coefficient : pd.DataFrame
        Collisional ionization rate coefficients from estimators.
    """

    outputs = ("collisional_ionization_rate_coefficient",)


class CollisionalExcitationRateCoefficient(Input):
    """
    Attributes
    ----------
    collisional_excitation_rate_coefficient : pd.DataFrame
        Collisional excitation rate coefficients from estimators.
    """

    outputs = ("collisional_excitation_rate_coefficient",)


class CollisionalDeexcitationRateCoefficient(Input):
    """
    Attributes
    ----------
    collisional_deexcitation_rate_coefficient : pd.DataFrame
        Collisional de-excitation rate coefficients from estimators.
    """

    outputs = ("collisional_deexcitation_rate_coefficient",)


class PhotoionizationRateEstimator(Input):
    """
    Attributes
    ----------
    photoionization_rate_estimator : pd.DataFrame
        Photoionization rate estimator.
    """

    outputs = ("photoionization_rate_estimator",)


class StimulatedRecombinationRateEstimator(Input):
    """
    Attributes
    ----------
    stimulated_recombination_rate_estimator : pd.DataFrame
        Stimulated recombination rate estimator.
    """

    outputs = ("stimulated_recombination_rate_estimator",)


class FreeFreeHeatingEstimator(Input):
    """
    Attributes
    ----------
    free_free_heating_estimator : pandas.DataFrame
        Free-free heating estimator from the Monte Carlo transport.
    """

    outputs = ("free_free_heating_estimator",)


class BoundFreeHeatingEstimator(Input):
    """
    Attributes
    ----------
    bound_free_heating_estimator : pandas.DataFrame
        Bound-free heating estimator from the Monte Carlo transport.
    """

    outputs = ("bound_free_heating_estimator",)


class StimulatedRecombinationCoolingEstimator(Input):
    """
    Attributes
    ----------
    stimulated_recombination_cooling_estimator : pandas.DataFrame
        Stimulated recombination cooling estimator from Monte Carlo transport.
    """

    outputs = ("stimulated_recombination_cooling_estimator",)


class PreviousIonNumberDensity(PreviousIterationProperty):
    """
    Attributes
    ----------
    previous_ion_number_density : pd.DataFrame
        Previous iteration ion number density.
    """

    outputs = ("previous_ion_number_density",)

    def set_initial_value(self, kwargs):
        self._set_initial_value(None)


class PreviousLevelNumberDensity(PreviousIterationProperty):
    """
    Attributes
    ----------
    previous_level_number_density : pd.DataFrame
        Previous iteration level number density.
    """

    outputs = ("previous_level_number_density",)

    def set_initial_value(self, kwargs):
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


class HydrogenContinuumLevelBoltzmannFactor(ProcessingPlasmaProperty):
    """
    Calculate level Boltzmann factors for the hydrogen continuum plasma.

    Attributes
    ----------
    hydrogen_continuum_level_boltzmann_factor : pandas.DataFrame
        Level population proportionality factors indexed by atomic number,
        ion number, and level number with columns corresponding to zones.
    """

    outputs = ("hydrogen_continuum_level_boltzmann_factor",)

    def __init__(self, plasma_parent):
        super().__init__(plasma_parent)
        self._update_inputs()

    def calculate_equilibrium_level_boltzmann_factor(
        self,
        atomic_data,
        species_list,
        t_electrons,
        dilute_planckian_radiation_field,
        general_level_boltzmann_factor,
        previous_electron_densities,
    ):
        for species in species_list:
            logger.info("Calculating rates for species %s", species)
            species_slice = (species[0], species[1], slice(None), slice(None))
            radiative_transitions = atomic_data.lines.loc[species_slice, :]
            radiative_rate_solver = RadiativeRatesSolver(radiative_transitions)

            col_strength_temperatures = atomic_data.collision_data_temperatures
            col_strengths = atomic_data.yg_data.loc[species_slice, :]

            collisional_rate_solver = ThermalCollisionalRateSolver(
                atomic_data.levels,
                radiative_transitions,
                col_strength_temperatures,
                col_strengths,
                "cmfgen",
            )

            rate_matrix_solver = LevelRateMatrix(atomic_data.levels)

            # A fake electron distribution. Will eventually be a direct input
            # to the plasma property.
            electron_distribution = ThermalElectronEnergyDistribution(
                0 * u.erg,
                t_electrons * u.K,
                previous_electron_densities * u.g / u.cm**3,
            )

            radiative_rates_df = radiative_rate_solver.solve(
                dilute_planckian_radiation_field
            )
            collisional_rates_df = collisional_rate_solver.solve(
                electron_distribution.temperature
            )

            rate_matrix = rate_matrix_solver.solve(
                radiative_rates_df,
                collisional_rates_df,
                electron_distribution.number_density,
            )

            solver = LevelPopulationSolver(rate_matrix, atomic_data.levels)

            level_pops = solver.solve()

            general_level_boltzmann_factor.loc[species] = (
                level_pops.loc[species]
            ).values
        return general_level_boltzmann_factor

    def calculate(
        self,
        atomic_data,
        nlte_data,
        t_electrons,
        dilute_planckian_radiation_field,
        general_level_boltzmann_factor,
        previous_electron_densities,
        j_blues,
    ):
        if j_blues is None:
            return general_level_boltzmann_factor

        estimated_radiation_field = EstimatedRadiationField(j_blues)
        return self.calculate_equilibrium_level_boltzmann_factor(
            atomic_data,
            nlte_data.nlte_species,
            t_electrons,
            estimated_radiation_field,
            general_level_boltzmann_factor,
            previous_electron_densities,
        )


class HydrogenContinuumPartitionFunction(ProcessingPlasmaProperty):
    """
    Calculate partition functions from hydrogen continuum level factors.

    Attributes
    ----------
    hydrogen_continuum_partition_function : pandas.DataFrame
        Partition function indexed by atomic number and ion number with
        columns corresponding to zones.
    """

    outputs = ("hydrogen_continuum_partition_function",)

    @staticmethod
    def calculate(hydrogen_continuum_level_boltzmann_factor):
        return hydrogen_continuum_level_boltzmann_factor.groupby(
            level=["atomic_number", "ion_number"]
        ).sum()


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
        number_density,
        previous_electron_densities,
        levels,
        hydrogen_continuum_partition_function,
        thermal_g_electron,
        beta_electron,
        ionization_data,
        thermal_lte_level_boltzmann_factor,
        thermal_lte_partition_function,
    ):
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
            self.plasma_parent, electron_densities=previous_electron_densities
        ).calculate(
            thermal_phi_lte,
            thermal_lte_partition_function,
            number_density,
            previous_electron_densities,
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


class HydrogenContinuumIonPopulations(ProcessingPlasmaProperty):
    """
    Solve hydrogen continuum ion populations and electron densities.

    Attributes
    ----------
    ion_number_density : pandas.DataFrame
        Hydrogen ion number density indexed by atomic number and ion number
        with columns corresponding to zones.
    electron_densities : pandas.Series
        Electron number density in each zone.
    """

    outputs = ("ion_number_density", "electron_densities")

    def __init__(self, plasma_parent, photo_ion_cross_sections):
        super().__init__(plasma_parent)
        self._update_inputs()
        self.photoionization_data = photo_ion_cross_sections

    def calculate_lte_hydrogen_ion_number_density(
        self,
        t_electrons,
        previous_electron_densities,
        elemental_number_density,
        lte_ion_number_density,
        phi,
        partition_function,
    ):
        electron_dist = ThermalElectronEnergyDistribution(
            0 * u.erg,
            t_electrons * u.K,
            previous_electron_densities.values
            * u.g
            / u.cm**3,  # convert series to quantity
        )

        ion_rate_matrix_solver = LTEIonRateMatrix()

        solver = LTEIonPopulationSolver(
            ion_rate_matrix_solver, elemental_number_density
        )

        ion_population, electron_density = solver.solve(
            electron_dist,
            lte_ion_number_density,
            phi,
            partition_function,
            charge_conservation=False,
        )

        return ion_population, electron_density

    def calculate_analytic_equilibrium_hydrogen_ion_populations(
        self,
        t_electrons,
        previous_electron_densities,
        radiation_field,
        elemental_number_density,
        lte_level_number_density,
        lte_ion_number_density,
        partition_function,
        level_boltzmann_factor,
    ):
        photoionization_rate_solver = AnalyticPhotoionizationRateSolver(
            self.photoionization_data
        )
        collisional_rate_solver = CollisionalIonizationRateSolver(
            self.photoionization_data
        )
        solver = AnalyticEquilibriumIonPopulationSolver(
            EquilibriumIonRateMatrix(),
            photoionization_rate_solver,
            collisional_rate_solver,
            elemental_number_density,
        )
        electron_dist = ThermalElectronEnergyDistribution(
            0 * u.erg,
            t_electrons * u.K,
            previous_electron_densities.values * u.cm**-3,
        )
        return solver.solve(
            radiation_field,
            electron_dist,
            lte_level_number_density,
            lte_level_number_density,
            lte_ion_number_density,
            lte_ion_number_density,
            partition_function,
            level_boltzmann_factor,
            charge_conservation=False,
        )

    def calculate_estimated_equilibrium_hydrogen_ion_populations(
        self,
        t_electrons,
        previous_electron_densities,
        radiation_field_estimators,
        elemental_number_density,
        lte_level_number_density,
        level_number_density,
        lte_ion_number_density,
        ion_number_density,
        partition_function,
        level_boltzmann_factor,
        level2continuum_edge_idx,
    ):
        photoionization_rate_solver = EstimatedPhotoionizationRateSolver(
            self.photoionization_data, level2continuum_edge_idx
        )

        collisional_rate_solver = CollisionalIonizationRateSolver(
            self.photoionization_data
        )

        ion_rate_matrix_solver = EquilibriumIonRateMatrix()

        solver = EstimatedEquilibriumIonPopulationSolver(
            ion_rate_matrix_solver,
            photoionization_rate_solver,
            collisional_rate_solver,
            elemental_number_density,
        )

        electron_dist = ThermalElectronEnergyDistribution(
            0 * u.erg,
            t_electrons * u.K,
            previous_electron_densities.values
            * u.g
            / u.cm**3,  # convert series to quantity
        )

        ion_number_density, electron_densities = solver.solve(
            radiation_field_estimators,
            electron_dist,
            lte_level_number_density,
            level_number_density,
            lte_ion_number_density,
            ion_number_density,
            partition_function,
            level_boltzmann_factor,
            charge_conservation=False,
        )

        return ion_number_density, electron_densities

    def calculate(
        self,
        atomic_data,
        dilute_planckian_radiation_field,
        t_electrons,
        previous_electron_densities,
        previous_level_number_density,
        previous_ion_number_density,
        number_density,
        phi,
        hydrogen_continuum_partition_function,
        lte_level_number_density,
        lte_ion_number_density,
        hydrogen_continuum_level_boltzmann_factor,
        iteration,
        photoionization_rate_estimator,
        stimulated_recombination_rate_estimator,
    ):
        # No-estimator initialization uses dilute-LTE excitation with nebular
        # ionization and charge neutrality. Later iterations use estimators.

        if iteration == 0:
            hydrogen_ion_number_density, _ = IonNumberDensity(None).calculate(
                phi,
                hydrogen_continuum_partition_function,
                number_density,
            )
            charges = hydrogen_ion_number_density.index.get_level_values(
                "ion_number"
            )
            electron_number_density = hydrogen_ion_number_density.multiply(
                charges, axis=0
            ).sum(axis=0)
        else:
            radiation_field_estimators = {
                "photoionization_rate_estimator": photoionization_rate_estimator,
                "stimulated_recombination_rate_estimator": (
                    stimulated_recombination_rate_estimator
                ),
            }

            hydrogen_ion_number_density, electron_number_density = (
                self.calculate_estimated_equilibrium_hydrogen_ion_populations(
                    t_electrons,
                    previous_electron_densities,
                    radiation_field_estimators,
                    number_density,
                    lte_level_number_density,
                    previous_level_number_density,
                    lte_ion_number_density,
                    previous_ion_number_density,
                    hydrogen_continuum_partition_function,
                    hydrogen_continuum_level_boltzmann_factor,
                    atomic_data.level2continuum_edge_idx,
                )
            )

        hydrogen_ion_number_density.columns = number_density.columns
        electron_number_density.index = number_density.columns
        return hydrogen_ion_number_density, electron_number_density


class HydrogenContinuumLevelNumberDensity(ProcessingPlasmaProperty):
    """
    Calculate level populations from hydrogen continuum ion populations.

    Attributes
    ----------
    level_number_density : pandas.DataFrame
        Level number density indexed by atomic number, ion number, and level
        number with columns corresponding to zones.
    """

    outputs = ("level_number_density",)

    def calculate(
        self,
        hydrogen_continuum_level_boltzmann_factor,
        ion_number_density,
        levels,
        hydrogen_continuum_partition_function,
    ):
        return LevelNumberDensity(self.plasma_parent).calculate(
            hydrogen_continuum_level_boltzmann_factor,
            ion_number_density,
            levels,
            hydrogen_continuum_partition_function,
        )


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
        ion_ratio = (
            ion_number_density.loc[(1, 1)] / ion_number_density.loc[(1, 0)]
        )
        return ion_ratio
