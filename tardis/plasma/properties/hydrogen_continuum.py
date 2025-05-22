import logging

import astropy.units as u

from tardis.plasma.electron_energy_distribution import (
    ThermalElectronEnergyDistribution,
)
from tardis.plasma.equilibrium.ion_populations import IonPopulationSolver
from tardis.plasma.equilibrium.rate_matrix import IonRateMatrix
from tardis.plasma.equilibrium.rates import (
    AnalyticPhotoionizationRateSolver,
    CollisionalIonizationRateSolver,
)
from tardis.plasma.properties.base import ProcessingPlasmaProperty
from tardis.plasma.properties.ion_population import IonNumberDensity
from tardis.plasma.properties.level_population import (
    LevelNumberDensity,
)

logger = logging.getLogger(__name__)

__all__ = ["HydrogenContinuumPropertes"]


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
        ion_threshold,
    ):
        return self.calculate_with_n_electron(
            thermal_phi_lte,
            thermal_lte_partition_function,
            number_density,
            electron_densities,
            block_ids,
            ion_threshold,
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


class HydrogenContinuumPropertes(ProcessingPlasmaProperty):
    outputs = ("ion_number_density", "electron_densities")

    def __init__(self, plasma_parent, photo_ion_cross_sections):
        super().__init__(plasma_parent)
        self._update_inputs()
        self.photoionization_data = photo_ion_cross_sections

    def calculate(
        self,
        t_electrons,
        previous_electron_densities,
        dilute_planckian_radiation_field,
        number_density,
        lte_level_number_density,
        level_number_density,
        lte_ion_number_density,
        ion_number_density,
    ):
        electron_dist = ThermalElectronEnergyDistribution(
            0 * u.erg,
            t_electrons * u.K,
            previous_electron_densities * u.g / u.cm**3,
        )

        photoionization_rate_solver = AnalyticPhotoionizationRateSolver(
            self.photoionization_data
        )

        collisional_rate_solver = CollisionalIonizationRateSolver(
            self.photoionization_data
        )

        ion_rate_matrix_solver = IonRateMatrix(
            photoionization_rate_solver, collisional_rate_solver
        )

        solver = IonPopulationSolver(ion_rate_matrix_solver, [(1, 0), (1, 1)])

        fractional_ion_population, fractional_electron_density = solver.solve(
            dilute_planckian_radiation_field,
            electron_dist,
            lte_level_number_density,
            level_number_density,
            lte_ion_number_density,
            ion_number_density,
            charge_conservation=True,
        )

        ion_number_density = fractional_ion_population * number_density
        electron_densities = (
            fractional_electron_density * previous_electron_densities
        )

        return ion_number_density, electron_densities
