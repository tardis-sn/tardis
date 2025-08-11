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
from tardis.plasma.properties.base import (
    Input,
    ProcessingPlasmaProperty,
)
from tardis.plasma.properties.ion_population import (
    IonNumberDensity,
    PhiSahaNebular,
    ThermalPhiSahaLTE,
)
from tardis.plasma.properties.level_population import (
    LevelNumberDensity,
)
from tardis.plasma.properties.partition_function import LevelBoltzmannFactorNLTE

logger = logging.getLogger(__name__)

__all__ = ["HydrogenContinuumProperties", "Iteration"]


class Iteration(Input):
    """
    Attributes
    ----------
    iteration : int
        Current iteration number.
    """

    outputs = ("iteration",)


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


class HydrogenContinuumProperties(ProcessingPlasmaProperty):
    outputs = (
        "level_number_density",
        "stimulated_emission_factor",
        "plus other things",
    )

    def __init__(self, plasma_parent, photo_ion_cross_sections):
        super().__init__(plasma_parent)
        self._update_inputs()
        self.photoionization_data = photo_ion_cross_sections

    def calculate_hydrogen_ion_populations(
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

        solver = IonPopulationSolver(ion_rate_matrix_solver)

        fractional_ion_population, fractional_electron_density = solver.solve(
            dilute_planckian_radiation_field,
            electron_dist,
            elemental_number_density,
            lte_level_number_density,
            level_number_density,
            lte_ion_number_density,
            ion_number_density,
            partition_function,
            boltzmann_factor,
            charge_conservation=False,
        )

        ion_number_density = fractional_ion_population * number_density
        electron_densities = (
            fractional_electron_density * previous_electron_densities
        )

        return ion_number_density, electron_densities

    def calculate_lte_quantities(
        self,
        lte_level_number_density_solver,
        number_density,
        electron_number_density,
        levels,
        block_ids,
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

        lte_ion_number_density = LTEIonNumberDensity(
            self.plasma_parent, electron_densities=electron_number_density
        ).calculate(
            thermal_phi_lte,
            thermal_lte_partition_function,
            number_density,
            electron_number_density,
            block_ids,
        )

        lte_level_number_density = lte_level_number_density_solver.calculate(
            thermal_lte_level_boltzmann_factor,
            lte_ion_number_density,
            levels,
            thermal_lte_partition_function,
        )

        return lte_ion_number_density, lte_level_number_density

    def calculate(
        self,
        atomic_data,
        ionization_data,
        nlte_data,
        levels,
        t_electrons,
        dilute_planckian_radiation_field,
        general_level_boltzmann_factor,
        thermal_lte_level_boltzmann_factor,
        thermal_lte_partition_function,
        previous_electron_densities,
        number_density,
        g,
        thermal_g_electron,
        block_ids,
        t_rad,
        w,
        zeta_data,
        delta,
        g_electron,
        beta_rad,
        beta_electron,
        iteration,
    ):
        # set up initial values for major properties- level number density, ion
        # number density, electron number density, boltzmann factor, partition function, phi
        # and thermal LTE (electron temperature) versions of level number density,
        # ion number density, phi, boltzmann factor, partition function

        level_boltzmann_factor_solver = LevelBoltzmannFactorNLTE(
            self.plasma_parent
        )

        level_boltzmann_factor = (
            level_boltzmann_factor_solver._calculate_general(
                atomic_data,
                nlte_data,
                t_electrons,  # changes every loop
                dilute_planckian_radiation_field,
                None,
                general_level_boltzmann_factor,  # changes every loop
                previous_electron_densities,
                g,
            )
        )

        partition_function = level_boltzmann_factor.groupby(
            level=["atomic_number", "ion_number"]
        ).sum()

        phi = PhiSahaNebular.calculate(
            t_rad,
            w,
            zeta_data,
            t_electrons,  # changes every loop
            delta,  # changes every loop
            g_electron,
            beta_rad,
            partition_function,  # changes every loop
            ionization_data,
        )

        ion_number_density = IonNumberDensity(self.plasma_parent).calculate(
            phi,  # changes every loop
            partition_function,  # changes every loop
            number_density,
        )

        lte_level_number_density_solver = LTELevelNumberDensity(
            self.plasma_parent
        )

        level_number_density_solver = LevelNumberDensity(self.plasma_parent)

        level_number_density = level_number_density_solver.calculate(
            level_boltzmann_factor,  # changes every loop
            ion_number_density,  # changes every loop
            levels,
            partition_function,  # changes every loop
        )

        lte_level_number_density, lte_ion_number_density = (
            self.calculate_lte_quantities(
                lte_level_number_density_solver,
                number_density,
                previous_electron_densities,
                levels,
                block_ids,
                thermal_g_electron,  # changes every loop
                beta_electron,  # changes every loop
                ionization_data,
                thermal_lte_level_boltzmann_factor,  # changes every loop
                thermal_lte_partition_function,  # changes every loop
            )
        )

        level_number_density_converged = False

        electron_number_density = previous_electron_densities

        # Need to solve iteratively to get a consistent solution between
        # hydrogen and other elements. First iteration of the TARDIS simulation
        # uses the analytic solution for the hydrogen ion number density.
        # After that use estimators.

        while not level_number_density_converged:
            hydrogen_ion_number_density, electron_number_density = (
                self.calculate_hydrogen_ion_populations(
                    t_electrons,
                    previous_electron_densities,
                    dilute_planckian_radiation_field,
                    number_density,
                    lte_level_number_density,
                    level_number_density,
                    lte_ion_number_density,
                    ion_number_density,
                )
            )

            ion_number_density, electron_number_density = IonNumberDensity(
                self.plasma_parent, electron_densities=electron_number_density
            ).calculate(phi, partition_function, number_density)

            lte_level_number_density, lte_ion_number_density = (
                self.calculate_lte_quantities(lte_level_number_density_solver)
            )

            level_number_density = level_number_density_solver.calculate(
                level_boltzmann_factor,
                ion_number_density,
                levels,
                partition_function,
            )

            # update properties before next loop

            level_boltzmann_factor = (
                level_boltzmann_factor_solver._calculate_general(
                    atomic_data,
                    nlte_data,
                    t_electrons,
                    dilute_planckian_radiation_field,
                    None,
                    general_level_boltzmann_factor,
                    electron_number_density,
                    g,
                )
            )

            partition_function = level_boltzmann_factor.groupby(
                level=["atomic_number", "ion_number"]
            ).sum()

            level_number_density_converged = True

        stimulated_emission_factor = (
            self.plasma_parent.stimulated_emission_factor
        )

        return level_number_density, stimulated_emission_factor
