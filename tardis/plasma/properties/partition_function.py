import logging

import astropy.units as u
import numpy as np
import pandas as pd

from tardis.plasma.electron_energy_distribution import (
    ThermalElectronEnergyDistribution,
)
from tardis.plasma.equilibrium.level_populations import LevelPopulationSolver
from tardis.plasma.equilibrium.rate_matrix import RateMatrix
from tardis.plasma.equilibrium.rates import (
    RadiativeRatesSolver,
    ThermalCollisionalRateSolver,
)
from tardis.plasma.exceptions import PlasmaConfigError
from tardis.plasma.properties.base import ProcessingPlasmaProperty

logger = logging.getLogger(__name__)

__all__ = [
    "LevelBoltzmannFactorLTE",
    "LevelBoltzmannFactorDiluteLTE",
    "LevelBoltzmannFactorNoNLTE",
    "LevelBoltzmannFactorNLTE",
    "PartitionFunction",
    "ThermalLevelBoltzmannFactorLTE",
    "ThermalLTEPartitionFunction",
]


class LevelBoltzmannFactorLTE(ProcessingPlasmaProperty):
    """
    Attributes
    ----------
    general_level_boltzmann_factor : Pandas DataFrame, dtype float
        Level population proportionality values.
        Evaluated at the radiation temperature.
        Indexed by atomic number, ion number, level number.
        Columns corresponding to zones. Does not consider
        NLTE.
    """

    outputs = ("general_level_boltzmann_factor",)
    latex_name = ("bf_{i,j,k}",)
    latex_formula = (
        r"g_{i,j,k}e^{\dfrac{-\epsilon_{i,j,k}}{k_{\textrm{B}}T_{\textrm{rad}}}}",
    )

    @staticmethod
    def calculate(excitation_energy, g, beta_rad, levels):
        exponential = np.exp(np.outer(excitation_energy.values, -beta_rad))
        level_boltzmann_factor_array = g.values[np.newaxis].T * exponential
        level_boltzmann_factor = pd.DataFrame(
            level_boltzmann_factor_array,
            index=levels,
            columns=np.arange(len(beta_rad)),
            dtype=np.float64,
        )
        return level_boltzmann_factor


class ThermalLevelBoltzmannFactorLTE(LevelBoltzmannFactorLTE):
    """
    Attributes
    ----------
    thermal_lte_level_boltzmann_factor : Pandas DataFrame, dtype float
        Level population proportionality values for LTE.
        Evaluated at the temperature of the
        electron gas (thermal). Indexed
        by atomic number, ion number, level number.
        Columns corresponding to zones.
    """

    outputs = ("thermal_lte_level_boltzmann_factor",)
    latex_name = (r"bf_{i,j,k}^{\textrm{LTE}}(T_e)",)
    latex_formula = (
        r"g_{i,j,k}e^{\dfrac{-\epsilon_{i,j,k}}{k_{\textrm{B}}T_{\textrm{electron}}}}",
    )

    @staticmethod
    def calculate(excitation_energy, g, beta_electron, levels):
        return super(
            ThermalLevelBoltzmannFactorLTE, ThermalLevelBoltzmannFactorLTE
        ).calculate(excitation_energy, g, beta_electron, levels)


class LevelBoltzmannFactorDiluteLTE(ProcessingPlasmaProperty):
    """
    Attributes
    ----------
    general_level_boltzmann_factor : Pandas DataFrame, dtype float
        Level population proportionality values. Indexed
        by atomic number, ion number, level number.
        Columns corresponding to zones. Dilute radiation
        field means non-metastable level values are
        multiplied by an additional factor W. Does not
        consider NLTE.
    """

    outputs = ("general_level_boltzmann_factor",)
    latex_name = ("bf_{i,j,k}",)
    latex_formula = (
        r"Wg_{i,j,k}e^{\dfrac{-\epsilon_{i,j,k}}{k_{\textrm{B}}T_{\textrm{rad}}}}",
    )

    def calculate(
        self, levels, g, excitation_energy, beta_rad, w, metastability
    ):
        level_boltzmann_factor = LevelBoltzmannFactorLTE.calculate(
            excitation_energy, g, beta_rad, levels
        )
        level_boltzmann_factor[~metastability] *= w
        return level_boltzmann_factor


class LevelBoltzmannFactorNoNLTE(ProcessingPlasmaProperty):
    """
    Attributes
    ----------
    level_boltzmann_factor : Pandas DataFrame, dtype float
        Returns general_level_boltzmann_factor as this
        property is included if NLTE is not used.
    """

    outputs = ("level_boltzmann_factor",)

    @staticmethod
    def calculate(general_level_boltzmann_factor):
        return general_level_boltzmann_factor


class LevelBoltzmannFactorNLTE(ProcessingPlasmaProperty):
    """
    Attributes
    ----------
    level_boltzmann_factor : Pandas DataFrame, dtype float
        Returns general_level_boltzmann_factor but
        updated for those species treated in NLTE.
    """

    outputs = ("level_boltzmann_factor",)

    def calculate(self):
        raise AttributeError(
            "This attribute is not defined on the parent class."
            "Please use one of the subclasses."
        )

    @staticmethod
    def from_config(nlte_conf):
        if nlte_conf.classical_nebular and not nlte_conf.coronal_approximation:
            return LevelBoltzmannFactorNLTEClassic
        elif (
            nlte_conf.coronal_approximation and not nlte_conf.classical_nebular
        ):
            return LevelBoltzmannFactorNLTECoronal
        elif nlte_conf.coronal_approximation and nlte_conf.classical_nebular:
            raise PlasmaConfigError(
                "Both coronal approximation and classical nebular specified in the config."
            )
        else:
            return LevelBoltzmannFactorNLTEGeneral

    def __init__(self, plasma_parent):
        """
        Selects appropriate 'calculate' function based on NLTE config
        options selected.
        """
        super().__init__(plasma_parent)

        self._update_inputs()

    def _main_nlte_calculation(
        self,
        atomic_data,
        nlte_data,
        t_electrons,
        dilute_planckian_radiation_field,
        beta_sobolevs,
        general_level_boltzmann_factor,
        previous_electron_densities,
        g,
    ):
        """
        The core of the NLTE calculation, used with all possible config.
        options.
        """
        for species in nlte_data.nlte_species:
            logger.info(f"Calculating rates for species {species}")
            species_slice = (species[0], species[1], slice(None), slice(None))
            radiative_transitions = atomic_data.lines.loc[species_slice, :]
            radiative_rate_solver = RadiativeRatesSolver(radiative_transitions)

            col_strength_temperatures = atomic_data.collision_data_temperatures
            col_strengths = atomic_data.collision_data.loc[species_slice, :]

            collisional_rate_solver = ThermalCollisionalRateSolver(
                atomic_data.levels,
                radiative_transitions,
                col_strength_temperatures,
                col_strengths,
                "chianti",
            )
            rate_solvers = [
                (radiative_rate_solver, "radiative"),
                (collisional_rate_solver, "electron"),
            ]

            rate_matrix_solver = RateMatrix(rate_solvers, atomic_data.levels)

            electron_distribution = ThermalElectronEnergyDistribution(
                0,
                t_electrons * u.K,
                previous_electron_densities * u.g / u.cm**3,
            )

            rate_matrix = rate_matrix_solver.solve(
                dilute_planckian_radiation_field, electron_distribution
            )

            solver = LevelPopulationSolver(rate_matrix, atomic_data.levels)

            level_pops = solver.solve()

            pd.testing.assert_index_equal(
                general_level_boltzmann_factor.loc[species].index,
                level_pops.loc[species].index,
            )

            general_level_boltzmann_factor.loc[species] = (
                level_pops.loc[species]
                * g.loc[species][0]
                / level_pops.loc[species].iloc[0]
            ).values
        return general_level_boltzmann_factor

    def _calculate_classical_nebular(
        self,
        atomic_data,
        nlte_data,
        t_electrons,
        dilute_planckian_radiation_field,
        previous_beta_sobolev,
        general_level_boltzmann_factor,
        previous_electron_densities,
        g,
    ):
        """
        Performs NLTE calculations using the classical nebular treatment.
        All beta sobolev values taken as 1.
        """
        beta_sobolevs = 1.0

        general_level_boltzmann_factor = self._main_nlte_calculation(
            atomic_data,
            nlte_data,
            t_electrons,
            dilute_planckian_radiation_field,
            beta_sobolevs,
            general_level_boltzmann_factor,
            previous_electron_densities,
            g,
        )
        return general_level_boltzmann_factor

    def _calculate_coronal_approximation(
        self,
        atomic_data,
        nlte_data,
        t_electrons,
        dilute_planckian_radiation_field,
        previous_beta_sobolev,
        general_level_boltzmann_factor,
        previous_electron_densities,
        g,
    ):
        """
        Performs NLTE calculations using the coronal approximation.
        All beta sobolev values taken as 1 and j_blues taken as 0.
        """
        beta_sobolevs = 1.0
        j_blues = 0.0
        general_level_boltzmann_factor = self._main_nlte_calculation(
            atomic_data,
            nlte_data,
            t_electrons,
            dilute_planckian_radiation_field,
            beta_sobolevs,
            general_level_boltzmann_factor,
            previous_electron_densities,
            g,
        )
        return general_level_boltzmann_factor

    def _calculate_general(
        self,
        atomic_data,
        nlte_data,
        t_electrons,
        dilute_planckian_radiation_field,
        previous_beta_sobolev,
        general_level_boltzmann_factor,
        previous_electron_densities,
        g,
    ):
        """
        Full NLTE calculation without approximations.
        """
        if previous_beta_sobolev is None:
            beta_sobolevs = 1.0
        else:
            beta_sobolevs = previous_beta_sobolev

        general_level_boltzmann_factor = self._main_nlte_calculation(
            atomic_data,
            nlte_data,
            t_electrons,
            dilute_planckian_radiation_field,
            beta_sobolevs,
            general_level_boltzmann_factor,
            previous_electron_densities,
            g,
        )
        return general_level_boltzmann_factor


class LevelBoltzmannFactorNLTECoronal(LevelBoltzmannFactorNLTE):
    calculate = LevelBoltzmannFactorNLTE._calculate_coronal_approximation


class LevelBoltzmannFactorNLTEClassic(LevelBoltzmannFactorNLTE):
    calculate = LevelBoltzmannFactorNLTE._calculate_classical_nebular


class LevelBoltzmannFactorNLTEGeneral(LevelBoltzmannFactorNLTE):
    calculate = LevelBoltzmannFactorNLTE._calculate_general


class PartitionFunction(ProcessingPlasmaProperty):
    """
    Attributes
    ----------
    partition_function : Pandas DataFrame, dtype float
        Indexed by atomic number, ion number.
        Columns are zones.
    """

    outputs = ("partition_function",)
    latex_name = ("Z_{i,j}",)
    latex_formula = (r"\sum_{k}bf_{i,j,k}",)

    def calculate(self, level_boltzmann_factor):
        return level_boltzmann_factor.groupby(
            level=["atomic_number", "ion_number"]
        ).sum()


class ThermalLTEPartitionFunction(PartitionFunction):
    """
    Attributes
    ----------
    thermal_lte_partition_function : Pandas DataFrame, dtype float
        Indexed by atomic number, ion number.
        Columns are zones.
    """

    outputs = ("thermal_lte_partition_function",)
    latex_name = (r"Z_{i,j}(T_\mathrm{e}",)

    def calculate(self, thermal_lte_level_boltzmann_factor):
        return super().calculate(thermal_lte_level_boltzmann_factor)
