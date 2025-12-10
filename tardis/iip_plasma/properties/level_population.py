import logging

import numpy as np
import pandas as pd

from tardis.iip_plasma.continuum.util import *
from tardis.iip_plasma.properties.base import (
    PreviousIterationProperty,
    ProcessingPlasmaProperty,
)

logger = logging.getLogger(__name__)

__all__ = [
    "DepartureCoefficient",
    "LTELevelNumberDensity",
    "LevelNumberDensity",
    "PhiLucy",
    "PreviousDepartureCoefficient",
]


class LevelNumberDensity(ProcessingPlasmaProperty):
    """
    Attributes
    ----------
    level_number_density : Pandas DataFrame, dtype float
                           Index atom number, ion number, level number. Columns are zones.
    """

    outputs = ("level_number_density",)
    latex_name = ("N_{i,j,k}",)
    latex_formula = ("N_{i,j}\\dfrac{bf_{i,j,k}}{Z_{i,j}}",)

    def calculate(self):
        pass

    def __init__(self, plasma_parent, helium_treatment="dilute-lte"):
        """
        Calculates the level populations with the Boltzmann equation in LTE.
        """
        super(LevelNumberDensity, self).__init__(plasma_parent)
        if hasattr(self.plasma_parent, "plasma_properties_dict"):
            if "HeliumNLTE" in self.plasma_parent.plasma_properties_dict.keys():
                helium_treatment = "recomb-nlte"
        if helium_treatment == "recomb-nlte":
            self.calculate = self._calculate_helium_nlte
        elif helium_treatment == "dilute-lte":
            self.calculate = self._calculate_dilute_lte
        self._update_inputs()

        self.initialize_indices = True

    def _initialize_indices(self, levels, partition_function):
        indexer = pd.Series(
            np.arange(partition_function.shape[0]),
            index=partition_function.index,
        )
        self._ion2level_idx = indexer.loc[levels.droplevel(2)].values

    def _calculate_dilute_lte(
        self,
        level_boltzmann_factor,
        ion_number_density,
        levels,
        partition_function,
    ):
        """
        Reduces non-metastable level populations by a factor of W compared to LTE in the case of dilute-lte excitation.
        """
        if self.initialize_indices:
            self._initialize_indices(levels, partition_function)
            self.initialize_indices = False
        partition_function_broadcast = partition_function.values[
            self._ion2level_idx
        ]
        level_population_fraction = (
            level_boltzmann_factor.values / partition_function_broadcast
        )
        ion_number_density_broadcast = ion_number_density.values[
            self._ion2level_idx
        ]
        level_number_density = (
            level_population_fraction * ion_number_density_broadcast
        )
        # print "Calculating new lvl_number_densities"
        return pd.DataFrame(
            level_number_density, index=level_boltzmann_factor.index
        )

    def _calculate_helium_nlte(
        self,
        level_boltzmann_factor,
        ion_number_density,
        levels,
        partition_function,
        helium_population,
    ):
        """
        If one of the two helium NLTE methods is used, this updates the helium level populations to the appropriate
        values.
        """
        level_number_density = self._calculate_dilute_lte(
            level_boltzmann_factor,
            ion_number_density,
            levels,
            partition_function,
        )
        if helium_population is not None:
            level_number_density.loc[2].update(helium_population)
        return level_number_density


class LTELevelNumberDensity(LevelNumberDensity):
    outputs = ("lte_level_number_density",)
    latex_name = ("N_{i,j,k}^*",)

    def _calculate_dilute_lte(
        self,
        lte_level_boltzmann_factor_Te,
        lte_ion_number_density,
        levels,
        lte_partition_function_Te,
    ):
        return super(LTELevelNumberDensity, self)._calculate_dilute_lte(
            lte_level_boltzmann_factor_Te,
            lte_ion_number_density,
            levels,
            lte_partition_function_Te,
        )


class DepartureCoefficient(ProcessingPlasmaProperty):
    outputs = ("b",)

    def calculate(
        self,
        phi_lucy,
        electron_densities,
        ion_number_density,
        level_number_density,
    ):
        index = get_ion_multi_index(phi_lucy.index)
        lte_level_number_density = phi_lucy.multiply(
            (ion_number_density.loc[index] * electron_densities).values
        )
        return level_number_density / lte_level_number_density


class PreviousDepartureCoefficient(PreviousIterationProperty):
    outputs = ("previous_b",)

    def set_initial_value(self, kwargs):
        initial_value = np.ones(len(kwargs["abundance"].columns))
        initial_value = pd.DataFrame({(1, 0, 0): pd.Series(initial_value)}).T
        self._set_initial_value(initial_value)


class PhiLucy(ProcessingPlasmaProperty):
    outputs = ("phi_lucy",)
    latex_name = ("\\Phi_{i,\\kappa}",)

    def calculate(
        self, phi_Te, lte_level_boltzmann_factor_Te, lte_partition_function_Te
    ):
        boltzmann_factor = self._prepare_boltzmann_factor(
            lte_level_boltzmann_factor_Te
        )
        phi_saha_index = get_ion_multi_index(boltzmann_factor.index)
        partition_function_index = get_ion_multi_index(
            boltzmann_factor.index, next_higher=False
        )
        phi_saha = phi_Te.loc[phi_saha_index].values
        partition_function = lte_partition_function_Te.loc[
            partition_function_index
        ].values
        return boltzmann_factor / (phi_saha * partition_function)

    @staticmethod
    def _prepare_boltzmann_factor(boltzmann_factor):
        atomic_number = boltzmann_factor.index.get_level_values(0)
        ion_number = boltzmann_factor.index.get_level_values(1)
        selected_ions_mask = atomic_number != ion_number
        return boltzmann_factor[selected_ions_mask]
