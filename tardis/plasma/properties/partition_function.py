import logging

import numpy as np
from numpy.linalg.linalg import LinAlgError
import pandas as pd

from tardis.plasma.properties.base import ProcessingPlasmaProperty
from tardis.plasma.exceptions import PlasmaConfigError

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
        super(LevelBoltzmannFactorNLTE, self).__init__(plasma_parent)

        self._update_inputs()

    def _main_nlte_calculation(
        self,
        atomic_data,
        nlte_data,
        t_electrons,
        j_blues,
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
            number_of_levels = atomic_data.levels.energy.loc[species].count()
            lnl = nlte_data.lines_level_number_lower[species]
            lnu = nlte_data.lines_level_number_upper[species]
            (lines_index,) = nlte_data.lines_idx[species]

            try:
                j_blues_filtered = j_blues.iloc[lines_index]
            except AttributeError:
                j_blues_filtered = j_blues
            try:
                beta_sobolevs_filtered = beta_sobolevs.iloc[lines_index]
            except AttributeError:
                beta_sobolevs_filtered = beta_sobolevs
            A_uls = nlte_data.A_uls[species]
            B_uls = nlte_data.B_uls[species]
            B_lus = nlte_data.B_lus[species]
            r_lu_index = lnu * number_of_levels + lnl
            r_ul_index = lnl * number_of_levels + lnu
            r_ul_matrix = np.zeros(
                (number_of_levels, number_of_levels, len(t_electrons)),
                dtype=np.float64,
            )
            r_ul_matrix_reshaped = r_ul_matrix.reshape(
                (number_of_levels ** 2, len(t_electrons))
            )
            r_ul_matrix_reshaped[r_ul_index] = (
                A_uls[np.newaxis].T + B_uls[np.newaxis].T * j_blues_filtered
            )
            r_ul_matrix_reshaped[r_ul_index] *= beta_sobolevs_filtered
            r_lu_matrix = np.zeros_like(r_ul_matrix)
            r_lu_matrix_reshaped = r_lu_matrix.reshape(
                (number_of_levels ** 2, len(t_electrons))
            )
            r_lu_matrix_reshaped[r_lu_index] = (
                B_lus[np.newaxis].T * j_blues_filtered * beta_sobolevs_filtered
            )
            if atomic_data.collision_data is None:
                collision_matrix = np.zeros_like(r_ul_matrix)
            else:
                if previous_electron_densities is None:
                    collision_matrix = np.zeros_like(r_ul_matrix)
                else:
                    collision_matrix = (
                        nlte_data.get_collision_matrix(species, t_electrons)
                        * previous_electron_densities.values
                    )
            rates_matrix = r_lu_matrix + r_ul_matrix + collision_matrix
            for i in range(number_of_levels):
                rates_matrix[i, i] = -rates_matrix[:, i].sum(axis=0)
            rates_matrix[0, :, :] = 1.0
            x = np.zeros(rates_matrix.shape[0])
            x[0] = 1.0
            for i in range(len(t_electrons)):
                try:
                    level_boltzmann_factor = np.linalg.solve(
                        rates_matrix[:, :, i], x
                    )
                except LinAlgError as e:
                    if e.message == "Singular matrix":
                        raise ValueError(
                            "SingularMatrixError during solving of the "
                            "rate matrix. Does the atomic data contain "
                            "collision data?"
                        )
                    else:
                        raise e
                general_level_boltzmann_factor[i].loc[species] = (
                    level_boltzmann_factor
                    * g.loc[species][0]
                    / level_boltzmann_factor[0]
                )
        return general_level_boltzmann_factor

    def _calculate_classical_nebular(
        self,
        t_electrons,
        lines,
        atomic_data,
        nlte_data,
        general_level_boltzmann_factor,
        j_blues,
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
            j_blues,
            beta_sobolevs,
            general_level_boltzmann_factor,
            previous_electron_densities,
            g,
        )
        return general_level_boltzmann_factor

    def _calculate_coronal_approximation(
        self,
        t_electrons,
        lines,
        atomic_data,
        nlte_data,
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
            j_blues,
            beta_sobolevs,
            general_level_boltzmann_factor,
            previous_electron_densities,
            g,
        )
        return general_level_boltzmann_factor

    def _calculate_general(
        self,
        t_electrons,
        lines,
        atomic_data,
        nlte_data,
        general_level_boltzmann_factor,
        j_blues,
        previous_beta_sobolev,
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
            j_blues,
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
        return super(ThermalLTEPartitionFunction, self).calculate(
            thermal_lte_level_boltzmann_factor
        )
