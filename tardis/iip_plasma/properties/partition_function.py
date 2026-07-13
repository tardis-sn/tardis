import logging

import numpy as np
import pandas as pd
from numba import njit
from scipy.optimize import least_squares, root

from tardis.iip_plasma.exceptions import (
    PlasmaConfigError,
    PlasmaNLTEExcitationError,
)
from tardis.iip_plasma.properties.base import ProcessingPlasmaProperty
from tardis.opacities.tau_sobolev import calculate_beta_sobolev
from tardis.transport.montecarlo import njit_dict_no_fastmath

logger = logging.getLogger(__name__)

__all__ = [
    "LTEPartitionFunction",
    "LTEPartitionFunctionTe",
    "LevelBoltzmannFactorDiluteLTE",
    "LevelBoltzmannFactorLTE",
    "LevelBoltzmannFactorLTECont",
    "LevelBoltzmannFactorLTETe",
    "LevelBoltzmannFactorNLTE",
    "LevelBoltzmannFactorNoNLTE",
    "PartitionFunction",
]


@njit(njit_dict_no_fastmath)
def _calculate_beta_sobolevs_from_selected_lines(
    level_density_values: np.ndarray,
    line_indices: np.ndarray,
    lines_lower_level_index: np.ndarray,
    lines_upper_level_index: np.ndarray,
    g_lower: np.ndarray,
    g_upper: np.ndarray,
    meta_stable_upper: np.ndarray,
    nlte_lines_mask: np.ndarray,
    tau_coefficient: np.ndarray,
) -> np.ndarray:
    """Calculate Sobolev escape probabilities for selected lines.

    Parameters
    ----------
    level_density_values : numpy.ndarray
        Level number densities indexed by level position.
    line_indices : numpy.ndarray
        Line positions for which beta Sobolev values should be calculated.
    lines_lower_level_index : numpy.ndarray
        Level-density positions of the lower line levels.
    lines_upper_level_index : numpy.ndarray
        Level-density positions of the upper line levels.
    g_lower : numpy.ndarray
        Statistical weights of the lower line levels.
    g_upper : numpy.ndarray
        Statistical weights of the upper line levels.
    meta_stable_upper : numpy.ndarray
        Boolean mask marking lines with metastable upper levels.
    nlte_lines_mask : numpy.ndarray
        Boolean mask marking NLTE lines.
    tau_coefficient : numpy.ndarray
        Precomputed multiplicative coefficient for Sobolev optical depth.

    Returns
    -------
    numpy.ndarray
        Beta Sobolev values for ``line_indices``.
    """
    beta_sobolev = np.empty(line_indices.shape[0], dtype=np.float64)
    for i in range(line_indices.shape[0]):
        line_index = line_indices[i]
        n_lower = level_density_values[lines_lower_level_index[line_index]]
        n_upper = level_density_values[lines_upper_level_index[line_index]]

        stimulated_emission_factor = 1 - (
            g_lower[line_index] * n_upper / (g_upper[line_index] * n_lower)
        )
        if (
            n_lower == 0.0
            or np.isneginf(stimulated_emission_factor)
            or (
                meta_stable_upper[line_index] and stimulated_emission_factor < 0
            )
            or (nlte_lines_mask[line_index] and stimulated_emission_factor < 0)
        ):
            stimulated_emission_factor = 0.0

        tau_sobolev = (
            tau_coefficient[line_index] * n_lower * stimulated_emission_factor
        )
        if tau_sobolev > 1e3:
            beta_sobolev[i] = tau_sobolev**-1
        elif tau_sobolev < 1e-4:
            beta_sobolev[i] = 1 - 0.5 * tau_sobolev
        else:
            beta_sobolev[i] = (1 - np.exp(-tau_sobolev)) / tau_sobolev
    return beta_sobolev[:, np.newaxis]


class LevelBoltzmannFactorLTE(ProcessingPlasmaProperty):
    """
    Attributes
    ----------
    general_level_boltzmann_factor : Pandas DataFrame, dtype float
                             Level population proportionality values. Indexed
                             by atomic number, ion number, level number.
                             Columns corresponding to zones. Does not consider
                             NLTE.
    """

    outputs = ("general_level_boltzmann_factor",)
    latex_name = ("bf_{i,j,k}",)
    latex_formula = (
        "g_{i,j,k}e^{\\dfrac{-\\epsilon_{i,j,k}}{k_{\
        \\textrm{B}}T_{\\textrm{rad}}}}",
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


class LevelBoltzmannFactorLTETe(LevelBoltzmannFactorLTE):
    """
    Attributes
    ----------
    level_boltzmann_factor_LTE_Te : Pandas DataFrame, dtype float
                             Level population proportionality values for LTE.
                             Evaluated at the kinetic temperature T_e. Indexed
                             by atomic number, ion number, level number.
                             Columns corresponding to zones.
    """

    outputs = ("lte_level_boltzmann_factor_Te",)
    latex_name = ("bf_{i,j,k}^{\\textrm{LTE}}(T_e)",)
    latex_formula = (
        "g_{i,j,k}e^{\\dfrac{-\\epsilon_{i,j,k}}{k_{\
        \\textrm{B}}T_{\\textrm{electron}}}}",
    )

    @staticmethod
    def calculate(excitation_energy, g, beta_electron, levels):
        return super(
            LevelBoltzmannFactorLTETe, LevelBoltzmannFactorLTETe
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
        "Wg_{i,j,k}e^{\\dfrac{-\\epsilon_{i,j,k}}{k_{\
        \\textrm{B}}T_{\\textrm{rad}}}}",
    )

    def calculate(
        self, levels, g, excitation_energy, beta_rad, w, metastability
    ):
        # print "Calculating LevelBoltzmannFactorDiluteLTE\n"
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

    outputs = ("level_boltzmann_factor", "ion_ratio")

    def calculate(self):
        pass

    def __init__(
        self,
        plasma_parent,
        classical_nebular=False,
        coronal_approximation=False,
    ):
        """
        Selects appropriate 'calculate' function based on NLTE config
        options selected.
        """
        super(LevelBoltzmannFactorNLTE, self).__init__(plasma_parent)
        if classical_nebular == True and coronal_approximation == False:
            self.calculate = self._calculate_classical_nebular
        elif coronal_approximation == True and classical_nebular == False:
            self.calculate = self._calculate_coronal_approximation
        elif coronal_approximation == True and classical_nebular == True:
            raise PlasmaConfigError(
                "Both coronal approximation and classical"
                "nebular specified in the config."
            )
        elif self.plasma_parent.continuum_treatment:
            self.calculate = self._calculate_with_continuum
        else:
            self.calculate = self._calculate_general
        self._beta_sobolev_inputs = None
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
    ):
        """
        The core of the NLTE calculation, used with all possible config.
        options.
        """
        for species in self.plasma_parent.nlte_species:
            # logger.info('Calculating rates for species %s', species)
            rates_matrix, x, lines_idx, r_lu_matrix, r_ul_matrix = (
                self._setup_bb_rates(
                    species,
                    atomic_data,
                    nlte_data,
                    t_electrons,
                    j_blues,
                    beta_sobolevs,
                    previous_electron_densities,
                )
            )

            for i in range(len(t_electrons)):
                level_boltzmann_factor = np.linalg.solve(
                    rates_matrix[:, :, i], x
                )
                general_level_boltzmann_factor[i].loc[species] = (
                    level_boltzmann_factor
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
        lte_j_blues,
        previous_electron_densities,
    ):
        """
        Performs NLTE calculations using the classical nebular treatment.
        All beta sobolev values taken as 1.
        """
        beta_sobolevs = np.ones((len(lines), len(t_electrons)))
        if len(j_blues) == 0:
            j_blues = lte_j_blues
        else:
            j_blues = pd.DataFrame(
                j_blues, index=lines.index, columns=range(len(t_electrons))
            )
        general_level_boltzmann_factor = self._main_nlte_calculation(
            atomic_data,
            nlte_data,
            t_electrons,
            j_blues,
            beta_sobolevs,
            general_level_boltzmann_factor,
            previous_electron_densities,
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
    ):
        """
        Performs NLTE calculations using the coronal approximation.
        All beta sobolev values taken as 1 and j_blues taken as 0.
        """
        beta_sobolevs = np.ones((len(lines), len(t_electrons)))
        j_blues = np.zeros((len(lines), len(t_electrons)))
        general_level_boltzmann_factor = self._main_nlte_calculation(
            atomic_data,
            nlte_data,
            t_electrons,
            j_blues,
            beta_sobolevs,
            general_level_boltzmann_factor,
            previous_electron_densities,
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
        lte_j_blues,
        previous_electron_densities,
    ):
        """
        Full NLTE calculation without approximations.
        """
        if previous_beta_sobolev is None:
            beta_sobolevs = np.ones((len(lines), len(t_electrons)))
        else:
            beta_sobolevs = previous_beta_sobolev
        if len(j_blues) == 0:
            j_blues = lte_j_blues
        else:
            j_blues = pd.DataFrame(
                j_blues, index=lines.index, columns=range(len(t_electrons))
            )
        general_level_boltzmann_factor = self._main_nlte_calculation(
            atomic_data,
            nlte_data,
            t_electrons,
            j_blues,
            beta_sobolevs,
            general_level_boltzmann_factor,
            previous_electron_densities,
        )
        return general_level_boltzmann_factor

    def _calculate_with_continuum(
        self,
        t_electrons,
        lines,
        atomic_data,
        nlte_data,
        general_level_boltzmann_factor,
        j_blues,
        previous_beta_sobolev,
        lte_j_blues,
        previous_electron_densities,
        gamma,
        previous_ion_number_density,
        alpha_stim,
        alpha_sp,
        coll_exc_coeff,
        coll_deexc_coeff,
        coll_ion_coeff,
        coll_recomb_coeff,
    ):
        """
        Full NLTE calculation without approximations.
        """
        if previous_beta_sobolev is None:
            beta_sobolevs = np.ones((len(lines), len(t_electrons)))
        else:
            beta_sobolevs = previous_beta_sobolev
        if len(j_blues) == 0:
            j_blues = lte_j_blues
        else:
            j_blues = pd.DataFrame(
                j_blues, index=lines.index, columns=range(len(t_electrons))
            )
        general_level_boltzmann_factor, ion_ratio = (
            self._main_nlte_cont_calculation(
                atomic_data,
                nlte_data,
                t_electrons,
                j_blues,
                beta_sobolevs,
                general_level_boltzmann_factor,
                previous_electron_densities,
                gamma,
                previous_ion_number_density,
                alpha_stim,
                alpha_sp,
                coll_exc_coeff,
                coll_deexc_coeff,
                coll_ion_coeff,
                coll_recomb_coeff,
            )
        )
        return general_level_boltzmann_factor, ion_ratio

    def _main_nlte_cont_calculation(
        self,
        atomic_data,
        nlte_data,
        t_electrons,
        j_blues,
        beta_sobolevs,
        general_level_boltzmann_factor,
        previous_electron_densities,
        gamma,
        previous_ion_number_density,
        alpha_stim,
        alpha_sp,
        coll_exc_coeff,
        coll_deexc_coeff,
        coll_ion_coeff,
        coll_recomb_coeff,
    ):
        """
        The core of the NLTE calculation, used with all possible config.
        options.
        """
        for species in self.plasma_parent.nlte_species:
            # logger.info('Calculating rates for species %s', species)
            (
                rates_matrix,
                x,
                lines_idx,
                r_lu_matrix,
                r_ul_matrix,
                r_lu_index,
                r_ul_index,
            ) = self._setup_bb_rates(
                species,
                atomic_data,
                nlte_data,
                t_electrons,
                j_blues,
                beta_sobolevs,
                previous_electron_densities,
                coll_exc_coeff,
                coll_deexc_coeff,
            )

            no_of_levels = atomic_data.levels.energy.loc[species].count()
            # general_level_boltzmann_factor_old = general_level_boltzmann_factor.copy(deep=True)
            ion_ratio = np.zeros_like(t_electrons)

            if previous_ion_number_density is not None and gamma is not None:
                gamma1 = gamma.loc[species]
                alpha_sp1 = alpha_sp.loc[species]
                alpha_stim1 = alpha_stim.loc[species]
                coll_ion_coeff1 = coll_ion_coeff.loc[species]
                coll_recomb_coeff1 = coll_recomb_coeff.loc[species]
                coll_deexc_coeff_species = coll_deexc_coeff.loc[species]
                coll_exc_coeff_species = coll_exc_coeff.loc[species]
                collision_rate_matrices = self._setup_collision_rate_matrices(
                    coll_exc_coeff_species,
                    coll_deexc_coeff_species,
                    no_of_levels,
                    previous_electron_densities,
                )

            for i in range(len(t_electrons)):
                if (
                    previous_ion_number_density is not None
                    and gamma is not None
                ):
                    n_e = previous_electron_densities[i]

                    coll_ion = coll_ion_coeff1[i].values * n_e
                    rates_matrix[:, :, i] -= np.diag(
                        gamma1[i].values + coll_ion
                    )

                    gamma_vec = gamma1[i].values + coll_ion
                    next_ion_index = (species[0], species[1] + 1)
                    next_ion_density = previous_ion_number_density.loc[
                        next_ion_index, i
                    ]
                    next_ion_density /= previous_ion_number_density.loc[
                        species, i
                    ]

                    alpha = alpha_stim1[i].values + alpha_sp1[i].values
                    coll_recomb = coll_recomb_coeff1[i].values
                    # TODO: setup x outside of loop
                    x = (alpha + coll_recomb * n_e) * n_e  # * next_ion_density

                    collision_rates = collision_rate_matrices[:, :, i]
                    rates_matrix[:, :, i] += collision_rates

                    # x[0] = 0.
                    rates_matrix[0, :, i] = 1.0

                    remaining_rates = (
                        -np.diag(gamma1[i].values + coll_ion) + collision_rates
                    )
                    # TODO: passing self is horrible
                    # args = {'self': self, 'x_vec': x, 'remaining_rates': remaining_rates}
                    phis = (
                        self.plasma_parent.phi.loc[(species[0],), i].values
                        / n_e
                    )
                    level_density_base_values = (
                        self.plasma_parent.level_number_density[0].to_numpy(
                            copy=True
                        )
                    )
                    species_level_positions = self._get_species_level_positions(
                        species
                    )

                    args = (
                        self,
                        x,
                        remaining_rates,
                        lines_idx,
                        r_lu_matrix[:, :, i],
                        r_ul_matrix[:, :, i],
                        r_lu_index,
                        r_ul_index,
                        self.plasma_parent.number_density.loc[species[0], i],
                        gamma_vec,
                        species,
                        phis,
                        level_density_base_values,
                        species_level_positions,
                    )

                    initial = (
                        self.plasma_parent.level_number_density[i]
                        .loc[species]
                        .values
                    ).copy()
                    initial /= initial.sum()
                    initial = np.hstack([initial, np.array(next_ion_density)])

                    res = root(self._rate_equations, x0=initial, args=args)

                    if (res.x[:-1] < 0).sum():
                        logger.warning(
                            "Resorting to Lsq Solver for excitation in shell %d",
                            i,
                        )
                        ubound = [1.0] * (len(initial) - 1) + [np.inf]
                        res = least_squares(
                            self._rate_equations,
                            x0=initial,
                            args=args,
                            bounds=(0.0, ubound),
                        )

                        if (res.x[:-1] < 0).sum():
                            raise PlasmaNLTEExcitationError(
                                "level_boltzmann_factor is negative - aborting"
                            )

                    level_boltzmann_factor = res.x[:-1]
                    ion_ratio[i] = res.x[-1]

                    general_level_boltzmann_factor.loc[species, i] = (
                        level_boltzmann_factor
                    )

        return general_level_boltzmann_factor, ion_ratio

    @staticmethod
    def _rate_equations(
        x,
        self,
        x_vec,
        remaining_rates,
        lines_idx,
        r_lu_matrix,
        r_ul_matrix,
        r_lu_index,
        r_ul_index,
        number_density,
        gamma_vec,
        species,
        phis,
        level_density_base_values,
        species_level_positions,
    ):
        ion_numbers = self._get_ion_numbers(
            species, phis, number_density, phi_nlte=x[-1]
        )
        level_density_nlte_species = ion_numbers[species[1]] * x[:-1]
        level_density_values = level_density_base_values.copy()
        level_density_values[species_level_positions] = (
            level_density_nlte_species
        )

        beta_sobolev = self._calculate_beta_sobolevs_from_values(
            level_density_values, lines_idx
        )
        radiative_rates = self._setup_radiative_rates(
            beta_sobolev,
            r_ul_matrix.copy(),
            r_lu_matrix.copy(),
            r_lu_index,
            r_ul_index,
        )

        # print "Set radiative rates to zero"
        # radiative_rates=0
        rates_matrix = radiative_rates + remaining_rates

        rates_matrix[0, :] = 1.0
        x_tot = -x_vec.sum()
        x_vec[0] = 0.0

        no_of_levels = rates_matrix.shape[0]
        augmented_rates_matrix = np.empty((no_of_levels + 1, no_of_levels + 1))
        augmented_rates_matrix[:no_of_levels, :no_of_levels] = rates_matrix
        augmented_rates_matrix[:no_of_levels, no_of_levels] = x_vec
        # gamma_vec[0] = 0.0
        augmented_rates_matrix[no_of_levels, :no_of_levels] = gamma_vec
        augmented_rates_matrix[no_of_levels, no_of_levels] = x_tot
        num_vec = np.zeros(no_of_levels + 1)
        num_vec[0] = 1.0

        func = np.dot(augmented_rates_matrix, x) - num_vec
        # jac = np.ones((20,20))

        # level_boltzmann_factor = np.linalg.solve(rates_matrix, x_vec)
        return func

    @staticmethod
    def _get_ion_numbers(species, phis, number_density, phi_nlte):
        number_ratio_matrix = np.diag(np.ones(species[0]), k=1)
        number_ratio_matrix[-1] = 1.0
        diag_indices = np.diag_indices(len(number_ratio_matrix) - 1)
        number_ratio_matrix[diag_indices] = -phis
        number_ratio_matrix[species[1], species[1]] = -phi_nlte
        number_conservation_vec = np.zeros(len(number_ratio_matrix))
        number_conservation_vec[-1] = number_density
        ion_numbers = np.linalg.solve(
            number_ratio_matrix, number_conservation_vec
        )
        return ion_numbers

    def _setup_radiative_rates(
        self,
        beta_sobolev,
        r_ul_matrix,
        r_lu_matrix,
        r_lu_index,
        r_ul_index,
    ):
        # number_of_levels = 20
        # r_ul_matrix_reshaped = r_ul_matrix.reshape(number_of_levels**2)
        # r_lu_matrix_reshaped = r_lu_matrix.reshape(number_of_levels**2)

        r_ul_matrix_reshaped = r_ul_matrix.reshape(-1)
        r_lu_matrix_reshaped = r_lu_matrix.reshape(-1)

        # TODO: ? is this the right thing to do
        selected_beta_sobolev = beta_sobolev.ravel()
        r_ul_matrix_reshaped[r_ul_index] *= selected_beta_sobolev
        r_lu_matrix_reshaped[r_lu_index] *= selected_beta_sobolev

        rates_matrix = r_lu_matrix + r_ul_matrix

        # TODO: Revert
        # rates_matrix[0,:] = 0.0
        # rates_matrix[:,0] = 0.0

        # rates_matrix[0,:3] = 0.0
        # rates_matrix[:3,0] = 0.0

        # rates_matrix[0, 5:] = 0.0
        # rates_matrix[5:,0] = 0.0

        # rates_matrix[0,:5] = 0.0
        # rates_matrix[:5,0] = 0.0

        for i in range(r_ul_matrix.shape[0]):
            rates_matrix[i, i] = -rates_matrix[:, i].sum(axis=0)

        return rates_matrix

    def _caculate_beta_sobolevs(self, level_number_density):
        pl = self.plasma_parent

        tau_sobolev_calc = self.plasma_parent.plasma_properties_dict[
            "TauSobolev"
        ].calculate
        # beta_sobolev_calc = self.plasma_parent.plasma_properties_dict['BetaSobolev'].calculate
        stim_emission_factor_calc = self.plasma_parent.plasma_properties_dict[
            "StimulatedEmissionFactor"
        ].calculate

        stim_emission_factor = stim_emission_factor_calc(
            pl.g,
            level_number_density,
            pl.lines_lower_level_index,
            pl.lines_upper_level_index,
            pl.metastability,
            pl.lines,
        )

        tau_sobolev = tau_sobolev_calc(
            pl.lines,
            level_number_density,
            pl.lines_lower_level_index,
            pl.time_explosion,
            stim_emission_factor,
            pl.j_blues,
            pl.f_lu,
            pl.wavelength_cm,
        )

        # import ipdb; ipdb.set_trace()
        beta_sobolev = np.zeros_like(tau_sobolev.values)

        beta_sobolev = calculate_beta_sobolev(tau_sobolev)

        return beta_sobolev.values

    def _get_species_level_positions(
        self, species: tuple[int, int]
    ) -> np.ndarray:
        """Return cached level-number-density positions for a species.

        Parameters
        ----------
        species : tuple[int, int]
            Atomic number and ion number identifying the species.

        Returns
        -------
        numpy.ndarray
            Integer positions in the level number density index for the
            requested species.
        """
        if not hasattr(self, "_species_level_positions"):
            self._species_level_positions = {}
        if species not in self._species_level_positions:
            level_index = self.plasma_parent.level_number_density.index
            self._species_level_positions[species] = np.flatnonzero(
                (level_index.get_level_values("atomic_number") == species[0])
                & (level_index.get_level_values("ion_number") == species[1])
            )
        return self._species_level_positions[species]

    def _get_cached_beta_sobolev_inputs(self) -> dict[str, np.ndarray | float]:
        """Return cached inputs shared by beta Sobolev array calculations.

        Returns
        -------
        dict[str, numpy.ndarray or float]
            Arrays and scalar coefficients derived from plasma state that do
            not change during an individual beta Sobolev solve.
        """
        if self._beta_sobolev_inputs is None:
            pl = self.plasma_parent
            g_values = pl.g.to_numpy()
            lines_lower_level_index = pl.lines_lower_level_index
            lines_upper_level_index = pl.lines_upper_level_index
            atomic_numbers = pl.lines.index.get_level_values(
                "atomic_number"
            ).to_numpy()
            ion_numbers = pl.lines.index.get_level_values(
                "ion_number"
            ).to_numpy()
            nlte_lines_mask = np.zeros(len(pl.lines), dtype=bool)
            for atomic_number, ion_number in self.plasma_parent.nlte_species:
                nlte_lines_mask |= (atomic_numbers == atomic_number) & (
                    ion_numbers == ion_number
                )
            tau_sobolev_property = pl.plasma_properties_dict["TauSobolev"]
            time_explosion = pl.time_explosion
            if hasattr(time_explosion, "to_value"):
                time_explosion = time_explosion.to_value("s")
            self._beta_sobolev_inputs = {
                "g_lower": g_values[lines_lower_level_index],
                "g_upper": g_values[lines_upper_level_index],
                "lines_lower_level_index": lines_lower_level_index,
                "lines_upper_level_index": lines_upper_level_index,
                "meta_stable_upper": pl.metastability.to_numpy()[
                    lines_upper_level_index
                ],
                "nlte_lines_mask": nlte_lines_mask,
                "tau_coefficient": (
                    tau_sobolev_property.sobolev_coefficient
                    * pl.f_lu.to_numpy()
                    * pl.wavelength_cm.to_numpy()
                    * time_explosion
                ),
            }
        return self._beta_sobolev_inputs

    def _calculate_beta_sobolevs_from_values(
        self,
        level_density_values: np.ndarray,
        line_indices: np.ndarray | slice | None = None,
    ) -> np.ndarray:
        """Calculate beta Sobolev values from level-density arrays.

        Parameters
        ----------
        level_density_values : numpy.ndarray
            Level number densities aligned with the plasma level index.
        line_indices : numpy.ndarray or slice or None, optional
            Line positions to calculate. If None, all lines are calculated.

        Returns
        -------
        numpy.ndarray
            Beta Sobolev values for the requested line positions.
        """
        cached_inputs = self._get_cached_beta_sobolev_inputs()
        if line_indices is None:
            line_indices = np.arange(len(cached_inputs["tau_coefficient"]))
        elif isinstance(line_indices, tuple):
            line_indices = line_indices[0]
        beta_sobolev = _calculate_beta_sobolevs_from_selected_lines(
            level_density_values,
            line_indices,
            cached_inputs["lines_lower_level_index"],
            cached_inputs["lines_upper_level_index"],
            cached_inputs["g_lower"],
            cached_inputs["g_upper"],
            cached_inputs["meta_stable_upper"],
            cached_inputs["nlte_lines_mask"],
            cached_inputs["tau_coefficient"],
        )
        if np.any(np.isnan(beta_sobolev)) or np.any(
            np.isinf(np.abs(beta_sobolev))
        ):
            raise ValueError(
                "Some beta_sobolevs are nan, inf, -inf in beta_sobolevs."
                " Something went wrong!"
            )
        return beta_sobolev

    @staticmethod
    def _setup_bb_rates(
        species,
        atomic_data,
        nlte_data,
        t_electrons,
        j_blues,
        beta_sobolevs,
        previous_electron_densities,
        coll_exc_coeff=None,
        coll_deexc_coeff=None,
    ):
        j_blues = j_blues.values
        number_of_levels = atomic_data.levels.energy.loc[species].count()
        lnl = nlte_data.lines_level_number_lower[species]
        lnu = nlte_data.lines_level_number_upper[species]
        lines_index = nlte_data.lines_idx[species]
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
            (
                number_of_levels**2,
                len(t_electrons),
            )
        )
        r_ul_matrix_reshaped[r_ul_index] = (
            A_uls[np.newaxis].T + B_uls[np.newaxis].T * j_blues[lines_index]
        )
        # TODO: Revert
        # r_ul_matrix_reshaped[r_ul_index] *= beta_sobolevs[lines_index]
        r_lu_matrix = np.zeros_like(r_ul_matrix)
        r_lu_matrix_reshaped = r_lu_matrix.reshape(
            (
                number_of_levels**2,
                len(t_electrons),
            )
        )
        # r_lu_matrix_reshaped[r_lu_index] = B_lus[np.newaxis].T * \
        #        j_blues[lines_index] * beta_sobolevs[lines_index]
        r_lu_matrix_reshaped[r_lu_index] = (
            B_lus[np.newaxis].T * j_blues[lines_index]
        )
        if atomic_data.has_collision_data:
            if previous_electron_densities is None:
                collision_matrix = r_ul_matrix.copy()
                collision_matrix.fill(0.0)
            else:
                collision_matrix = (
                    nlte_data.get_collision_matrix(species, t_electrons)
                    * previous_electron_densities.values
                )
        elif coll_exc_coeff is not None:
            collision_matrix = r_ul_matrix.copy()
            collision_matrix.fill(0.0)
        else:
            collision_matrix = r_ul_matrix.copy()
            collision_matrix.fill(0.0)
        rates_matrix = r_lu_matrix + r_ul_matrix + collision_matrix
        for i in range(number_of_levels):
            rates_matrix[i, i] = -rates_matrix[:, i].sum(axis=0)
        rates_matrix[0, :, :] = 1.0
        x = np.zeros(rates_matrix.shape[0])
        x[0] = 1.0
        # import ipdb; ipdb.set_trace()
        return (
            rates_matrix,
            x,
            lines_index,
            r_lu_matrix,
            r_ul_matrix,
            r_lu_index,
            r_ul_index,
        )

    def _setup_collision_rate_matrix(
        self, coll_exc_coeff, coll_deexc_coeff, no_of_levels, shell, n_e
    ):
        collision_rate_matrices = self._setup_collision_rate_matrices(
            coll_exc_coeff[[shell]],
            coll_deexc_coeff[[shell]],
            no_of_levels,
            np.array([n_e]),
        )
        return collision_rate_matrices[:, :, 0]

    def _setup_collision_rate_matrices(
        self,
        coll_exc_coeff,
        coll_deexc_coeff,
        no_of_levels,
        electron_densities,
    ):
        """Build collisional rate matrices for all selected shells.

        Parameters
        ----------
        coll_exc_coeff : pandas.DataFrame
            Collisional excitation coefficients indexed by lower and upper
            level number.
        coll_deexc_coeff : pandas.DataFrame
            Collisional de-excitation coefficients indexed by upper and lower
            level number.
        no_of_levels : int
            Number of levels in the species rate-equation system.
        electron_densities : numpy.ndarray
            Electron densities for the shells being solved.

        Returns
        -------
        numpy.ndarray
            Collisional rate matrices with shape
            ``(no_of_levels, no_of_levels, len(electron_densities))``.
        """
        electron_densities = np.asarray(electron_densities)
        coll_deexc_coeff = coll_deexc_coeff.copy(deep=False)
        coll_deexc_coeff.index = coll_deexc_coeff.index.swaplevel(0, 1)

        index = list(self._get_rate_index(no_of_levels))
        coll_excitation_rates = (
            coll_exc_coeff.reindex(index)
            .fillna(0)
            .to_numpy()
            .reshape((no_of_levels, no_of_levels, len(electron_densities)))
        )
        coll_deexcitation_rates = (
            coll_deexc_coeff.reindex(index)
            .fillna(0)
            .to_numpy()
            .reshape((no_of_levels, no_of_levels, len(electron_densities)))
        )

        diagonal = np.zeros((no_of_levels, len(electron_densities)))
        coll_exc_lvl_sum = coll_exc_coeff.groupby(level=0).sum()
        coll_deexc_lvl_sum = coll_deexc_coeff.groupby(level=0).sum()
        diagonal[coll_exc_lvl_sum.index.astype(int)] += coll_exc_lvl_sum.values
        diagonal[coll_deexc_lvl_sum.index.astype(int)] += (
            coll_deexc_lvl_sum.values
        )

        collision_rate_matrices = (
            coll_deexcitation_rates + coll_excitation_rates
        ) * electron_densities
        for shell in range(len(electron_densities)):
            collision_rate_matrices[:, :, shell] -= np.diag(
                diagonal[:, shell] * electron_densities[shell]
            )

        # x = np.zeros(collision_rate_matrix.shape[0])
        # x[0] = 1
        # collision_rate_matrix[0,:] = 1
        # level_boltzmann_factor = np.linalg.solve(collision_rate_matrix, x)
        return collision_rate_matrices

    @staticmethod
    def _get_rate_index(no_of_levels):
        x_ind, y_ind = np.indices((no_of_levels, no_of_levels))
        index = zip(y_ind.flatten(), x_ind.flatten())
        return index


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
    latex_formula = ("\\sum_{k}bf_{i,j,k}",)

    @staticmethod
    def calculate(level_boltzmann_factor):
        return level_boltzmann_factor.groupby(
            level=["atomic_number", "ion_number"]
        ).sum()


class LevelBoltzmannFactorLTECont(LevelBoltzmannFactorLTE):
    outputs = ("lte_level_boltzmann_factor",)


class LTEPartitionFunction(PartitionFunction):
    outputs = ("lte_partition_function",)

    def calculate(self, lte_level_boltzmann_factor):
        return super(LTEPartitionFunction, self).calculate(
            lte_level_boltzmann_factor
        )


class LTEPartitionFunctionTe(PartitionFunction):
    outputs = ("lte_partition_function_Te",)

    def calculate(self, lte_level_boltzmann_factor_Te):
        return super(LTEPartitionFunctionTe, self).calculate(
            lte_level_boltzmann_factor_Te
        )
