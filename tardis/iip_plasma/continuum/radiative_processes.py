import logging

import astropy.units as units
import numpy as np
import pandas as pd
from astropy import constants as const
from scipy.integrate import cumulative_trapezoid, trapezoid

from tardis.iip_plasma.continuum.base import (
    BoundFreeEnergyMixIn,
    PhysicalContinuumProcess,
)
from tardis.iip_plasma.continuum.constants import continuum_constants as cconst
from tardis.plasma.exceptions import PlasmaException
from tardis.util.base import intensity_black_body

logger = logging.getLogger(__name__)


class RadiativeIonization(PhysicalContinuumProcess, BoundFreeEnergyMixIn):
    """
    Represents the process of radiative ionization.

    Attributes
    ----------
    input: `tardis.iip_plasma.continuum.input_data.ContinuumInputData`-object
        The common input data object.
    rate_coefficient: pd.DataFrame
        Multiplying the rate coefficient with the number densities of the interacting particles gives the rate
        per unit volume of the transition.

    Class Attributes
    ----------------
    name: str
        The name used in setattr(object, name, value).
    cooling: bool
        True if the physical process contributes to the cooling of the plasma. Enables calculation of cooling_rate.
    macro_atom_transitions: str
        The type of transitions in the macro atom.
    """

    name = "radiative_ionization"
    cooling = False
    macro_atom_transitions = "continuum"

    def __init__(self, input_data):
        super(RadiativeIonization, self).__init__(input_data)

    def _calculate_rate_coefficient(self, **kwargs):
        rate_coefficient_dilute_bb = (
            self._calculate_rate_coefficient_dilute_blackbody()
        )
        if not self.has_estimators:
            logger.info(
                "Calculating photoionization rate from dilute-blackbody radiation field model"
            )
            rate_coefficient = rate_coefficient_dilute_bb
        else:
            logger.info("Calculating photoionization rate from MC estimators")
            rate_coefficient = self._calculate_rate_coefficient_from_estimator()

            no_of_bad_elements = self._check_for_low_statistics(
                self.estimators["photo_ion_statistics"]
            )
            if self.replace_values_with_low_statistics and (
                no_of_bad_elements != 0
            ):
                logger.info(
                    "Replacing {} photoionization rates with values based on the "
                    "radiation field model".format(no_of_bad_elements)
                )

                rate_coefficient = self._calculate_rate_coefficient_combination(
                    rate_coefficient, rate_coefficient_dilute_bb
                )
        return rate_coefficient

    def _calculate_rate_coefficient_dilute_blackbody(self):
        """
        Calculates the rate coefficient for photoionization based on a dilute-blackbody radiation-field model.
        Stimulated recombinations are treated as negative photoionizations.

        Returns
        -------
        corrected_photoion_coeff: pd.DataFrame

        Notes
        -----
        .. math::

            \tilde \gamma_i = 4 \pi \int_{\nu_i}^{\infty} \! \frac{\tilde a_{i \kappa}}{h \nu} J_{\nu}
             \, \mathrm{d}\nu \\
            \tilde a_{i\kappa}(\nu)= a_{i\kappa}(\nu)\left(1-\frac{n_{\kappa} n_i^*}{n_i n_{\kappa}^*}
            e^{-h\nu/kT} \right)

        """
        j_nus = self._calculate_j_nus()
        stimulated_emission_correction = (
            self._calculate_stimulated_emission_correction()
        )
        corrected_photoion_coeff = j_nus.multiply(
            4.0
            * np.pi
            * self.photoionization_data["x_sect"]
            / self.photoionization_data["nu"]
            / const.h.cgs.value,
            axis=0,
        )
        corrected_photoion_coeff = corrected_photoion_coeff.multiply(
            stimulated_emission_correction
        )
        corrected_photoion_coeff.insert(
            0, "nu", self.photoionization_data["nu"]
        )
        corrected_photoion_coeff = corrected_photoion_coeff.groupby(
            level=[0, 1, 2]
        )
        tmp = {}
        for i in range(self.no_of_shells):
            tmp[i] = corrected_photoion_coeff.apply(
                lambda sub: trapezoid(sub[i], sub["nu"])
            )

        corrected_photoion_coeff = pd.DataFrame(tmp)
        return corrected_photoion_coeff

    def _calculate_rate_coefficient_from_estimator(self):
        index = self._get_estimator_index()
        lte_nonlte_level_pop_ratio = self._get_lte_nonlte_level_pop_ratio(index)
        corrected_photoion_coeff = (
            self.photo_ion_estimator
            - lte_nonlte_level_pop_ratio * self.stim_recomb_estimator
        )
        if (corrected_photoion_coeff < 0).any().any():
            raise PlasmaException(
                "Negative values in _calculate_rate_coefficient_from_estimator. Try raising the number of montecarlo packets."
            )
        corrected_photoion_coeff = pd.DataFrame(
            corrected_photoion_coeff,
            index=index,
            columns=np.arange(self.no_of_shells),
        )
        return corrected_photoion_coeff

    def _calculate_rate_coefficient_combination(
        self, rate_coeff_estimator, rate_coeff_dilute_bb, min_counts=100
    ):
        combined_rate_coeff = rate_coeff_estimator.where(
            self.estimators["photo_ion_statistics"] > min_counts,
            other=rate_coeff_dilute_bb,
        )
        return combined_rate_coeff

    def _calculate_j_nus(self):
        nus = self.photoionization_data["nu"].values
        j_nus = self.ws * intensity_black_body(nus[np.newaxis].T, self.t_rads)
        return pd.DataFrame(
            j_nus,
            index=self.photoionization_data.index,
            columns=np.arange(self.no_of_shells),
        )

    def _calculate_boltzmann_factor(self, nu):
        u0s = self._calculate_u0s(nu)
        return np.exp(-u0s)

    @staticmethod
    def _check_for_low_statistics(estimator_statistics, count_threshold=100):
        min_count = estimator_statistics.min()
        no_of_bad_elements = (estimator_statistics < count_threshold).sum()
        if no_of_bad_elements != 0:
            logger.warning(
                "{} MC estimators have been updated less than {} times, with a minimum of"
                " {} updates".format(
                    no_of_bad_elements, count_threshold, min_count
                )
            )
        return no_of_bad_elements

    def _calculate_stimulated_emission_correction(self):
        nu = self.photoionization_data["nu"].values
        boltzmann_factor = self._calculate_boltzmann_factor(nu)
        lte_nonlte_level_pop_ratio = self._get_lte_nonlte_level_pop_ratio(
            self.photoionization_data.index
        )
        correction_factor = 1.0 - lte_nonlte_level_pop_ratio * boltzmann_factor
        return correction_factor

    def _get_lte_nonlte_level_pop_ratio(self, index):
        level_pop_lte = self._get_lte_level_pop(index)
        level_pop = self._get_level_pop(index)
        continuum_pop = self._get_ion_number_density(index)
        continuum_pop_lte = self._get_lte_ion_number_density(index)
        ratio = (continuum_pop / continuum_pop_lte) * (
            level_pop_lte / level_pop
        )
        return ratio

    def _get_estimator_index(self):
        index = pd.MultiIndex.from_tuples(
            self.input.atom_data.continuum_data.multi_index_nu_sorted
        )
        index.names = ["atomic_number", "ion_number", "level_number"]
        return index

    @property
    def level_lower_energy(self):
        return self._get_level_energy(self.rate_coefficient.index)


class RadiativeRecombination(PhysicalContinuumProcess, BoundFreeEnergyMixIn):
    """
    Represents the process of radiative recombination.

    Attributes
    ----------
    input: `tardis.iip_plasma.continuum.input_data.ContinuumInputData`-object
        The common input data object.
    rate_coefficient: pd.DataFrame
        Multiplying the rate coefficient with the number densities of the interacting particles gives the rate
        per unit volume of the transition.
    cooling_rate: pd.DataFrame
        The rate per unit volume at which heat is radiated by spontaneous free-bound transitions.

    Class Attributes
    ----------------
    name: str
        The name used in setattr(object, name, value).
    cooling: bool
        True if the physical process contributes to the cooling of the plasma. Enables calculation of cooling_rate.
    """

    name = "radiative_recombination"

    def __init__(self, input_data):
        super(RadiativeRecombination, self).__init__(input_data)

    def _calculate_cooling_rate(self):
        sp_recombination_coeff_E = self._calculate_rate_coefficient(
            modified=True
        )
        fb_cooling_rate = sp_recombination_coeff_E - self.rate_coefficient
        fb_cooling_rate = fb_cooling_rate.multiply(
            const.h.cgs.value * self.nu_i, axis=0
        )
        fb_cooling_rate = fb_cooling_rate.multiply(
            self.electron_densities, axis=1
        )
        ion_number_density = self._get_ion_number_density(fb_cooling_rate.index)
        fb_cooling_rate = fb_cooling_rate.multiply(ion_number_density)
        continuum_edge_idx = self._get_continuum_edge_idx(fb_cooling_rate.index)
        fb_cooling_rate.set_index(continuum_edge_idx, inplace=True)

        continuum_edge_idx_alt = self._get_continuum_edge_idx(
            self.input.sp_fb_cooling_rates.index
        )
        fb_cooling_rate_alt = self.input.sp_fb_cooling_rates.set_index(
            continuum_edge_idx_alt
        )
        if len(fb_cooling_rate_alt) == 0:
            print("First iteration use usual fb cooling rate")
            fb_cooling_rate_alt = fb_cooling_rate
        # print "Fb cooling rate ratio"
        # print fb_cooling_rate_alt / fb_cooling_rate
        # return fb_cooling_rate
        return fb_cooling_rate_alt

    def _calculate_rate_coefficient(self, modified=False):
        """
        Calculates the rate coefficient for spontaneous recombination.

        Parameters
        ----------
        modified: bool, optional
            Switches between calculation of normal rate coefficient and a modified version, which
             is needed for calculating cooling rates.

        Returns
        -------
        recomb_coeff: pd.DataFrame
            The rate coefficient for spontaneous recombination.

        """
        if modified == False:
            recomb_coeff = (
                8
                * np.pi
                * self.photoionization_data["x_sect"]
                * (self.photoionization_data["nu"]) ** 2
                / (const.c.cgs.value) ** 2
            ).values
        else:
            recomb_coeff = (
                8
                * np.pi
                * self.photoionization_data["x_sect"]
                * (self.photoionization_data["nu"]) ** 3
                / (const.c.cgs.value) ** 2
            ).values

        recomb_coeff = recomb_coeff[:, np.newaxis]
        boltzmann_factor = np.exp(
            -self.photoionization_data.nu.values[np.newaxis].T
            / self.t_electrons
            * (const.h.cgs.value / const.k_B.cgs.value)
        )
        recomb_coeff = pd.DataFrame(
            boltzmann_factor * recomb_coeff,
            index=self.photoionization_data.index,
        )

        if modified:
            em = recomb_coeff.copy(deep=True)
            em.insert(0, "nu", self.photoionization_data["nu"])
            em = em.groupby(level=[0, 1, 2])
            em_tmp = []
            for i in range(self.no_of_shells):
                emissivities = em.apply(
                    lambda sub: pd.DataFrame(
                        cumulative_trapezoid(sub[i], sub["nu"], initial=0.0),
                        columns=[i],
                    )
                )
                emissivities.index = emissivities.index.droplevel(3)

                emissivities = (
                    emissivities / emissivities.groupby(level=[0, 1, 2]).max()
                )

                em_tmp.append(emissivities)
            self.bf_emissivities = pd.concat(em_tmp, axis=1)
            self.bf_emissivities_mc = self.bf_emissivities.loc[
                self.input.atom_data.continuum_data.multi_index_nu_sorted
            ]
            self.bf_emissivities_mc = np.asfortranarray(
                self.bf_emissivities_mc.values
            )
            # self.bf_emissivities_mc = self.bf_emissivities.values

        recomb_coeff = recomb_coeff.divide(self.electron_densities, axis=1)
        recomb_coeff.insert(0, "nu", self.photoionization_data["nu"])
        recomb_coeff = recomb_coeff.groupby(level=[0, 1, 2])
        tmp = {}
        for i in range(self.no_of_shells):
            tmp[i] = recomb_coeff.apply(
                lambda sub: trapezoid(sub[i], sub["nu"])
            )
            if modified == True:
                tmp[i] /= self.nu_i
        recomb_coeff = pd.DataFrame(tmp)
        recomb_coeff = recomb_coeff.multiply(
            self._get_lte_level_pop(recomb_coeff.index)
        )
        ion_number_density = self._get_ion_number_density(
            recomb_coeff.index, dtype="dataframe"
        )
        recomb_coeff = recomb_coeff.divide(ion_number_density.values)
        if not modified:
            recomb_coeff = self.input.alpha_sp
        return recomb_coeff

    @property
    def deactivation_probabilities(self):
        energy_difference = self.nu_i * const.h.cgs * units.erg.to(units.eV)
        return self.rate_coefficient.multiply(energy_difference, axis=0)

    @property
    def internal_jump_probabilities(self):
        energy_difference = self.nu_i * const.h.cgs * units.erg.to(units.eV)
        index0 = pd.MultiIndex.from_arrays(
            [
                self.rate_coefficient.index.get_level_values(0),
                self.rate_coefficient.index.get_level_values(1),
            ]
        )
        energy0 = (
            self.nu_0.loc[index0].values * const.h.cgs * units.erg.to(units.eV)
        )
        energy_lower = energy0 - energy_difference
        return self.rate_coefficient.multiply(energy_lower, axis=0)


class RadiativeExcitation(PhysicalContinuumProcess):
    name = "radiative_excitation"
    cooling = False
    macro_atom_transitions = "up"

    @property
    def internal_jump_probabilities(self):
        return (
            self.input.radiative_transition_probabilities_prep.loc[
                self.transition_up_filter
            ]
            * cconst.c_einstein
        )


class RadiativeDeexcitation(PhysicalContinuumProcess):
    name = "radiative_deexcitation"
    cooling = False
    macro_atom_transitions = "down"

    @property
    def internal_jump_probabilities(self):
        return (
            self.input.radiative_transition_probabilities_prep.loc[
                self.transition_down_filter
            ]
            * cconst.c_einstein
        )

    @property
    def deactivation_probabilities(self):
        filter = self.transition_deactivation_filter
        deactivation_probabilities = (
            self.input.radiative_transition_probabilities_prep.loc[filter]
            * cconst.c_einstein
        )
        deactivation_probabilities.insert(
            0, "lines_idx", self.macro_atom_data.loc[filter, "lines_idx"].values
        )
        return deactivation_probabilities


class FreeFree(PhysicalContinuumProcess):
    """
    Represents free-free transitions.

    Attributes
    ----------
    input: `tardis.iip_plasma.continuum.input_data.ContinuumInputData`-object
        The common input data object.
    cooling_rate: np.ndarray
        The rate per unit volume at which heat is converted into radiant energy by ff-emissions.
    chi_ff_factor: np.ndarray
        Used in the calculation of free-free opacities in the montecarlo run.

    Class Attributes
    ----------------
    name: str
        The name used in setattr(object, name, value).
    cooling: bool
        True if the physical process contributes to the cooling of the plasma. Enables calculation of cooling_rate.
    """

    name = "free_free"

    def __init__(self, input_data, **kwargs):
        super(FreeFree, self).__init__(input_data, **kwargs)
        self.chi_ff_factor = self._calculate_chi_ff_factor()

    def _calculate_cooling_rate(self, **kwargs):
        # TODO: value for Gaunt factor (Lucy: = 1; Osterbrock recommendation for nebular conditions: = 1.3 )
        factor = (
            self.ion_number_density.mul(
                np.square(self.input.ion_charges), axis=0
            )
            .sum()
            .values
        )
        cooling_rate = (
            cconst.C0_ff
            * self.electron_densities
            * np.sqrt(self.t_electrons)
            * factor
        )
        return cooling_rate

    def _calculate_chi_ff_factor(self):
        ionic_charge_squared = np.square(self._get_ionic_charge())
        ff_gaunt_factor = self._get_ff_gaunt_factor(
            self.ion_number_density.index
        )
        chi_ff_helper = (
            3.69255e8 * self.electron_densities / np.sqrt(self.t_electrons)
        )
        # Likely to be Eq. 6.1.8 in http://personal.psu.edu/rbc3/A534/lec6.pdf
        # see also FF_OPAC_CONST in tardis/opacities/opacities.py
        chi_ff_factor = (
            self.ion_number_density.multiply(
                ionic_charge_squared * ff_gaunt_factor, axis=0
            )
            .sum()
            .values
        )
        chi_ff_factor *= chi_ff_helper
        return chi_ff_factor

    def _get_ionic_charge(self):
        return self.ion_number_density.index.get_level_values(1).values

    def _get_ff_gaunt_factor(self, ion_index):
        return np.ones(ion_index.values.shape[0])
