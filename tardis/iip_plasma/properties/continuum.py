import logging

import numpy as np
import pandas as pd
from astropy import constants as const
from scipy.integrate import simpson, trapezoid
from scipy.interpolate import PchipInterpolator

from tardis.iip_plasma.properties.base import ProcessingPlasmaProperty
from tardis.util.base import intensity_black_body

logger = logging.getLogger(__name__)

__all__ = [
    "BfHeatingRateCoeff",
    "BfHeatingRateCoeffEstim",
    "CollDeexcRateCoeff",
    "CollExcCooling",
    "CollExcRateCoeff",
    "CollIonRateCoeff",
    "CollRecombRateCoeff",
    "Logger",
    "PhotoIonRateCoeff",
    "SpontRecombRateCoeff",
    "StimRecombCoolingRateCoeff",
    "StimRecombRateCoeff",
    "ThermalBalanceTest",
    "Yg",
    "YgInterpolator",
]


def get_estimator_index(photo_ion_index_sorted):
    index = pd.MultiIndex.from_tuples(photo_ion_index_sorted)
    index.names = ["atomic_number", "ion_number", "level_number"]
    return index


def calculate_rate_coefficient_from_estimator(
    estimator, photo_ion_index_sorted
):
    no_of_shells = estimator.shape[1]
    index = get_estimator_index(photo_ion_index_sorted)
    rate_coeff = pd.DataFrame(
        estimator, index=index, columns=np.arange(no_of_shells), copy=True
    )
    return rate_coeff


class SpontRecombRateCoeff(ProcessingPlasmaProperty):
    """
    Attributes
    ----------
    alpha_sp : Pandas DataFrame, dtype float
               The rate coefficient for spontaneous recombination.
    """

    outputs = ("alpha_sp",)
    latex_name = ("\\alpha^{\\textrm{sp}}",)
    # TODO
    latex_formula = ("",)

    def calculate(self, photo_ion_cross_sections, t_electrons, phi_lucy):
        # print "Calculate alpha_sp \n"
        x_sect = photo_ion_cross_sections["x_sect"]
        nu = photo_ion_cross_sections["nu"]

        alpha_sp = (
            8 * np.pi * x_sect * nu**2 / (const.c.cgs.value) ** 2
        ).values
        alpha_sp = alpha_sp[:, np.newaxis]
        boltzmann_factor = np.exp(
            -nu.values[np.newaxis].T
            / t_electrons
            * (const.h.cgs.value / const.k_B.cgs.value)
        )
        recomb_coeff = pd.DataFrame(boltzmann_factor * alpha_sp, index=nu.index)
        recomb_coeff.insert(0, "nu", nu)
        recomb_coeff = recomb_coeff.groupby(level=[0, 1, 2])
        tmp = {}
        for i in range(len(t_electrons)):
            tmp[i] = recomb_coeff.apply(
                lambda sub: trapezoid(sub[i], sub["nu"])
            )
        alpha_sp = pd.DataFrame(tmp)
        phi_lucy = phi_lucy.loc[alpha_sp.index]

        # TODO: Revert
        if self.plasma_parent.niter > self.plasma_parent.niter_ly:
            # print "Set alpha_sp to zero"
            alpha_sp.loc[(1, 0, 0)] = 0.0
        return alpha_sp.multiply(phi_lucy)


class PhotoIonRateCoeff(ProcessingPlasmaProperty):
    """
    Attributes
    ----------
    gamma : Pandas DataFrame, dtype float
               The rate coefficient for radiative ionization.
    """

    outputs = (
        "gamma",
        "gamma_Ly_estim",
        "gamma_Ly_analytic",
        "gamma_Ly_analytic_evolution",
        "gamma_Ly_estim_evolution",
    )
    latex_name = ("\\gamma", "", "")
    # TODO
    latex_formula = ("",)

    def calculate(
        self,
        photo_ion_estimator,
        photo_ion_index_sorted,
        previous_t_electrons,
        photo_ion_cross_sections,
        previous_b,
    ):
        # print "Calculate gamma \n"
        if photo_ion_estimator is not None:
            gamma = calculate_rate_coefficient_from_estimator(
                photo_ion_estimator, photo_ion_index_sorted
            )
            gamma_Ly_old = gamma.loc[(1, 0, 0)].copy(deep=True).values
            gamma_Ly_new = self._update_Lymann_continuum(
                gamma,
                photo_ion_cross_sections,
                previous_t_electrons,
                previous_b,
            )
            gamma_Ly_new = gamma_Ly_new[0]

            # gamma.loc[(1,0,0)] = gamma_Ly_new
            # print "\n"
            # print "gamma_anal/gamma_estim:"
            # print gamma_Ly_new / gamma_Ly_old
            # print "\n\n"
            self.gamma_Ly_analytic_evolution.append(gamma_Ly_new)
            self.gamma_Ly_estim_evolution.append(gamma_Ly_old)
            # print self.plasma_parent.niter, self.plasma_parent.niter_ly
            if self.plasma_parent.niter > self.plasma_parent.niter_ly:
                # print "Setting gamma to zero"
                gamma.loc[(1, 0, 0)] = 0.0
        else:
            gamma = None
            gamma_Ly_new = None
            gamma_Ly_old = None
            self.gamma_Ly_estim_evolution = []
            self.gamma_Ly_analytic_evolution = []
        return (
            gamma,
            gamma_Ly_old,
            gamma_Ly_new,
            self.gamma_Ly_analytic_evolution,
            self.gamma_Ly_estim_evolution,
        )

    def calculate_from_radiation_field_model(
        self, photo_ion_cross_sections, w, t_rad
    ):
        j_nus = self._calculate_j_nus(photo_ion_cross_sections, w, t_rad)
        photoion_coeff = j_nus.multiply(
            4.0
            * np.pi
            * photo_ion_cross_sections["x_sect"]
            / photo_ion_cross_sections["nu"]
            / const.h.cgs.value,
            axis=0,
        )
        photoion_coeff.insert(0, "nu", photo_ion_cross_sections["nu"])
        photoion_coeff = photoion_coeff.groupby(level=[0, 1, 2])
        tmp = {}
        for i in range(len(w)):
            tmp[i] = photoion_coeff.apply(
                lambda sub: trapezoid(sub[i], sub["nu"])
            )
        photoion_coeff = pd.DataFrame(tmp)

        # if photo_ion_estimator is not None:
        #    print "Bf_estim/bf:", photo_ion_estimator/photoion_coeff
        return photoion_coeff

    def _update_Lymann_continuum(
        self, gamma, photo_ion_cross_sections, t_electrons, previous_b
    ):
        x_sect_Ly = photo_ion_cross_sections.loc[(1, 0, 0)]
        # print "Previous b:", previous_b.loc[(1,0,0)]
        W = (previous_b.loc[(1, 0, 0)].values) ** -1
        gamma_Ly = self.calculate_from_radiation_field_model(
            x_sect_Ly, W, t_electrons
        )
        return gamma_Ly.values

    @staticmethod
    def _calculate_j_nus(photoionization_data, ws, t_rads):
        nus = photoionization_data["nu"].values
        j_nus = ws * intensity_black_body(nus[np.newaxis].T, t_rads)
        return pd.DataFrame(j_nus, index=photoionization_data.index)


class BfHeatingRateCoeffEstim(ProcessingPlasmaProperty):
    """
    Attributes
    ----------
    gamma : Pandas DataFrame, dtype float
               The rate coefficient for bound-free heating.
    """

    outputs = ("bf_heating_coeff_estim",)
    latex_name = ("\\h_{\\textrm{bf}}}",)
    # TODO
    latex_formula = ("",)

    def calculate(self, bf_heating_estimator, photo_ion_index_sorted):
        if bf_heating_estimator is not None:
            bf_heating_coeff = calculate_rate_coefficient_from_estimator(
                bf_heating_estimator, photo_ion_index_sorted
            )
            if self.plasma_parent.niter > self.plasma_parent.niter_ly:
                bf_heating_coeff.loc[(1, 0, 0)] = 0.0
        else:
            bf_heating_coeff = None
        return bf_heating_coeff


class BfHeatingRateCoeff(ProcessingPlasmaProperty):
    """
    Attributes
    ----------
    gamma : Pandas DataFrame, dtype float
               The rate coefficient for bound-free heating.
    """

    outputs = ("bf_heating_coeff",)
    latex_name = ("\\h_{\\textrm{bf}}}",)
    # TODO
    latex_formula = ("",)

    def calculate(self, bf_heating_estimator, photo_ion_index_sorted):
        if bf_heating_estimator is not None:
            bf_heating_coeff = calculate_rate_coefficient_from_estimator(
                bf_heating_estimator, photo_ion_index_sorted
            )
            if self.plasma_parent.niter > self.plasma_parent.niter_ly:
                bf_heating_coeff.loc[(1, 0, 0)] = 0.0
        else:
            bf_heating_coeff = None
        return bf_heating_coeff

    def calculate_from_radiation_field_model(
        self, photo_ion_cross_sections, w, t_rad
    ):
        j_nus = self._calculate_j_nus(photo_ion_cross_sections, w, t_rad)
        nu_i = photo_ion_cross_sections["nu"].groupby(level=[0, 1, 2]).first()
        nu_is = nu_i.loc[photo_ion_cross_sections.index]
        photoion_heating_coeff = j_nus.multiply(
            4.0
            * np.pi
            * photo_ion_cross_sections["x_sect"]
            * (1 - nu_is / photo_ion_cross_sections["nu"]),
            axis=0,
        )
        photoion_heating_coeff.insert(0, "nu", photo_ion_cross_sections["nu"])
        photoion_heating_coeff = photoion_heating_coeff.groupby(level=[0, 1, 2])
        tmp = {}
        for i in range(len(w)):
            tmp[i] = photoion_heating_coeff.apply(
                lambda sub: simpson(sub[i], sub["nu"], even="first")
            )
        photoion_heating_coeff = pd.DataFrame(tmp)

        return photoion_heating_coeff

    @staticmethod
    def _calculate_j_nus(photoionization_data, ws, t_rads):
        nus = photoionization_data["nu"].values
        j_nus = ws * intensity_black_body(nus[np.newaxis].T, t_rads)
        return pd.DataFrame(j_nus, index=photoionization_data.index)


class StimRecombCoolingRateCoeff(ProcessingPlasmaProperty):
    """
    Attributes
    ----------
    gamma : Pandas DataFrame, dtype float
               The rate coefficient for stimulated recombination cooling.
    """

    outputs = ("stim_recomb_cooling_coeff",)
    latex_name = ("\\alpha^{\\textrm{stim, E}}}",)
    # TODO
    latex_formula = ("",)

    def calculate(self, stim_recomb_cooling_estimator, photo_ion_index_sorted):
        if stim_recomb_cooling_estimator is not None:
            coeff = calculate_rate_coefficient_from_estimator(
                stim_recomb_cooling_estimator, photo_ion_index_sorted
            )
            if self.plasma_parent.niter > self.plasma_parent.niter_ly:
                coeff.loc[(1, 0, 0)] = 0.0
        else:
            coeff = None
        return coeff


class StimRecombRateCoeff(ProcessingPlasmaProperty):
    """
    Attributes
    ----------
    gamma : Pandas DataFrame, dtype float
               The rate coefficient for radiative ionization.
    """

    outputs = ("alpha_stim",)
    latex_name = ("\\alpha^{\\textrm{stim}}",)
    # TODO
    latex_formula = ("",)

    def calculate(
        self, stim_recomb_estimator, photo_ion_index_sorted, phi_lucy
    ):
        # print "Calculate alpha_stim \n"
        if stim_recomb_estimator is not None:
            alpha_stim = calculate_rate_coefficient_from_estimator(
                stim_recomb_estimator, photo_ion_index_sorted
            )
            phi_lucy = phi_lucy.loc[alpha_stim.index]
            alpha_stim = alpha_stim.multiply(phi_lucy)
            if self.plasma_parent.niter > self.plasma_parent.niter_ly:
                alpha_stim.loc[(1, 0, 0)] = 0.0
        else:
            alpha_stim = None
        return alpha_stim

    def calculate_from_radiation_field_model(
        self,
        photo_ion_cross_sections,
        w,
        t_rad,
        stim_recomb_estimator,
        t_electrons,
        phi_lucy,
    ):
        j_nus = self._calculate_j_nus(photo_ion_cross_sections, w, t_rad)
        photoion_coeff = j_nus.multiply(
            4.0
            * np.pi
            * photo_ion_cross_sections["x_sect"]
            / photo_ion_cross_sections["nu"]
            / const.h.cgs.value,
            axis=0,
        )
        boltzmann_factor = np.exp(
            -photo_ion_cross_sections.nu.values[np.newaxis].T
            / t_electrons
            * (const.h.cgs.value / const.k_B.cgs.value)
        )
        photoion_coeff = photoion_coeff * boltzmann_factor
        photoion_coeff.insert(0, "nu", photo_ion_cross_sections["nu"])
        photoion_coeff = photoion_coeff.groupby(level=[0, 1, 2])
        tmp = {}
        for i in range(len(w)):
            tmp[i] = photoion_coeff.apply(
                lambda sub: trapezoid(sub[i], sub["nu"])
            )
        photoion_coeff = pd.DataFrame(tmp)

        # if stim_recomb_estimator is not None:
        #    print "stim_recomb_estim/stim_recomb:", stim_recomb_estimator/photoion_coeff

        phi_lucy = phi_lucy.loc[photoion_coeff.index]
        alpha_stim = photoion_coeff.multiply(phi_lucy)

        return alpha_stim

    @staticmethod
    def _calculate_j_nus(photoionization_data, ws, t_rads):
        nus = photoionization_data["nu"].values
        j_nus = ws * intensity_black_body(nus[np.newaxis].T, t_rads)
        return pd.DataFrame(j_nus, index=photoionization_data.index)

    # def _calculate_boltzmann_factor(self, nu):
    #     u0s = self._calculate_u0s(nu)
    #     return np.exp(-u0s)
    #
    # def _calculate_u0s(self, nu, t_electrons):
    #     u0s = nu[np.newaxis].T / t_electrons * (const.h.cgs.value / const.k_B.cgs.value)
    #     return u0s


class CollIonRateCoeff(ProcessingPlasmaProperty):
    """
    Attributes
    ----------
    coll_ion_coeff : Pandas DataFrame, dtype float
               The rate coefficient for collisional ionization.
               Multiply with the electron density and the level number density
               to obtain the total rate.
    """

    outputs = ("coll_ion_coeff",)
    latex_name = ("\\c_{\\textrm{i,}\\kappa}",)
    # TODO
    latex_formula = ("",)

    def calculate(self, photo_ion_cross_sections, t_electrons):
        """
        Calculates the rate coefficient for ionization by thermal electrons using the Seaton
        approximation (Equation 5-79 of Mihalas 1978).

        Returns
        -------
        recomb_coeff: pd.DataFrame
            The rate coefficient for collisional ionization.

        """
        # print "Calculate coll_ion \n"
        collion_coeff = photo_ion_cross_sections.groupby(
            level=[0, 1, 2]
        ).first()
        factor = self._calculate_factor(
            collion_coeff["nu"].values, t_electrons, index=collion_coeff.index
        )
        collion_coeff = 1.55e13 * collion_coeff["x_sect"]
        ion_number = collion_coeff.index.get_level_values(1).values
        collion_coeff[ion_number == 0] *= 0.1
        collion_coeff[ion_number == 1] *= 0.2
        collion_coeff[ion_number >= 2] *= 0.3
        collion_coeff = factor.multiply(collion_coeff, axis=0)
        collion_coeff = collion_coeff.divide(np.sqrt(t_electrons), axis=1)
        return collion_coeff

    def _calculate_factor(self, nu_ijk, t_electrons, index):
        u0s = self._calculate_u0s(nu_ijk, t_electrons)
        factor = 1.0 / u0s * np.exp(-u0s)
        factor = pd.DataFrame(
            factor, index=index, columns=np.arange(len(t_electrons))
        )
        return factor

    @staticmethod
    def _calculate_u0s(nu, t_electrons):
        u0s = (
            nu[np.newaxis].T
            / t_electrons
            * (const.h.cgs.value / const.k_B.cgs.value)
        )
        return u0s


class CollRecombRateCoeff(ProcessingPlasmaProperty):
    """
    Attributes
    ----------
    coll_recomb_coeff : Pandas DataFrame, dtype float
               The rate coefficient for collisional recombination.
               Multiply with the electron density squared and the ion number density
               to obtain the total rate.
    """

    outputs = ("coll_recomb_coeff",)
    latex_name = ("\\c_{\\kappa\\textrm{i,}}",)
    # TODO
    latex_formula = ("",)

    def calculate(self, phi_lucy, coll_ion_coeff):
        # print "Calculate coll_recomb \n"
        return coll_ion_coeff.multiply(phi_lucy.loc[coll_ion_coeff.index])


class YgInterpolator(ProcessingPlasmaProperty):
    outputs = ("Yg_interp", "yg_allowed_index", "yg_forbidden_index")
    latex_name = ("Yg_interp",)

    def calculate(self, yg_data, T_Yg, lines_multi_index):
        allowed_mask = yg_data.index.isin(lines_multi_index)
        yg_allowed_index = yg_data.index[allowed_mask]
        yg_forbidden_index = yg_data.index[~allowed_mask]
        interp = PchipInterpolator(T_Yg, yg_data, axis=1, extrapolate=True)
        return interp, yg_allowed_index, yg_forbidden_index


class Yg(ProcessingPlasmaProperty):
    outputs = ("yg",)
    latex_name = ("Y/g",)

    def calculate(self, Yg_interp, Yg_index, t_electrons):
        yg = Yg_interp(t_electrons)
        return pd.DataFrame(yg, index=Yg_index)


class CollExcRateCoeff(ProcessingPlasmaProperty):
    outputs = ("coll_exc_coeff", "nu_lines_coll")
    latex_name = ("\\c_{\\kappa\\textrm{i,}}", "\\nu_\\textrm{lines}_coll")
    # TODO
    latex_formula = ("",)

    def calculate(
        self,
        lines,
        excitation_energy,
        t_electrons,
        yg,
        yg_allowed_index,
        yg_forbidden_index,
    ):
        I_H = cconst.I_H

        # lines_filtered = lines.query('atomic_number == 1')
        lines_filtered = lines
        f_lu = lines_filtered.f_lu.values
        nu_lines = lines_filtered.nu.values

        coll_excitation_coeff = pd.DataFrame(
            14.5
            * cconst.c0_reg
            * f_lu
            * (I_H / (const.h.cgs.value * nu_lines)) ** 2
        )
        coll_excitation_coeff = coll_excitation_coeff.dot(
            pd.DataFrame(np.sqrt(t_electrons)).T
        )
        u0 = (
            nu_lines[np.newaxis].T
            / t_electrons
            * (const.h.cgs.value / const.k_B.cgs.value)
        )
        gamma = 0.276 * self.exp1_times_exp(u0)
        gamma[gamma < 0.2] = 0.2
        factor = pd.DataFrame(u0 * np.exp(-u0) * gamma)
        coll_excitation_coeff = coll_excitation_coeff.multiply(factor, axis=0)

        coll_excitation_coeff.set_index(lines_filtered.index, inplace=True)

        # Calculate rates based on PB04
        ind = yg.index
        ll_index = pd.MultiIndex.from_arrays(
            [
                ind.get_level_values(0),
                ind.get_level_values(1),
                ind.get_level_values(2),
            ]
        )
        lu_index = pd.MultiIndex.from_arrays(
            [
                ind.get_level_values(0),
                ind.get_level_values(1),
                ind.get_level_values(3),
            ]
        )
        e_lu = excitation_energy.loc[lu_index]
        e_ll = excitation_energy.loc[ll_index]
        hnu = e_lu.values - e_ll.values
        hnu = pd.DataFrame(hnu, index=ind)
        # hnu = pd.DataFrame(const.h.cgs.value * nu_lines, index=index)
        # boltzman_factor = pd.DataFrame(np.exp(-hnu.values / (const.k_B.cgs.value
        #    * t_electrons.T)), index=hnu.index)
        boltzman_factor = pd.DataFrame(
            np.exp(-hnu.values / (const.k_B.cgs.value * t_electrons.T)),
            index=yg.index,
        )
        q_ij = (
            8.629e-6 / np.sqrt(t_electrons) * yg * boltzman_factor.loc[yg.index]
        )
        # coll_excitation_coeff.loc[q_ij.index] = q_ij
        coll_excitation_coeff.loc[yg_allowed_index] = q_ij.loc[yg_allowed_index]
        # coll_excitation_coeff = pd.concat(
        #    [coll_excitation_coeff, q_ij.loc[yg_forbidden_index]]
        # )
        # TODO: Also use forbidden transitions

        return coll_excitation_coeff, nu_lines

    @staticmethod
    def exp1_times_exp(x):
        """
        Product of the Exponential integral E1 and an exponential.
        This function calculates the product of the Exponential integral E1
        and an exponential in a way that also works for large values.

        Parameters
        ----------
        x : array_like
            Input values.

        Returns
        -------
        array_like
            Output array.
        """
        f = exp1(x) * np.exp(x)
        # Use Laurent series for large values to avoid infinite exponential
        mask = x > 500
        f[mask] = (x**-1 - x**-2 + 2 * x**-3 - 6 * x**-4)[mask]
        return f


class CollDeexcRateCoeff(ProcessingPlasmaProperty):
    outputs = ("coll_deexc_coeff",)

    def calculate(self, lte_level_boltzmann_factor_Te, coll_exc_coeff):
        atom_number = coll_exc_coeff.index.get_level_values(0).values
        ion_number = coll_exc_coeff.index.get_level_values(1).values
        level_number_lower = coll_exc_coeff.index.get_level_values(2).values
        level_number_upper = coll_exc_coeff.index.get_level_values(3).values
        index_lower = pd.MultiIndex.from_arrays(
            [atom_number, ion_number, level_number_lower]
        )
        index_upper = pd.MultiIndex.from_arrays(
            [atom_number, ion_number, level_number_upper]
        )

        n_lower_prop = lte_level_boltzmann_factor_Te.loc[index_lower].values
        n_upper_prop = lte_level_boltzmann_factor_Te.loc[index_upper].values
        lte_level_pop_ratio = n_lower_prop / n_upper_prop

        coll_deexc_coeff = coll_exc_coeff * lte_level_pop_ratio
        return coll_deexc_coeff


class CollExcCooling(ProcessingPlasmaProperty):
    outputs = ("coll_exc_cooling", "coll_deexc_heating")

    def calculate(
        self,
        coll_exc_coeff,
        coll_deexc_coeff,
        level_number_density,
        electron_densities,
        nu_lines_coll,
    ):
        coll_exc_cooling_coeff = coll_exc_coeff.multiply(
            nu_lines_coll * const.h.cgs.value, axis=0
        )
        coll_exc_cooling_coeff = coll_exc_cooling_coeff.multiply(
            electron_densities, axis=1
        )

        coll_deexc_heating_coeff = coll_deexc_coeff.multiply(
            nu_lines_coll * const.h.cgs.value, axis=0
        )
        coll_deexc_heating_coeff = coll_deexc_heating_coeff.multiply(
            electron_densities, axis=1
        )

        lower_index = coll_exc_coeff.index.droplevel(level=3)
        upper_index = coll_exc_coeff.index.droplevel(level=2)

        coll_exc_cooling = coll_exc_cooling_coeff.multiply(
            level_number_density.loc[lower_index].values
        )
        coll_deexc_heating = coll_deexc_heating_coeff.multiply(
            level_number_density.loc[upper_index].values
        )

        return coll_exc_cooling.sum(), coll_deexc_heating.sum()


from tardis.iip_plasma.properties.general import BetaRadiation, GElectron
from tardis.iip_plasma.properties.ion_population import PhiSahaElectrons
from tardis.iip_plasma.properties.level_population import PhiLucy
from tardis.iip_plasma.properties.partition_function import (
    LevelBoltzmannFactorLTETe,
    PartitionFunction,
)
from scipy.special import exp1

from tardis.iip_plasma.continuum.constants import continuum_constants as cconst
from tardis.iip_plasma.continuum.util import get_ion_multi_index


class ThermalBalanceTest(ProcessingPlasmaProperty):
    outputs = ("fractional_heating", "sp_fb_cooling_rates")
    latex_name = ("",)
    # TODO
    latex_formula = ("",)

    def __init__(self, plasma_parent, T_min=2500, T_max=1.0e5):
        super(ThermalBalanceTest, self).__init__(plasma_parent)
        self.T_min = T_min
        self.T_max = T_max
        self._t_diff = None
        self._counter = -1
        # self._iterations = -1
        self.max_it = 25
        self.t_electrons_evolution = None
        self.t_electrons_evolution_all = []
        self.frac_heat_all = []
        self.level_number_previous = None
        self.heating_rates = None
        self._converged_t_electrons = []

    def calculate(
        self,
        photo_ion_cross_sections,
        t_rad,
        t_electrons,
        previous_t_electrons,
        excitation_energy,
        g,
        levels,
        ionization_data,
        ff_heating_estimator,
        bf_heating_coeff,
        stim_recomb_cooling_coeff,
        coll_deexc_heating_estimator,
        ion_number_density,
        level_number_density,
        electron_densities,
        lines,
        w,
        time_explosion,
        previous_b,
        coll_exc_cooling,
        coll_deexc_heating,
    ):
        # new_t_electrons = np.zeros_like(electron_densities)

        # self._iterations += 1
        # print 'Iterations:', self._iterations
        fractional_heating = np.zeros_like(t_electrons)
        sp_fb_cooling_rates = pd.DataFrame(columns=ion_number_density.columns)
        if ff_heating_estimator is not None:
            # if self.t_electrons_evolution is None:
            # self.t_electrons_evolution = np.zeros((self.max_it + 2, len(new_t_electrons)))
            # self.heating_rates = np.zeros((self.max_it + 2, len(new_t_electrons)))
            # self.coll_exc_heating_evolution = np.zeros((self.max_it + 2, len(new_t_electrons)))
            # self.coll_t_electrons = np.zeros((self.max_it + 2, len(new_t_electrons)))
            # self.delta_heat_estim = np.zeros((self.max_it + 2, len(new_t_electrons)))
            # self.delta_heat = np.zeros((self.max_it + 2, len(new_t_electrons)))
            # self.frac_heat = np.zeros((self.max_it + 2, len(new_t_electrons)))

            # conv_status = np.ones_like(t_electrons) * 0.8
            # if self._t_diff is None:
            #    self._t_diff = np.ones_like(t_electrons) * 100.

            for shell in range(len(electron_densities)):
                t_old = t_electrons[shell]
                # if not self._conv_status[shell]:
                args = (
                    ff_heating_estimator,
                    bf_heating_coeff,
                    stim_recomb_cooling_coeff,
                    coll_deexc_heating_estimator,
                    electron_densities,
                    ion_number_density,
                    level_number_density,
                    excitation_energy,
                    g,
                    levels,
                    ionization_data,
                    photo_ion_cross_sections,
                    lines,
                    shell,
                    t_rad,
                    w,
                    time_explosion,
                    previous_b,
                    previous_t_electrons,
                    coll_exc_cooling,
                    coll_deexc_heating,
                )

                heat_mid, frac_heating, sp_fb_cooling_rate = (
                    self.heating_function(t_old, *args)
                )
                fractional_heating[shell] = frac_heating
                sp_fb_cooling_rates[shell] = sp_fb_cooling_rate

        # logger.info("Fractional Heating: {}".format(fractional_heating[::2]))
        return fractional_heating, sp_fb_cooling_rates
        # gradient = (heat_max - heat_min)/(2. * self._t_diff[shell])
        # delta_T = - 0.4 * heat_mid / gradient
        # if np.fabs(delta_T) > 0.1 * t_old:
        #    delta_T = 0.1 * t_old * np.sign(delta_T)
        # delta_T = - 0.10 * heat_mid / gradient
        # if np.fabs(delta_T) > 200.:
        #    delta_T = 200. * np.sign(delta_T)
        # conv_status[shell] = np.fabs(delta_T)/t_old

        # delta_heat_estim = self.heating_function(t_old + delta_T, *args)[0] - heat_mid
        # self.delta_heat_estim[self._counter][shell] = delta_heat_estim
        # self.frac_heat[self._counter][shell] = frac_heating

        # if self._counter > 0:
        #    delta_heat = heat_mid - self.heating_rates[self._counter - 1][shell]
        #    self.delta_heat[self._counter - 1][shell] = delta_heat

        # else:
        #    delta_T = 0.
        #   conv_status[shell] = 0.0
        # heat_mid = 0.0

        # self.heating_rates[self._counter][shell] = heat_mid
        # new_t_electrons[shell] = t_old + delta_T
        # if new_t_electrons[shell] < self.T_min:
        #    new_t_electrons[shell] = self.T_min

        # self._t_diff[shell] = 0.1 * np.fabs(delta_T)

        # updated_link_t_rad_t_electron = new_t_electrons / t_rad

        # self.t_electrons_evolution[self._counter] = new_t_electrons
        # self._conv_status = conv_status < 1.e-4
        # self._conv_status = conv_status < 1.e-4

        # if self.level_number_previous is not None:
        #    l_conv = np.fabs(level_number_density - self.level_number_previous)/level_number_density
        #    print l_conv.max()
        # self.level_number_previous = level_number_density

        # logger.info("Thermal balance conv: {} it: {}".format(conv_status[::2], self._counter))
        # logger.info("new_T: {}".format(new_t_electrons[::2]))
        # if not np.all(conv_status < 1.e-4):
        #    if not self._counter > self.max_it:
        # self._counter += 1
        #        self.plasma_parent.update(link_t_rad_t_electron=new_link_t_rad_t_electron)

        # self._t_diff = None
        # self._counter = 0
        # self.level_number_previous = None

        # TODO: Solve this in a better way
        # if new_t_electrons.sum() == 0:
        #    print "Zero t_electrons"
        #    updated_link_t_rad_t_electron = t_electrons / t_rad

        # if self.t_electrons_evolution is not None:
        #    if (self.t_electrons_evolution>0).sum():
        #        self.t_electrons_evolution_all.append(self.t_electrons_evolution)
        #        self.t_electrons_evolution = np.zeros_like(self.t_electrons_evolution)
        # if hasattr(self, 'frac_heat'):
        #    if (self.frac_heat>0).sum():
        #        self.frac_heat_all.append(self.frac_heat)
        #        self.frac_heat = np.zeros_like(self.frac_heat)

        # self._conv_status = np.zeros(t_electrons.shape, dtype=bool)

        # self._converged_t_electrons.append(new_t_electrons)

    @staticmethod
    def _calculate_phi_lucy(
        t_electrons, excitation_energy, g, levels, ionization_data
    ):
        beta_electron = BetaRadiation.calculate(np.array([t_electrons]))
        level_boltzmann_factor = LevelBoltzmannFactorLTETe.calculate(
            excitation_energy, g, beta_electron, levels
        )
        partition_function = PartitionFunction.calculate(level_boltzmann_factor)
        g_electron_Te = GElectron.calculate(beta_electron)
        phi_Te = PhiSahaElectrons.calculate(
            g_electron_Te, beta_electron, partition_function, ionization_data
        )
        phi_lucy = PhiLucy(plasma_parent=None).calculate(
            phi_Te, level_boltzmann_factor, partition_function
        )
        return phi_lucy[0]

    def _calculate_ff_heating_balance(
        self,
        t_electron,
        ff_heating_estimator,
        electron_densities,
        ion_number_density,
        shell,
    ):
        ff_heating_factor = self._calculate_ff_heating_factor(
            ion_number_density[shell], electron_densities[shell]
        )
        # TODO
        # heating_balance = (ff_heating_estimator[shell] / np.sqrt(t_electron) - cconst.C0_ff * np.sqrt(t_electron)) \
        #                  * ff_heating_factor
        ff_heating = (
            ff_heating_estimator[shell]
            / np.sqrt(t_electron)
            * ff_heating_factor
        )
        ff_cooling = cconst.C0_ff * np.sqrt(t_electron) * ff_heating_factor
        return ff_heating, ff_cooling
        # return heating_balance

    @staticmethod
    def _calculate_ff_heating_factor(ion_number_density, electron_densities):
        ionic_charge_squared = np.square(
            ion_number_density.index.get_level_values(1).values
        )
        factor = (
            electron_densities
            * ion_number_density.multiply(ionic_charge_squared, axis=0).sum()
        )
        return factor

    def _calculate_bf_heating_rate(
        self,
        bf_heating_coeff,
        level_number_density,
        shell,
        t_electron,
        photo_ion_cross_sections,
        b,
        previous_t_electrons,
    ):
        bf_heating = bf_heating_coeff[shell]
        # self._update_Lymann_continuum(bf_heating, photo_ion_cross_sections, previous_t_electrons[shell], b.loc[(1,0,0), shell])
        return bf_heating.multiply(
            level_number_density.loc[bf_heating.index, shell]
        ).sum()

    def _calculate_adiabatic_cooling(
        self, electron_density, t_electron, time_explosion, shell
    ):
        # import pdb; pdb.set_trace()
        # v_inner = self.plasma_parent.v_inner[shell]
        # v_outer = self.plasma_parent.v_outer[shell]
        # volume = 4. * np.pi * (v_outer.value ** 3. - v_inner.value ** 3.) * time_explosion **3. / 3.
        ad_cool = (
            3.0
            * electron_density
            * const.k_B.cgs.value
            * t_electron
            / time_explosion
        )
        return ad_cool

    def _calculate_fb_cooling_rate(
        self,
        t_electron,
        stim_recomb_cooling_coeff,
        phi_lucy,
        electron_densities,
        ion_number_density,
        photo_ion_cross_sections,
        shell,
    ):
        alpha_E_sp = self._calculate_sp_recomb_heating_rate_coeff(
            t_electron, photo_ion_cross_sections
        )

        if self.plasma_parent.niter > self.plasma_parent.niter_ly:
            alpha_E_sp.loc[(1, 0, 0)] = 0.0

        recomb_cooling_coeff = alpha_E_sp
        level_cooling_rates = (
            recomb_cooling_coeff.multiply(
                phi_lucy.loc[recomb_cooling_coeff.index]
            )
            * electron_densities[shell]
        )
        ion_index = get_ion_multi_index(level_cooling_rates.index)
        level_cooling_rates = level_cooling_rates.multiply(
            ion_number_density.loc[ion_index, shell].values
        )
        ly_cooling_rate = level_cooling_rates.loc[(1, 0, 0)]
        cooling_rate = level_cooling_rates.sum()

        st_recomb_cooling_coeff = stim_recomb_cooling_coeff[shell]
        # st_recomb_cooling_coeff = self._calculate_stim_recomb_heating_rate_coeff(t_electron, photo_ion_cross_sections, t_rad[shell], w[shell])
        st_level_cooling_rates = (
            st_recomb_cooling_coeff.multiply(
                phi_lucy.loc[st_recomb_cooling_coeff.index]
            )
            * electron_densities[shell]
        )
        ion_index = get_ion_multi_index(st_level_cooling_rates.index)
        st_level_cooling_rates = st_level_cooling_rates.multiply(
            ion_number_density.loc[ion_index, shell].values
        )
        st_cooling_rate = st_level_cooling_rates.sum()

        return (
            cooling_rate + st_cooling_rate,
            ly_cooling_rate,
            level_cooling_rates,
        )

    def _calculate_coll_ion_heating_balance(
        self,
        t_electron,
        photo_ion_cross_sections,
        phi_lucy,
        electron_densities,
        level_number_density,
        ion_number_density,
        shell,
    ):
        electron_density = electron_densities[shell]
        coll_ion_rate_coeff = CollIonRateCoeff(plasma_parent=None).calculate(
            photo_ion_cross_sections, np.array([t_electron])
        )[0]
        index = coll_ion_rate_coeff.index
        ion_index = get_ion_multi_index(index)
        ion_density = ion_number_density.loc[ion_index, shell].values
        # heating_balance =\
        #    electron_density * coll_ion_rate_coeff * (phi_lucy.loc[index] * electron_density * ion_density -
        #        level_number_density.loc[index, shell])

        coll_ion_heating = (
            electron_density
            * coll_ion_rate_coeff
            * phi_lucy.loc[index]
            * electron_density
            * ion_density
        )

        coll_ion_cooling = (
            level_number_density.loc[index, shell]
            * electron_density
            * coll_ion_rate_coeff
        )
        nu_i = photo_ion_cross_sections["nu"].groupby(level=[0, 1, 2]).first()
        # heating_balance = (heating_balance * nu_i * const.h.cgs.value).sum()
        # return heating_balance
        return (coll_ion_heating * nu_i * const.h.cgs.value).sum(), (
            coll_ion_cooling * nu_i * const.h.cgs.value
        ).sum()

    def _calculate_sp_recomb_heating_rate_coeff(
        self, t_electron, photo_ion_cross_sections
    ):
        x_sect = photo_ion_cross_sections["x_sect"]
        nu = photo_ion_cross_sections["nu"]
        nu_i = photo_ion_cross_sections["nu"].groupby(level=[0, 1, 2]).first()
        alpha_sp_E = (
            8
            * np.pi
            * x_sect
            * nu**3
            * const.h.cgs.value
            / (const.c.cgs.value) ** 2
        )
        alpha_sp_E = alpha_sp_E.multiply(1.0 - nu_i / nu)
        boltzmann_factor = np.exp(
            -nu.values / t_electron * (const.h.cgs.value / const.k_B.cgs.value)
        )
        alpha_sp_E = pd.DataFrame(alpha_sp_E.multiply(boltzmann_factor))
        alpha_sp_E.insert(0, "nu", nu)
        alpha_sp_E = alpha_sp_E.groupby(level=[0, 1, 2]).apply(
            lambda sub: trapezoid(sub[0], sub["nu"])
        )
        return alpha_sp_E

    def _calculate_stim_recomb_heating_rate_coeff(
        self, t_electron, photo_ion_cross_sections, t_rad, w
    ):
        nus = photo_ion_cross_sections["nu"].values
        j_nus = w * intensity_black_body(nus, t_rad)
        x_sect = photo_ion_cross_sections["x_sect"]
        nu = photo_ion_cross_sections["nu"]
        nu_i = photo_ion_cross_sections["nu"].groupby(level=[0, 1, 2]).first()
        alpha_sp_E = (4.0 * np.pi * x_sect) * j_nus
        alpha_sp_E = alpha_sp_E.multiply(1.0 - nu_i / nu)
        boltzmann_factor = np.exp(
            -nu.values / t_electron * (const.h.cgs.value / const.k_B.cgs.value)
        )
        alpha_sp_E = pd.DataFrame(alpha_sp_E.multiply(boltzmann_factor))
        alpha_sp_E.insert(0, "nu", nu)
        alpha_sp_E = alpha_sp_E.groupby(level=[0, 1, 2]).apply(
            lambda sub: trapezoid(sub[0], sub["nu"])
        )
        return alpha_sp_E

    def heating_function(
        self,
        t_electron,
        ff_heating_estimator,
        bf_heating_coeff,
        stim_recomb_cooling_coeff,
        coll_deexc_heating_estimator,
        electron_densities,
        ion_number_density,
        level_number_density,
        excitation_energy,
        g,
        levels,
        ionization_data,
        photo_ion_cross_sections,
        lines,
        shell,
        t_rad,
        w,
        time_explosion,
        b,
        previous_t_electrons,
        coll_exc_cooling,
        coll_deexc_heating,
        mode="detailed",
    ):
        bf_heating = self._calculate_bf_heating_rate(
            bf_heating_coeff,
            level_number_density,
            shell,
            t_electron,
            photo_ion_cross_sections,
            b,
            previous_t_electrons,
        )

        # TODO: reuse these quantities in calculation of phi_lucy
        # beta_electron = BetaRadiation.calculate(np.array([t_electron]))
        # level_boltzmann_factor_Te = LevelBoltzmannFactorLTETe.calculate(excitation_energy, g, beta_electron, levels)

        phi_lucy = self._calculate_phi_lucy(
            t_electron, excitation_energy, g, levels, ionization_data
        )

        # ff_heating_balance = self._calculate_ff_heating_balance(
        #        t_electron, ff_heating_estimator, electron_densities, ion_number_density, shell)
        ff_heating, ff_cooling = self._calculate_ff_heating_balance(
            t_electron,
            ff_heating_estimator,
            electron_densities,
            ion_number_density,
            shell,
        )

        fb_cooling, ly_cooling_rate, sp_fb_cooling_rate = (
            self._calculate_fb_cooling_rate(
                t_electron,
                stim_recomb_cooling_coeff,
                phi_lucy,
                electron_densities,
                ion_number_density,
                photo_ion_cross_sections,
                shell,
            )
        )

        # coll_ion_heating_balance = \
        #    self._calculate_coll_ion_heating_balance(t_electron, photo_ion_cross_sections, phi_lucy,
        #        electron_densities, level_number_density, ion_number_density, shell)

        coll_ion_heating, coll_ion_cooling = (
            self._calculate_coll_ion_heating_balance(
                t_electron,
                photo_ion_cross_sections,
                phi_lucy,
                electron_densities,
                level_number_density,
                ion_number_density,
                shell,
            )
        )

        # coll_deexc_heating = coll_deexc_heating_estimator[shell]
        # coll_deexc_heating = 0.0

        # coll_exc_heating_balance = self._calculate_coll_exc_cooling_rate_regemorter(
        #    lines, electron_densities[shell], t_electron, level_number_density[shell], level_boltzmann_factor_Te)

        # coll_deexc_heating_old, coll_exc_cooling_old = self._calculate_coll_exc_cooling_rate_regemorter(
        #    lines, electron_densities[shell], t_electron,
        #    level_number_density[shell], level_boltzmann_factor_Te)

        #        print 'Heating:', (coll_deexc_heating_old -
        #                coll_deexc_heating[shell])/coll_deexc_heating_old
        #
        #        print 'Cooling:', (coll_exc_cooling_old -
        #                coll_exc_cooling[shell])/coll_exc_cooling_old
        #        #self.coll_exc_heating_evolution[self._counter][shell] = coll_exc_heating_balance
        # self.coll_t_electrons[self._counter][shell] = t_electron

        adiabatic_cooling = self._calculate_adiabatic_cooling(
            electron_densities[shell], t_electron, time_explosion, shell
        )

        total_cooling = (
            fb_cooling + ff_cooling + coll_ion_cooling + coll_exc_cooling[shell]
        )
        total_heating = (
            bf_heating
            + ff_heating
            + coll_ion_heating
            + coll_deexc_heating[shell]
        )

        # if self.plasma_parent.niter > self.plasma_parent.niter_ly:
        #    total_heating += ly_cooling_rate

        # print "Frac Adiabatic Cooling:", adiabatic_cooling/total_cooling
        # total_cooling += adiabatic_cooling

        total_heating_rate = total_heating - total_cooling
        conv = (total_heating - total_cooling) / total_cooling

        if not np.isfinite(conv):
            print("Nonfinte Heating rate:")
            print(f"Shell: {shell}, H: {total_heating}, Q: {total_cooling}")
            print("Setting conv to finite value 1e-16")
            conv = 1e-16
        # total_heating_rate = bf_heating + ff_heating_balance - fb_cooling + coll_ion_heating_balance + \
        #                     coll_deexc_heating + coll_exc_heating_balance
        # TODO
        if mode == "normal":
            retval = total_heating_rate
        elif mode == "detailed":
            retval = (total_heating_rate, conv, sp_fb_cooling_rate)
        else:
            raise ValueError

        return retval

    def calculate_from_radiation_field_model(
        self, photo_ion_cross_sections, w, t_rad
    ):
        j_nus = self._calculate_j_nus(photo_ion_cross_sections, w, t_rad)
        nu_i = photo_ion_cross_sections["nu"].groupby(level=[0, 1, 2]).first()
        nu_is = nu_i.loc[photo_ion_cross_sections.index]
        photoion_heating_coeff = j_nus.multiply(
            4.0
            * np.pi
            * photo_ion_cross_sections["x_sect"]
            * (1 - nu_is / photo_ion_cross_sections["nu"]),
            axis=0,
        )
        photoion_heating_coeff.insert(0, "nu", photo_ion_cross_sections["nu"])
        photoion_heating_coeff = photoion_heating_coeff.groupby(level=[0, 1, 2])
        tmp = {}
        for i in range(len(w)):
            tmp[i] = photoion_heating_coeff.apply(
                lambda sub: simps(sub[i], sub["nu"], even="first")
            )
        photoion_heating_coeff = pd.DataFrame(tmp)

        return photoion_heating_coeff

    @staticmethod
    def _calculate_j_nus(photoionization_data, ws, t_rads):
        nus = photoionization_data["nu"].values
        j_nus = ws * intensity_black_body(nus[np.newaxis].T, t_rads)
        return pd.DataFrame(j_nus, index=photoionization_data.index)

    def _update_Lymann_continuum(
        self, bf_heat_coeff, photo_ion_cross_sections, t_electron, b
    ):
        W = b**-1
        x_sect_Ly = photo_ion_cross_sections.loc[(1, 0, 0)]
        bf_heat_Ly = self.calculate_from_radiation_field_model(
            x_sect_Ly, np.array([W]), np.array(t_electron)
        )
        bf_heat_coeff.loc[(1, 0, 0)] = bf_heat_Ly.values


class Logger(ProcessingPlasmaProperty):
    outputs = ("logging",)
    latex_name = ("",)
    # TODO
    latex_formula = ("",)

    def calculate(
        self, ion_number_density, level_number_density, t_electrons, b, gamma
    ):
        if not hasattr(self, "t_electrons"):
            self.t_electrons = []
            self.ion_number_density = []
            self.ion_frac = []
            self.level_number_density = []
            self.b = []
            self.b0 = []
            self.gamma_sp = []

        self.t_electrons.append(t_electrons)
        self.ion_number_density.append(ion_number_density)
        self.level_number_density.append(level_number_density)
        self.b.append(b)
        self.b0.append(b.loc[(1, 0, 0)])
        self.gamma_sp.append(gamma)
        ion_frac = ion_number_density.loc[(1, 0)] / ion_number_density.sum()
        self.ion_frac.append(ion_frac)
        logdata = {
            "ion_frac": self.ion_frac,
            "t_electrons": self.t_electrons,
            "gamma": gamma,
            "b0": self.b0,
        }
        return logdata
