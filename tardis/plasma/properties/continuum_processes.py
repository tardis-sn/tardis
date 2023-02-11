import logging

import numpy as np
import pandas as pd

from numba import prange, njit
from tardis import constants as const

from tardis.plasma.exceptions import PlasmaException
from tardis.plasma.properties.base import (
    ProcessingPlasmaProperty,
    Input,
    TransitionProbabilitiesProperty,
)
from tardis.plasma.properties.j_blues import JBluesDiluteBlackBody

__all__ = [
    "SpontRecombRateCoeff",
    "StimRecombRateCoeff",
    "PhotoIonRateCoeff",
    "PhotoIonEstimatorsNormFactor",
    "PhotoIonRateCoeffEstimator",
    "StimRecombRateCoeffEstimator",
    "CorrPhotoIonRateCoeff",
    "BfHeatingRateCoeffEstimator",
    "StimRecombCoolingRateCoeffEstimator",
    "SpontRecombCoolingRateCoeff",
    "RawRecombTransProbs",
    "RawPhotoIonTransProbs",
    "CollDeexcRateCoeff",
    "CollExcRateCoeff",
    "RawCollisionTransProbs",
    "AdiabaticCoolingRate",
    "FreeFreeCoolingRate",
    "FreeBoundCoolingRate",
    "BoundFreeOpacity",
    "LevelNumberDensityLTE",
    "PhotoIonBoltzmannFactor",
    "FreeBoundEmissionCDF",
    "RawTwoPhotonTransProbs",
    "TwoPhotonEmissionCDF",
    "TwoPhotonFrequencySampler",
    "CollIonRateCoeffSeaton",
    "CollRecombRateCoeff",
    "RawCollIonTransProbs",
    "BoundFreeOpacityInterpolator",
    "FreeFreeOpacity",
    "ContinuumOpacityCalculator",
    "FreeFreeFrequencySampler",
    "FreeBoundFrequencySampler",
]

N_A = const.N_A.cgs.value
K_B = const.k_B.cgs.value
C = const.c.cgs.value
H = const.h.cgs.value
A0 = const.a0.cgs.value
M_E = const.m_e.cgs.value
E = const.e.esu.value
BETA_COLL = (H**4 / (8 * K_B * M_E**3 * np.pi**3)) ** 0.5
F_K = (
    16
    / (3.0 * np.sqrt(3))
    * np.sqrt((2 * np.pi) ** 3 * K_B / (H**2 * M_E**3))
    * (E**2 / C) ** 3
)  # See Eq. 19 in Sutherland, R. S. 1998, MNRAS, 300, 321
FF_OPAC_CONST = (
    (2 * np.pi / (3 * M_E * K_B)) ** 0.5 * 4 * E**6 / (3 * M_E * H * C)
)  # See Eq. 6.1.8 in http://personal.psu.edu/rbc3/A534/lec6.pdf

logger = logging.getLogger(__name__)

njit_dict = {"fastmath": False, "parallel": False}


@njit(**njit_dict)
def integrate_array_by_blocks(f, x, block_references):
    """
    Integrate a function over blocks.

    This function integrates a function `f` defined at locations `x`
    over blocks given in `block_references`.

    Parameters
    ----------
    f : numpy.ndarray, dtype float
        2D input array to integrate.
    x : numpy.ndarray, dtype float
        1D array with the sample points corresponding to the `f` values.
    block_references : numpy.ndarray, dtype int
        1D array with the start indices of the blocks to be integrated.

    Returns
    -------
    numpy.ndarray, dtype float
        2D array with integrated values.
    """
    integrated = np.zeros((len(block_references) - 1, f.shape[1]))
    for i in prange(f.shape[1]):  # columns
        for j in prange(len(integrated)):  # rows
            start = block_references[j]
            stop = block_references[j + 1]
            integrated[j, i] = np.trapz(f[start:stop, i], x[start:stop])
    return integrated


# It is currently not possible to use scipy.integrate.cumulative_trapezoid in
# numba. So here is my own implementation.
@njit(**njit_dict)
def numba_cumulative_trapezoid(f, x):
    """
    Cumulatively integrate f(x) using the composite trapezoidal rule.

    Parameters
    ----------
    f : numpy.ndarray, dtype float
        Input array to integrate.
    x : numpy.ndarray, dtype float
        The coordinate to integrate along.

    Returns
    -------
    numpy.ndarray, dtype float
        The result of cumulative integration of f along x
    """
    integ = (np.diff(x) * (f[1:] + f[:-1]) / 2.0).cumsum()
    return integ / integ[-1]


@njit(**njit_dict)
def cumulative_integrate_array_by_blocks(f, x, block_references):
    """
    Cumulatively integrate a function over blocks.

    This function cumulatively integrates a function `f` defined at
    locations `x` over blocks given in `block_references`.

    Parameters
    ----------
    f : numpy.ndarray, dtype float
        Input array to integrate. Shape is (N_freq, N_shells), where
        N_freq is the number of frequency values and N_shells is the number
        of computational shells.
    x : numpy.ndarray, dtype float
        The sample points corresponding to the `f` values. Shape is (N_freq,).
    block_references : numpy.ndarray, dtype int
        The start indices of the blocks to be integrated. Shape is (N_blocks,).

    Returns
    -------
    numpy.ndarray, dtype float
        Array with cumulatively integrated values. Shape is (N_freq, N_shells)
        same as f.
    """
    n_rows = len(block_references) - 1
    integrated = np.zeros_like(f)
    for i in prange(f.shape[1]):  # columns
        # TODO: Avoid this loop through vectorization of cumulative_trapezoid
        for j in prange(n_rows):  # rows
            start = block_references[j]
            stop = block_references[j + 1]
            integrated[start + 1 : stop, i] = numba_cumulative_trapezoid(
                f[start:stop, i], x[start:stop]
            )
    return integrated


def get_ion_multi_index(multi_index_full, next_higher=True):
    """
    Calculate the corresponding ion MultiIndex for a level MultiIndex.

    Parameters
    ----------
    multi_index_full : pandas.MultiIndex (atomic_number, ion_number,
                                          level_number)
    next_higher : bool, default True
        If True use ion number of next higher ion, else use ion_number from
        multi_index_full.

    Returns
    -------
    pandas.MultiIndex (atomic_number, ion_number)
       Ion MultiIndex for the given level MultiIndex.
    """
    atomic_number = multi_index_full.get_level_values(0)
    ion_number = multi_index_full.get_level_values(1)
    if next_higher is True:
        ion_number += 1
    return pd.MultiIndex.from_arrays([atomic_number, ion_number])


def get_ground_state_multi_index(multi_index_full):
    """
    Calculate the ground-state MultiIndex for the next higher ion.

    Parameters
    ----------
    multi_index_full : pandas.MultiIndex (atomic_number, ion_number,
                                          level_number)

    Returns
    -------
    pandas.MultiIndex (atomic_number, ion_number)
        Ground-state MultiIndex for the next higher ion.
    """
    atomic_number = multi_index_full.get_level_values(0)
    ion_number = multi_index_full.get_level_values(1) + 1
    level_number = np.zeros_like(ion_number)
    return pd.MultiIndex.from_arrays([atomic_number, ion_number, level_number])


def cooling_rate_series2dataframe(cooling_rate_series, destination_level_idx):
    """
    Transform cooling-rate Series to DataFrame.

    This function transforms a Series with cooling rates into
    an indexed DataFrame that can be used in MarkovChainTransProbs.

    Parameters
    ----------
    cooling_rate_series : pandas.Series, dtype float
        Cooling rates for a process with a single destination idx.
        Examples are adiabatic cooling or free-free cooling.
    destination_level_idx : str
        Destination idx of the cooling process; for example
        'adiabatic' for adiabatic cooling.

    Returns
    -------
    cooling_rate_frame : pandas.DataFrame, dtype float
        Indexed by source_level_idx, destination_level_idx, transition_type
        for the use in MarkovChainTransProbs.
    """
    index_names = [
        "source_level_idx",
        "destination_level_idx",
        "transition_type",
    ]
    index = pd.MultiIndex.from_tuples(
        [("k", destination_level_idx, -1)], names=index_names
    )
    cooling_rate_frame = pd.DataFrame(
        cooling_rate_series.values[np.newaxis], index=index
    )
    return cooling_rate_frame


def bf_estimator_array2frame(bf_estimator_array, level2continuum_idx):
    """
    Transform a bound-free estimator array to a DataFrame.

    This function transforms a bound-free estimator array with entries
    sorted by frequency to a multi-indexed DataFrame sorted by level.

    Parameters
    ----------
    bf_estimator_array : numpy.ndarray, dtype float
        Array of bound-free estimators (e.g., for the stimulated recombination rate)
        with entries sorted by the threshold frequency of the bound-free continuum.
    level2continuum_idx : pandas.Series, dtype int
        Maps a level MultiIndex (atomic_number, ion_number, level_number) to
        the continuum_idx of the corresponding bound-free continuum (which are
        sorted by decreasing frequency).

    Returns
    -------
    pandas.DataFrame, dtype float
        Bound-free estimators indexed by (atomic_number, ion_number, level_number).
    """
    bf_estimator_frame = pd.DataFrame(
        bf_estimator_array, index=level2continuum_idx.index
    ).sort_index()
    bf_estimator_frame.columns.name = "Shell No."
    return bf_estimator_frame


class IndexSetterMixin(object):
    @staticmethod
    def set_index(p, photo_ion_idx, transition_type=0, reverse=True):
        idx = photo_ion_idx.loc[p.index]
        transition_type = transition_type * np.ones_like(
            idx.destination_level_idx
        )
        transition_type = pd.Series(transition_type, name="transition_type")
        idx_arrays = [idx.source_level_idx, idx.destination_level_idx]
        if reverse:
            idx_arrays = idx_arrays[::-1]
        idx_arrays.append(transition_type)
        index = pd.MultiIndex.from_arrays(idx_arrays)
        if reverse:
            index.names = index.names[:-1][::-1] + [index.names[-1]]
        p = p.set_index(index, drop=True)
        return p


class SpontRecombRateCoeff(ProcessingPlasmaProperty):
    """
    Attributes
    ----------
    alpha_sp : pandas.DataFrame, dtype float
        The rate coefficient for spontaneous recombination.
    """

    outputs = ("alpha_sp",)
    latex_name = (r"\alpha^{\textrm{sp}}",)

    def calculate(
        self,
        photo_ion_cross_sections,
        t_electrons,
        photo_ion_block_references,
        photo_ion_index,
        phi_ik,
        boltzmann_factor_photo_ion,
    ):
        x_sect = photo_ion_cross_sections["x_sect"].values
        nu = photo_ion_cross_sections["nu"].values

        alpha_sp = 8 * np.pi * x_sect * nu**2 / C**2
        alpha_sp = alpha_sp[:, np.newaxis]
        alpha_sp = alpha_sp * boltzmann_factor_photo_ion
        alpha_sp = integrate_array_by_blocks(
            alpha_sp, nu, photo_ion_block_references
        )
        alpha_sp = pd.DataFrame(alpha_sp, index=photo_ion_index)
        return alpha_sp * phi_ik.loc[alpha_sp.index]


class SpontRecombCoolingRateCoeff(ProcessingPlasmaProperty):
    """
    Attributes
    ----------
    c_fb_sp : pandas.DataFrame, dtype float
        The rate coefficient for cooling by
        spontaneous recombination.
    """

    outputs = ("c_fb_sp",)
    latex_name = (r"c^{\textrm{sp}}_{\textrm{fb}}",)

    def calculate(
        self,
        photo_ion_cross_sections,
        t_electrons,
        photo_ion_block_references,
        photo_ion_index,
        phi_ik,
        nu_i,
        boltzmann_factor_photo_ion,
    ):
        x_sect = photo_ion_cross_sections["x_sect"].values
        nu = photo_ion_cross_sections["nu"].values
        factor = (1 - nu_i / photo_ion_cross_sections["nu"]).values
        alpha_sp = (8 * np.pi * x_sect * factor * nu**3 / C**2) * H
        alpha_sp = alpha_sp[:, np.newaxis]
        alpha_sp = alpha_sp * boltzmann_factor_photo_ion
        alpha_sp = integrate_array_by_blocks(
            alpha_sp, nu, photo_ion_block_references
        )
        alpha_sp = pd.DataFrame(alpha_sp, index=photo_ion_index)
        return alpha_sp * phi_ik.loc[alpha_sp.index]


class FreeBoundEmissionCDF(ProcessingPlasmaProperty):
    """
    Attributes
    ----------
    fb_emission_cdf : pandas.DataFrame, dtype float
        The cumulative distribution function (CDF) for the frequencies of
        energy packets emitted in free-bound transitions. The tabulated CDF
        is used to sample packet frequencies in the Monte Carlo simulation.
        We use the same CDF for free-bound emission from k- and i-packets
        (in contrast to ARTIS).
    """

    outputs = ("fb_emission_cdf",)
    latex_name = (r"P(\nu_{bf, emission}) \leq \nu)",)

    def calculate(
        self,
        photo_ion_cross_sections,
        t_electrons,
        photo_ion_block_references,
        photo_ion_index,
        nu_i,
        boltzmann_factor_photo_ion,
    ):
        x_sect = photo_ion_cross_sections["x_sect"].values
        nu = photo_ion_cross_sections["nu"].values
        # alpha_sp_E will be missing a lot of prefactors since we are only
        # interested in relative values here
        alpha_sp_E = nu**3 * x_sect
        alpha_sp_E = alpha_sp_E[:, np.newaxis]
        alpha_sp_E = alpha_sp_E * boltzmann_factor_photo_ion
        alpha_sp_E = cumulative_integrate_array_by_blocks(
            alpha_sp_E, nu, photo_ion_block_references
        )
        fb_emission_cdf = pd.DataFrame(
            alpha_sp_E, index=photo_ion_cross_sections.index
        )
        return fb_emission_cdf


class PhotoIonRateCoeff(ProcessingPlasmaProperty):
    """
    Attributes
    ----------
    gamma : pandas.DataFrame, dtype float
        The rate coefficient for radiative ionization.
    """

    outputs = ("gamma",)
    latex_name = (r"\gamma",)

    def calculate(
        self,
        photo_ion_cross_sections,
        gamma_estimator,
        photo_ion_norm_factor,
        photo_ion_block_references,
        photo_ion_index,
        t_rad,
        w,
        level2continuum_idx,
    ):
        # Used for initialization
        if gamma_estimator is None:
            gamma = self.calculate_from_dilute_bb(
                photo_ion_cross_sections,
                photo_ion_block_references,
                photo_ion_index,
                t_rad,
                w,
            )
        else:
            gamma_estimator = bf_estimator_array2frame(
                gamma_estimator, level2continuum_idx
            )
            gamma = gamma_estimator * photo_ion_norm_factor.value

        return gamma

    @staticmethod
    def calculate_from_dilute_bb(
        photo_ion_cross_sections,
        photo_ion_block_references,
        photo_ion_index,
        t_rad,
        w,
    ):
        nu = photo_ion_cross_sections["nu"]
        x_sect = photo_ion_cross_sections["x_sect"]
        j_nus = JBluesDiluteBlackBody.calculate(
            photo_ion_cross_sections, nu, t_rad, w
        )
        gamma = j_nus.multiply(4.0 * np.pi * x_sect / nu / H, axis=0)
        gamma = integrate_array_by_blocks(
            gamma.values, nu.values, photo_ion_block_references
        )
        gamma = pd.DataFrame(gamma, index=photo_ion_index)
        return gamma


class StimRecombRateCoeff(ProcessingPlasmaProperty):
    """
    Attributes
    ----------
    alpha_stim : pandas.DataFrame, dtype float
        The rate coefficient for stimulated recombination.
    """

    outputs = ("alpha_stim",)
    latex_name = (r"\alpha^{\textrm{stim}}",)

    def calculate(
        self,
        photo_ion_cross_sections,
        alpha_stim_estimator,
        photo_ion_norm_factor,
        photo_ion_block_references,
        photo_ion_index,
        t_rad,
        w,
        phi_ik,
        t_electrons,
        boltzmann_factor_photo_ion,
        level2continuum_idx,
    ):
        # Used for initialization
        if alpha_stim_estimator is None:
            alpha_stim = self.calculate_from_dilute_bb(
                photo_ion_cross_sections,
                photo_ion_block_references,
                photo_ion_index,
                t_rad,
                w,
                t_electrons,
                boltzmann_factor_photo_ion,
            )
        else:
            alpha_stim_estimator = bf_estimator_array2frame(
                alpha_stim_estimator, level2continuum_idx
            )
            alpha_stim = alpha_stim_estimator * photo_ion_norm_factor
        alpha_stim *= phi_ik.loc[alpha_stim.index]
        return alpha_stim

    @staticmethod
    def calculate_from_dilute_bb(
        photo_ion_cross_sections,
        photo_ion_block_references,
        photo_ion_index,
        t_rad,
        w,
        t_electrons,
        boltzmann_factor_photo_ion,
    ):
        nu = photo_ion_cross_sections["nu"]
        x_sect = photo_ion_cross_sections["x_sect"]
        j_nus = JBluesDiluteBlackBody.calculate(
            photo_ion_cross_sections, nu, t_rad, w
        )
        j_nus *= boltzmann_factor_photo_ion
        alpha_stim = j_nus.multiply(4.0 * np.pi * x_sect / nu / H, axis=0)
        alpha_stim = integrate_array_by_blocks(
            alpha_stim.values, nu.values, photo_ion_block_references
        )
        alpha_stim = pd.DataFrame(alpha_stim, index=photo_ion_index)
        return alpha_stim


class RawRecombTransProbs(TransitionProbabilitiesProperty, IndexSetterMixin):
    """
    Attributes
    ----------
    p_recomb : pandas.DataFrame, dtype float
        The unnormalized transition probabilities for
        spontaneous recombination.
    """

    outputs = ("p_recomb",)
    transition_probabilities_outputs = ("p_recomb",)
    latex_name = (r"p^{\textrm{recomb}}",)

    def calculate(self, alpha_sp, nu_i, energy_i, photo_ion_idx):
        p_recomb_deactivation = alpha_sp.multiply(nu_i, axis=0) * H
        p_recomb_deactivation = self.set_index(
            p_recomb_deactivation, photo_ion_idx, transition_type=-1
        )

        p_recomb_internal = alpha_sp.multiply(energy_i, axis=0)
        p_recomb_internal = self.set_index(
            p_recomb_internal, photo_ion_idx, transition_type=0
        )
        p_recomb = pd.concat([p_recomb_deactivation, p_recomb_internal])
        return p_recomb


class RawPhotoIonTransProbs(TransitionProbabilitiesProperty, IndexSetterMixin):
    """
    Attributes
    ----------
    p_photo_ion : pandas.DataFrame, dtype float
        The unnormalized transition probabilities for
        radiative ionization.
    """

    outputs = ("p_photo_ion",)
    transition_probabilities_outputs = ("p_photo_ion",)
    latex_name = (r"p^{\textrm{photo_ion}}",)

    def calculate(self, gamma_corr, energy_i, photo_ion_idx):
        p_photo_ion = gamma_corr.multiply(energy_i, axis=0)
        p_photo_ion = self.set_index(p_photo_ion, photo_ion_idx, reverse=False)
        return p_photo_ion


class CorrPhotoIonRateCoeff(ProcessingPlasmaProperty):
    """
    Attributes
    ----------
    gamma_corr : pandas.DataFrame, dtype float
        The rate coefficient for radiative ionization corrected for
        stimulated recombination.
    """

    outputs = ("gamma_corr",)
    latex_name = (r"\gamma_\mathrm{corr}",)

    def calculate(
        self,
        gamma,
        alpha_stim,
        electron_densities,
        ion_number_density,
        level_number_density,
    ):
        n_k_index = get_ion_multi_index(alpha_stim.index)
        n_k = ion_number_density.loc[n_k_index].values
        n_i = level_number_density.loc[alpha_stim.index].values
        gamma_corr = gamma - (alpha_stim * n_k / n_i).multiply(
            electron_densities
        )
        num_neg_elements = (gamma_corr < 0).sum().sum()
        if num_neg_elements:
            raise PlasmaException(
                "Negative values in CorrPhotoIonRateCoeff.  Try raising the number of montecarlo packets."
            )
        return gamma_corr


class PhotoIonEstimatorsNormFactor(ProcessingPlasmaProperty):
    outputs = ("photo_ion_norm_factor",)
    latex_name = (r"\frac{1}{t_\textrm{simulation volume h}}",)

    @staticmethod
    def calculate(time_simulation, volume):
        return (time_simulation * volume * H) ** -1


class PhotoIonRateCoeffEstimator(Input):
    """
    Attributes
    ----------
    gamma_estimator : pandas.DataFrame, dtype float
        Unnormalized MC estimator for the rate coefficient for radiative
        ionization.
    """

    outputs = ("gamma_estimator",)
    latex_name = (r"\gamma_\textrm{estim}",)


class StimRecombRateCoeffEstimator(Input):
    """
    Attributes
    ----------
    alpha_stim_estimator : pandas.DataFrame, dtype float
        Unnormalized MC estimator for the rate coefficient for stimulated
        recombination.
    """

    outputs = ("alpha_stim_estimator",)
    latex_name = (r"\alpha^{\textrm{stim}}_\textrm{estim}",)


class StimRecombCoolingRateCoeffEstimator(Input):
    """
    Attributes
    ----------
    stim_recomb_cooling_coeff_estimator : pandas.DataFrame, dtype float
        Unnormalized MC estimator for the stimulated recombination cooling rate
        coefficient.
    """

    outputs = ("stim_recomb_cooling_coeff_estimator",)


class BfHeatingRateCoeffEstimator(Input):
    """
    Attributes
    ----------
    bf_heating_coeff_estimator : pandas.DataFrame, dtype float
        Unnormalized MC estimator for the rate
        coefficient for bound-free heating.
    """

    outputs = ("bf_heating_coeff_estimator",)
    latex_name = (r"h_\textrm{bf}_\textrm{estim}",)


class CollExcRateCoeff(ProcessingPlasmaProperty):
    """
    Attributes
    ----------
    coll_exc_coeff : pandas.DataFrame, dtype float
        Rate coefficient for collisional excitation.
    """

    outputs = ("coll_exc_coeff",)
    latex_name = ("c_{lu}",)

    def calculate(self, yg_interp, yg_index, t_electrons, delta_E_yg):
        yg = yg_interp(t_electrons)
        boltzmann_factor = np.exp(
            -delta_E_yg.values[np.newaxis].T / (t_electrons * K_B)
        )
        q_ij = (
            BETA_COLL / np.sqrt(t_electrons) * yg * boltzmann_factor
        )  # see formula A2 in Przybilla, Butler 2004 - Apj 609, 1181
        return pd.DataFrame(q_ij, index=yg_index)


class CollDeexcRateCoeff(ProcessingPlasmaProperty):
    """
    Attributes
    ----------
    coll_deexc_coeff : pandas.DataFrame, dtype float
        Rate coefficient for collisional deexcitation.
    """

    outputs = ("coll_deexc_coeff",)
    latex_name = ("c_{ul}",)

    def calculate(self, thermal_lte_level_boltzmann_factor, coll_exc_coeff):
        level_lower_index = coll_exc_coeff.index.droplevel("level_number_upper")
        level_upper_index = coll_exc_coeff.index.droplevel("level_number_lower")

        n_lower_prop = thermal_lte_level_boltzmann_factor.loc[
            level_lower_index
        ].values
        n_upper_prop = thermal_lte_level_boltzmann_factor.loc[
            level_upper_index
        ].values

        coll_deexc_coeff = coll_exc_coeff * n_lower_prop / n_upper_prop
        return coll_deexc_coeff


class RawCollisionTransProbs(TransitionProbabilitiesProperty, IndexSetterMixin):
    """
    Attributes
    ----------
    p_coll : pandas.DataFrame, dtype float
        The unnormalized transition probabilities for
        collisional excitation.
    """

    outputs = ("p_coll",)
    transition_probabilities_outputs = ("p_coll",)
    latex_name = (r"p^{\textrm{coll}}",)

    def calculate(
        self,
        coll_exc_coeff,
        coll_deexc_coeff,
        yg_idx,
        electron_densities,
        delta_E_yg,
        atomic_data,
        level_number_density,
    ):
        p_deexc_deactivation = (coll_deexc_coeff * electron_densities).multiply(
            delta_E_yg.values, axis=0
        )
        p_deexc_deactivation = self.set_index(
            p_deexc_deactivation, yg_idx, reverse=True
        )
        p_deexc_deactivation = p_deexc_deactivation.groupby(level=[0]).sum()
        index_dd = pd.MultiIndex.from_product(
            [p_deexc_deactivation.index.values, ["k"], [0]],
            names=list(yg_idx.columns) + ["transition_type"],
        )
        p_deexc_deactivation = p_deexc_deactivation.set_index(index_dd)

        level_lower_index = coll_deexc_coeff.index.droplevel(
            "level_number_upper"
        )
        energy_lower = atomic_data.levels.energy.loc[level_lower_index]
        p_deexc_internal = (coll_deexc_coeff * electron_densities).multiply(
            energy_lower.values, axis=0
        )
        p_deexc_internal = self.set_index(
            p_deexc_internal, yg_idx, transition_type=0, reverse=True
        )

        p_exc_internal = (coll_exc_coeff * electron_densities).multiply(
            energy_lower.values, axis=0
        )
        p_exc_internal = self.set_index(
            p_exc_internal, yg_idx, transition_type=0, reverse=False
        )
        p_exc_cool = (coll_exc_coeff * electron_densities).multiply(
            delta_E_yg.values, axis=0
        )
        p_exc_cool = (
            p_exc_cool * level_number_density.loc[level_lower_index].values
        )
        p_exc_cool = self.set_index(p_exc_cool, yg_idx, reverse=False)
        p_exc_cool = p_exc_cool.groupby(level="destination_level_idx").sum()
        exc_cool_index = pd.MultiIndex.from_product(
            [["k"], p_exc_cool.index.values, [0]],
            names=list(yg_idx.columns) + ["transition_type"],
        )
        p_exc_cool = p_exc_cool.set_index(exc_cool_index)
        p_coll = pd.concat(
            [p_deexc_deactivation, p_deexc_internal, p_exc_internal, p_exc_cool]
        )
        return p_coll


class RawTwoPhotonTransProbs(TransitionProbabilitiesProperty, IndexSetterMixin):
    """
    Attributes
    ----------
    p_two_photon : pandas.DataFrame, dtype float
        The unnormalized transition probabilities for two photon decay.
    """

    outputs = ("p_two_photon",)
    transition_probabilities_outputs = ("p_two_photon",)

    def calculate(self, two_photon_data, two_photon_idx, density):
        no_shells = len(density)
        p_two_phot = two_photon_data.A_ul * two_photon_data.nu0 * H
        p_two_phot = pd.concat([p_two_phot] * no_shells, axis=1)
        # TODO: In principle there could be internal two photon transitions
        p_two_phot = self.set_index(
            p_two_phot,
            two_photon_idx,
            transition_type=-1,
            reverse=False,
        )
        p_two_phot.index = p_two_phot.index.set_levels(
            ["two-photon"], level="destination_level_idx"
        )
        return p_two_phot


class TwoPhotonEmissionCDF(ProcessingPlasmaProperty):
    """
    Attributes
    ----------
    two_photon_emission_cdf : pandas.DataFrame, dtype float
        The cumulative distribution function (CDF) for the frequencies of
        energy packets emitted in two photon transitions. The tabulated CDF
        is used to sample packet frequencies in the Monte Carlo simulation.
    """

    outputs = ("two_photon_emission_cdf",)

    def calculate(self, two_photon_data):
        bins = 500
        # The number of two photon transitions is very small
        # and the CDF has to be calculated only once.
        # There is no need to vectorize the calculation.
        emission_cdfs = []
        for index, row in two_photon_data.iterrows():
            alpha = row.alpha
            beta = row.beta
            gamma = row.gamma
            nu = np.linspace(0.0, row.nu0, bins)
            y = nu / row.nu0
            j_nu = self.calculate_j_nu(y, alpha, beta, gamma)

            cdf = np.zeros_like(nu)
            cdf[1:] = numba_cumulative_trapezoid(j_nu, nu)
            cdf /= cdf[-1]
            index_cdf = pd.MultiIndex.from_tuples([index] * bins)
            cdf = pd.DataFrame({"nu": nu, "cdf": cdf}, index=index_cdf)
            emission_cdfs.append(cdf)
        return pd.concat(emission_cdfs)

    @staticmethod
    def calculate_j_nu(y, alpha, beta, gamma):
        """
        Calculate two photon emissivity.

        This function calculates the two photon emissivity in the frequency
        scale based on Eq. 2 and Eq. 3 in Nussbaumer & Schmutz (1984). The
        emissivity is not normalized since it is only used to calculate
        relative emission probabilities.

        Parameters
        ----------
        y : numpy.ndarray, dtype float
            Emission frequency divided by that of the normal line
            transition corresponding to the two photon decay.
        alpha : float
            Fit coefficient.
        beta : float
            Fit coefficient.
        gamma : float
            Fit coefficient.

        Returns
        -------
        numpy.ndarray, dtype float
            Unnormalized two photon emissivity in the frequency scale.
        """
        ay = y * (1 - y) * (1 - (4 * y * (1 - y)) ** gamma)
        ay += alpha * (y * (1 - y)) ** beta * (4 * y * (1 - y)) ** gamma
        j_nu = ay * y
        return j_nu


class AdiabaticCoolingRate(TransitionProbabilitiesProperty):
    """
    Attributes
    ----------
    cool_rate_adiabatic : pandas.DataFrame, dtype float
        The adiabatic cooling rate of the electron gas.
    """

    outputs = ("cool_rate_adiabatic",)
    transition_probabilities_outputs = ("cool_rate_adiabatic",)
    latex_name = (r"C_{\textrm{adiabatic}}",)

    def calculate(self, electron_densities, t_electrons, time_explosion):
        cool_rate_adiabatic = (
            3.0 * electron_densities * K_B * t_electrons
        ) / time_explosion

        cool_rate_adiabatic = cooling_rate_series2dataframe(
            cool_rate_adiabatic, destination_level_idx="adiabatic"
        )
        return cool_rate_adiabatic


class FreeFreeCoolingRate(TransitionProbabilitiesProperty):
    """
    Attributes
    ----------
    cool_rate_ff : pandas.DataFrame, dtype float
        The free-free cooling rate of the electron gas.
    ff_cooling_factor : pandas.Series, dtype float
        Pre-factor needed in the calculation of the free-free cooling rate and
        the free-free opacity.

    Notes
    -----
    This implementation uses a free-free Gaunt factor of one for all species
    and ionization stages, which is an approximation.
    """

    outputs = ("cool_rate_ff", "ff_cooling_factor")
    transition_probabilities_outputs = ("cool_rate_ff",)
    latex_name = (r"C^{\textrm{ff}}",)

    def calculate(self, ion_number_density, electron_densities, t_electrons):
        ff_cooling_factor = self._calculate_ff_cooling_factor(
            ion_number_density, electron_densities
        )
        cool_rate_ff = F_K * np.sqrt(t_electrons) * ff_cooling_factor
        cool_rate_ff = cooling_rate_series2dataframe(
            cool_rate_ff, destination_level_idx="ff"
        )
        return cool_rate_ff, ff_cooling_factor.values

    @staticmethod
    def _calculate_ff_cooling_factor(ion_number_density, electron_densities):
        ion_charge = ion_number_density.index.get_level_values(1).values
        factor = (
            electron_densities
            * ion_number_density.multiply(ion_charge**2, axis=0).sum()
        )
        return factor


class FreeFreeOpacity(ProcessingPlasmaProperty):
    """
    Attributes
    ----------
    cool_rate_ff : pandas.DataFrame, dtype float
        The free-free cooling rate of the electron gas.
    ff_cooling_factor : pandas.Series, dtype float
        Pre-factor needed in the calculation of the free-free opacity.
    """

    outputs = ("chi_ff_calculator",)

    def calculate(self, t_electrons, ff_cooling_factor, electron_densities):
        ff_opacity_factor = ff_cooling_factor / np.sqrt(t_electrons)

        @njit(error_model="numpy", fastmath=True)
        def chi_ff(nu, shell):
            chi_ff = (
                FF_OPAC_CONST
                * ff_opacity_factor[shell]
                / nu**3
                * (1 - np.exp(-H * nu / (K_B * t_electrons[shell])))
            )
            return chi_ff

        return chi_ff


class FreeFreeFrequencySampler(ProcessingPlasmaProperty):
    """
    Attributes
    ----------
    nu_ff_sampler : float
        Frequency of the free-free emission process

    """

    outputs = ("nu_ff_sampler",)

    def calculate(self, t_electrons):
        @njit(error_model="numpy", fastmath=True)
        def nu_ff(shell):

            temperature = t_electrons[shell]
            zrand = np.random.random()
            return -K_B * temperature / H * np.log(zrand)

        return nu_ff


class FreeBoundFrequencySampler(ProcessingPlasmaProperty):
    """
    Attributes
    ----------
    nu_fb_sampler : float
        Frequency of the free-bounds emission process
    """

    outputs = ("nu_fb_sampler",)

    def calculate(
        self, photo_ion_cross_sections, fb_emission_cdf, level2continuum_idx
    ):

        phot_nus = photo_ion_cross_sections.nu.loc[level2continuum_idx.index]
        photo_ion_block_references = np.pad(
            phot_nus.groupby(level=[0, 1, 2], sort=False)
            .count()
            .values.cumsum(),
            [1, 0],
        )
        phot_nus = phot_nus.values
        emissivities = fb_emission_cdf.loc[level2continuum_idx.index].values

        @njit(error_model="numpy", fastmath=True)
        def nu_fb(shell, continuum_id):

            start = photo_ion_block_references[continuum_id]
            end = photo_ion_block_references[continuum_id + 1]
            phot_nus_block = phot_nus[start:end]
            em = emissivities[start:end, shell]

            zrand = np.random.random()
            idx = np.searchsorted(em, zrand, side="right")

            return phot_nus_block[idx] - (em[idx] - zrand) / (
                em[idx] - em[idx - 1]
            ) * (phot_nus_block[idx] - phot_nus_block[idx - 1])

        return nu_fb


class TwoPhotonFrequencySampler(ProcessingPlasmaProperty):
    """
    Attributes
    ----------
    nu_two_photon_sampler : float
        Frequency of the two-photon emission process
    """

    outputs = ("nu_two_photon_sampler",)

    def calculate(self, two_photon_emission_cdf):

        nus = two_photon_emission_cdf["nu"].values
        em = two_photon_emission_cdf["cdf"].values

        @njit(error_model="numpy", fastmath=True)
        def nu_two_photon():
            zrand = np.random.random()
            idx = np.searchsorted(em, zrand, side="right")

            return nus[idx] - (em[idx] - zrand) / (em[idx] - em[idx - 1]) * (
                nus[idx] - nus[idx - 1]
            )

        return nu_two_photon


class FreeBoundCoolingRate(TransitionProbabilitiesProperty):
    """
    Attributes
    ----------
    cool_rate_fb_total : pandas.DataFrame, dtype float
        The total free-bound cooling rate of the electron gas.
    cool_rate_fb : pandas.DataFrame, dtype float
        The individual free-bound cooling rates of the electron gas.
    p_fb_deactivation: pandas.DataFrame, dtype float
        Probabilities of free-bound cooling in a specific continuum
        (identified by its continuum_idx).
    """

    outputs = ("cool_rate_fb_tot", "cool_rate_fb", "p_fb_deactivation")
    transition_probabilities_outputs = ("cool_rate_fb_tot",)
    latex_name = (r"C^{\textrm{fb, tot}}", r"C^{\textrm{fb}}")

    def calculate(
        self,
        c_fb_sp,
        electron_densities,
        ion_number_density,
        level2continuum_idx,
    ):
        next_ion_stage_index = get_ion_multi_index(c_fb_sp.index)
        n_k = ion_number_density.loc[next_ion_stage_index]

        cool_rate_fb = c_fb_sp.multiply(electron_densities, axis=1) * n_k.values
        cool_rate_fb_tot = cooling_rate_series2dataframe(
            cool_rate_fb.sum(axis=0), "bf"
        )

        p_fb_deactivation = cool_rate_fb / cool_rate_fb_tot.values
        # TODO: this will be needed more often; put it in a function
        continuum_idx = level2continuum_idx.loc[p_fb_deactivation.index].values
        p_fb_deactivation = p_fb_deactivation.set_index(
            continuum_idx
        ).sort_index(ascending=True)
        p_fb_deactivation.index.name = "continuum_idx"
        return cool_rate_fb_tot, cool_rate_fb, p_fb_deactivation


class BoundFreeOpacity(ProcessingPlasmaProperty):
    """
    Attributes
    ----------
    chi_bf : pandas.DataFrame, dtype float
        Bound-free opacity corrected for stimulated emission.
    """

    outputs = ("chi_bf",)
    latex_name = (r"\chi^{\textrm{bf}}",)

    def calculate(
        self,
        photo_ion_cross_sections,
        t_electrons,
        phi_ik,
        level_number_density,
        lte_level_number_density,
        boltzmann_factor_photo_ion,
    ):
        x_sect = photo_ion_cross_sections["x_sect"].values

        n_i = level_number_density.loc[photo_ion_cross_sections.index]
        lte_n_i = lte_level_number_density.loc[photo_ion_cross_sections.index]
        chi_bf = (n_i - lte_n_i * boltzmann_factor_photo_ion).multiply(
            x_sect, axis=0
        )

        num_neg_elements = (chi_bf < 0).sum().sum()
        if num_neg_elements:
            raise PlasmaException("Negative values in bound-free opacity.")
        return chi_bf


class BoundFreeOpacityInterpolator(ProcessingPlasmaProperty):
    outputs = ("chi_bf_interpolator",)

    def calculate(
        self,
        chi_bf,
        photo_ion_cross_sections,
        get_current_bound_free_continua,
        level2continuum_idx,
    ):
        # Sort everything by descending frequeny (i.e. continuum_id)
        # TODO: Do we really gain much by sorting by continuum_id?
        phot_nus = photo_ion_cross_sections.nu.loc[level2continuum_idx.index]
        photo_ion_block_references = np.pad(
            phot_nus.groupby(level=[0, 1, 2], sort=False)
            .count()
            .values.cumsum(),
            [1, 0],
        )
        phot_nus = phot_nus.values
        chi_bf = chi_bf.loc[level2continuum_idx.index].values
        x_sect = photo_ion_cross_sections.x_sect.loc[
            level2continuum_idx.index
        ].values

        @njit(error_model="numpy", fastmath=True)
        def chi_bf_interpolator(nu, shell):
            """
            Interpolate the bound-free opacity.

            This function interpolates the tabulated bound-free opacities
            and cross-sections to new frequency values `nu`.

            Parameters
            ----------
            nu : float, dtype float
                Comoving frequency of the r-packet.
            shell : int, dtype float
                Current computational shell.

            Returns
            -------
            chi_bf_tot : float
                Total bound-free opacity at frequency `nu`.
            chi_bf_contributions : numpy.ndarray, dtype float
                Cumulative distribution function of the contributions of the
                individual bound free continua to the total bound-free opacity.
            current_continua : numpy.ndarray, dtype int
                Continuum ids for which absorption is possible for frequency `nu`.
            x_sect_bfs : numpy.ndarray, dtype float
                Photoionization cross-sections of all bound-free continua for
                which absorption is possible for frequency `nu`.
            """

            current_continua = get_current_bound_free_continua(nu)
            chi_bfs = np.zeros(len(current_continua))
            x_sect_bfs = np.zeros(len(current_continua))
            for i, continuum_id in enumerate(current_continua):
                start = photo_ion_block_references[continuum_id]
                end = photo_ion_block_references[continuum_id + 1]
                phot_nus_continuum = phot_nus[start:end]
                nu_idx = np.searchsorted(phot_nus_continuum, nu)
                interval = (
                    phot_nus_continuum[nu_idx] - phot_nus_continuum[nu_idx - 1]
                )
                high_weight = nu - phot_nus_continuum[nu_idx - 1]
                low_weight = phot_nus_continuum[nu_idx] - nu
                chi_bfs_continuum = chi_bf[start:end, shell]
                chi_bfs[i] = (
                    chi_bfs_continuum[nu_idx] * high_weight
                    + chi_bfs_continuum[nu_idx - 1] * low_weight
                ) / interval
                x_sect_bfs_continuum = x_sect[start:end]
                x_sect_bfs[i] = (
                    x_sect_bfs_continuum[nu_idx] * high_weight
                    + x_sect_bfs_continuum[nu_idx - 1] * low_weight
                ) / interval
            chi_bf_contributions = chi_bfs.cumsum()

            # If we are outside the range of frequencies
            # for which we have photo-ionization cross sections
            # we will have no local continuua and therefore
            # no bound-free interactions can occur
            # so we set the bound free opacity to zero
            if len(current_continua) == 0:
                chi_bf_tot = 0.0
            else:
                chi_bf_tot = chi_bf_contributions[-1]
                chi_bf_contributions /= chi_bf_tot

            return (
                chi_bf_tot,
                chi_bf_contributions,
                current_continua,
                x_sect_bfs,
            )

        return chi_bf_interpolator


class ContinuumOpacityCalculator(ProcessingPlasmaProperty):
    outputs = ("chi_continuum_calculator",)

    def calculate(self, chi_ff_calculator, chi_bf_interpolator):
        @njit(error_model="numpy", fastmath=True)
        def chi_continuum_calculator(nu, shell):
            (
                chi_bf_tot,
                chi_bf_contributions,
                current_continua,
                x_sect_bfs,
            ) = chi_bf_interpolator(nu, shell)
            chi_ff = chi_ff_calculator(nu, shell)
            return (
                chi_bf_tot,
                chi_bf_contributions,
                current_continua,
                x_sect_bfs,
                chi_ff,
            )

        return chi_continuum_calculator


class LevelNumberDensityLTE(ProcessingPlasmaProperty):
    """
    Attributes
    ----------
    lte_level_number_density : pandas.DataFrame, dtype float
    """

    outputs = ("lte_level_number_density",)
    latex_name = (r"n_{\textrm{i}}^*",)

    # TODO: only do this for continuum species
    def calculate(self, electron_densities, phi_ik, ion_number_density):
        next_higher_ion_index = get_ion_multi_index(
            phi_ik.index, next_higher=True
        )
        # TODO: Check that n_k is correct (and not n_k*)
        lte_level_number_density = (
            phi_ik * ion_number_density.loc[next_higher_ion_index].values
        ).multiply(electron_densities, axis=1)
        return lte_level_number_density


class PhotoIonBoltzmannFactor(ProcessingPlasmaProperty):
    """
    Attributes
    ----------
    boltzmann_factor_photo_ion : pandas.DataFrame, dtype float
    """

    outputs = ("boltzmann_factor_photo_ion",)

    def calculate(self, photo_ion_cross_sections, t_electrons):
        nu = photo_ion_cross_sections["nu"].values

        boltzmann_factor = np.exp(-nu[np.newaxis].T / t_electrons * (H / K_B))
        return boltzmann_factor


class CollIonRateCoeffSeaton(ProcessingPlasmaProperty):
    """
    Attributes
    ----------
    coll_ion_coeff : pandas.DataFrame, dtype float
        The rate coefficient for collisional ionization in the Seaton
        approximation. Multiply with the electron density and the
        level number density to obtain the total rate.

    Notes
    -----
    The rate coefficient for collisional ionization in the Seaton approximation
    is calculated according to Eq. 9.60 in [1].

    References
    ----------
    .. [1] Hubeny, I. and Mihalas, D., "Theory of Stellar Atmospheres". 2014.
    """

    outputs = ("coll_ion_coeff",)
    latex_name = (r"c_{\textrm{i,}\kappa}",)

    def calculate(self, photo_ion_cross_sections, t_electrons):
        photo_ion_cross_sections_threshold = photo_ion_cross_sections.groupby(
            level=[0, 1, 2]
        ).first()
        nu_i = photo_ion_cross_sections_threshold["nu"]
        factor = self._calculate_factor(nu_i, t_electrons)
        coll_ion_coeff = 1.55e13 * photo_ion_cross_sections_threshold["x_sect"]
        coll_ion_coeff = factor.multiply(coll_ion_coeff, axis=0)
        coll_ion_coeff = coll_ion_coeff.divide(np.sqrt(t_electrons), axis=1)

        ion_number = coll_ion_coeff.index.get_level_values("ion_number").values
        coll_ion_coeff[ion_number == 0] *= 0.1
        coll_ion_coeff[ion_number == 1] *= 0.2
        coll_ion_coeff[ion_number >= 2] *= 0.3
        return coll_ion_coeff

    def _calculate_factor(self, nu_i, t_electrons):
        u0s = self._calculate_u0s(nu_i.values, t_electrons)
        factor = np.exp(-u0s) / u0s
        factor = pd.DataFrame(factor, index=nu_i.index)
        return factor

    @staticmethod
    def _calculate_u0s(nu, t_electrons):
        u0s = nu[np.newaxis].T / t_electrons * (H / K_B)
        return u0s


class CollRecombRateCoeff(ProcessingPlasmaProperty):
    """
    Attributes
    ----------
    coll_recomb_coeff : pandas.DataFrame, dtype float
        The rate coefficient for collisional recombination.
        Multiply with the electron density squared and the ion number density
        to obtain the total rate.

    Notes
    -----
    The collisional recombination rate coefficient is calculated from the
    collisional ionization rate coefficient based on the requirement of detailed
    balance.
    """

    outputs = ("coll_recomb_coeff",)
    latex_name = (r"c_{\kappa\textrm{i,}}",)

    def calculate(self, phi_ik, coll_ion_coeff):
        return coll_ion_coeff.multiply(phi_ik.loc[coll_ion_coeff.index])


class RawCollIonTransProbs(TransitionProbabilitiesProperty, IndexSetterMixin):
    """
    Attributes
    ----------
    p_coll_ion : pandas.DataFrame, dtype float
        The unnormalized transition probabilities for
        collisional ionization.
    p_coll_recomb : pandas.DataFrame, dtype float
        The unnormalized transition probabilities for
        collisional recombination.
    cool_rate_coll_ion : pandas.DataFrame, dtype float
        The collisional ionization cooling rates of the electron gas.
    """

    outputs = ("p_coll_ion", "p_coll_recomb", "cool_rate_coll_ion")
    transition_probabilities_outputs = (
        "p_coll_ion",
        "p_coll_recomb",
        "cool_rate_coll_ion",
    )
    latex_name = (
        r"p^{\textrm{coll ion}}",
        r"p^{\textrm{coll recomb}}",
        r"C^{\textrm{ion}}",
    )

    def calculate(
        self,
        coll_ion_coeff,
        coll_recomb_coeff,
        nu_i,
        photo_ion_idx,
        electron_densities,
        energy_i,
        level_number_density,
    ):
        p_coll_ion = coll_ion_coeff.multiply(energy_i, axis=0)
        p_coll_ion = p_coll_ion.multiply(electron_densities, axis=1)
        p_coll_ion = self.set_index(p_coll_ion, photo_ion_idx, reverse=False)

        coll_recomb_rate = coll_recomb_coeff.multiply(
            electron_densities, axis=1
        )  # The full rate is obtained from this by multiplying by the
        # electron density and ion number density.
        p_recomb_deactivation = coll_recomb_rate.multiply(nu_i, axis=0) * H
        p_recomb_deactivation = self.set_index(
            p_recomb_deactivation, photo_ion_idx, transition_type=-1
        )
        p_recomb_deactivation = p_recomb_deactivation.groupby(level=[0]).sum()
        index_dd = pd.MultiIndex.from_product(
            [p_recomb_deactivation.index.values, ["k"], [0]],
            names=list(photo_ion_idx.columns) + ["transition_type"],
        )
        p_recomb_deactivation = p_recomb_deactivation.set_index(index_dd)

        p_recomb_internal = coll_recomb_rate.multiply(energy_i, axis=0)
        p_recomb_internal = self.set_index(
            p_recomb_internal, photo_ion_idx, transition_type=0
        )
        p_coll_recomb = pd.concat([p_recomb_deactivation, p_recomb_internal])

        cool_rate_coll_ion = (coll_ion_coeff * electron_densities).multiply(
            nu_i * H, axis=0
        )
        level_lower_index = coll_ion_coeff.index
        cool_rate_coll_ion = (
            cool_rate_coll_ion
            * level_number_density.loc[level_lower_index].values
        )
        cool_rate_coll_ion = self.set_index(
            cool_rate_coll_ion, photo_ion_idx, reverse=False
        )
        cool_rate_coll_ion = cool_rate_coll_ion.groupby(
            level="destination_level_idx"
        ).sum()
        ion_cool_index = pd.MultiIndex.from_product(
            [["k"], cool_rate_coll_ion.index.values, [0]],
            names=list(photo_ion_idx.columns) + ["transition_type"],
        )
        cool_rate_coll_ion = cool_rate_coll_ion.set_index(ion_cool_index)
        return p_coll_ion, p_coll_recomb, cool_rate_coll_ion
