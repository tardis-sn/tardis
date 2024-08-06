import numpy as np
import pandas as pd
from astropy import units as u
from scipy.special import exp1
from scipy.interpolate import PchipInterpolator

from tardis import constants as const


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


REGEMORTER_CONSTANT = (  #  Hubeny, I. and Mihalas, D., "Theory of Stellar Atmospheres". 2014. EQ 9.54 [below it]
    const.a0.cgs**2
    * np.pi
    * np.sqrt(8 * const.k_B.cgs / (np.pi * const.m_e.cgs))
)

HYDROGEN_IONIZATION_ENERGY = (
    13.598434005136003 * u.eV
).cgs  # taken from the classic TARDIS ionization data


class CollisionalCrossSections:
    def __init__(self, collision_cross_sections):
        self.collisional_cross_sections = collision_cross_sections

    def solve_collisional_cross_sections(self, temperature_electron):
        pass


N_A = const.N_A.cgs.value
K_B = const.k_B.cgs.value
C = const.c.cgs.value
H = const.h.cgs.value
A0 = const.a0.cgs.value
M_E = const.m_e.cgs.value
E = const.e.esu.value
BETA_COLL = (
    (const.h**4 / (8 * const.k_B * const.m_e**3 * np.pi**3)) ** 0.5
).cgs
F_K = (
    16
    / (3.0 * np.sqrt(3))
    * np.sqrt((2 * np.pi) ** 3 * K_B / (H**2 * M_E**3))
    * (E**2 / C) ** 3
)  # See Eq. 19 in Sutherland, R. S. 1998, MNRAS, 300, 321
FF_OPAC_CONST = (
    (2 * np.pi / (3 * M_E * K_B)) ** 0.5 * 4 * E**6 / (3 * M_E * H * C)
)  # See Eq. 6.1.8 in http://personal.psu.edu/rbc3/A534/lec6.pdf


class YGSolver(CollisionalCrossSections):
    """
    Attributes
    ----------
    yg_data : pandas.DataFrame
        Table of thermally averaged effective collision strengths
        (divided by the statistical weight of the lower level) Y_ij / g_i .
        Columns are temperatures.
    t_yg : numpy.ndarray
        Temperatures at which collision strengths are tabulated.
    yg_index : Pandas MultiIndex
    delta_E_yg : pandas.DataFrame
        Energy difference between upper and lower levels coupled by collisions.
    yg_idx : pandas.DataFrame
        Source_level_idx and destination_level_idx of collision transitions.
        Indexed by atomic_number, ion_number, level_number_lower,
        level_number_upper.
    """

    def __init__(self, yg_data, yg_temperature_data, delta_energies):
        yg_data.columns = yg_temperature_data
        self.yg_data = yg_data
        self.delta_energies = delta_energies

        yg_idx = pd.DataFrame(
            {
                "source_level_idx": source_idx.values,
                "destination_level_idx": destination_idx.values,
            },
            index=index,
        )
        self.yg_interpolator = PchipInterpolator(
            self.yg_data.columns, self.yg_data.values, axis=1, extrapolate=True
        )

    def solve(self, t_electrons):
        yg = self.yg_interpolator(t_electrons)

        boltzmann_factor = np.exp(
            -self.delta_energies.values[np.newaxis].T
            / (t_electrons * const.k_B).value
        )

        q_ij = (
            BETA_COLL / np.sqrt(t_electrons) * yg * boltzmann_factor
        )  # see formula A2 in Przybilla, Butler 2004 - Apj 609, 1181
        return pd.DataFrame(q_ij, index=self.delta_energies)


class YGRegemorterSolver:
    def __init__(self, transition_data) -> None:
        assert transition_data.index.names == [
            "atomic_number",
            "ion_number",
            "level_number_lower",
            "level_number_upper",
        ]
        assert {"f_lu", "nu"} - set(transition_data.columns) == set()

        assert np.all(
            transition_data.index.get_level_values("level_number_lower")
            < transition_data.index.get_level_values("level_number_upper")
        )
        self.transition_data = transition_data.sort_index()

    def solve(self, t_electrons):
        """
        Calculate collision strengths in the van Regemorter approximation.

        This function calculates thermally averaged effective collision
        strengths (divided by the statistical weight of the lower level)
        Y_ij / g_i using the van Regemorter approximation.

        Parameters
        ----------
        atomic_data : tardis.io.atom_data.AtomData
        t_electrons : numpy.ndarray
        continuum_interaction_species : pandas.MultiIndex

        Returns
        -------
        pandas.DataFrame
            Thermally averaged effective collision strengths
            (divided by the statistical weight of the lower level) Y_ij / g_i

        Notes
        -----
        See Eq. 9.58 in [2].

        References
        ----------
        .. [1] van Regemorter, H., “Rate of Collisional Excitation in Stellar
               Atmospheres.”, The Astrophysical Journal, vol. 136, p. 906, 1962.
               doi:10.1086/147445.
        .. [2] Hubeny, I. and Mihalas, D., "Theory of Stellar Atmospheres". 2014.
        """
        collision_cross_section = (
            self.transition_data.f_lu.values
            * (
                HYDROGEN_IONIZATION_ENERGY
                / (const.h * self.transition_data.nu.values * u.Hz)
            )
            ** 2
        )

        collision_cross_section = (
            14.5
            * REGEMORTER_CONSTANT
            * t_electrons.value
            * collision_cross_section[:, np.newaxis]
        )

        u0 = (
            const.h.cgs.value * self.transition_data.nu.values[np.newaxis].T
        ) / (t_electrons.value * const.k_B.cgs.value)
        gamma = 0.276 * exp1_times_exp(u0)
        gamma[gamma < 0.2] = 0.2
        collision_cross_section *= u0 * gamma / BETA_COLL
        collision_cross_section = pd.DataFrame(
            collision_cross_section.cgs.value,
            index=self.transition_data.index,
            columns=t_electrons.value,
        )

        return collision_cross_section


class CollisionalRatesSolver:
    def __init__(self, transition_data):
        pass


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
            BETA_COLL.value / np.sqrt(t_electrons) * yg * boltzmann_factor
        )  # see formula A2 in Przybilla, Butler 2004 - Apj 609, 1181
        return pd.DataFrame(q_ij, index=yg_index)


class CollDeexcRateCoeff:
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
