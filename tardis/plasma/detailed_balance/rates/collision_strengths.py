import numpy as np
import pandas as pd
from astropy import units as u
from scipy.special import exp1
from scipy.interpolate import PchipInterpolator, splrep, splev

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
K_B_EV = const.k_B.cgs.to("eV / K").value
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


def calculate_upsilon_g_2_collisional_rates(yg, t_electrons, delta_energies):
    boltzmann_factor = np.exp(
        -delta_energies.values[np.newaxis].T / (t_electrons * const.k_B).value
    )

    q_lu = (
        BETA_COLL.value / np.sqrt(t_electrons) * yg * boltzmann_factor
    )  # see formula A2 in Przybilla, Butler 2004 - Apj 609, 1181
    return pd.DataFrame(q_lu, index=delta_energies.index)


class UpsilonCMFGENSolver:
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

    def __init__(
        self,
        upsilon_temperatures,
        upsilon_g_data,
    ):
        self.upsilon_lu_data = upsilon_g_data

        # can produce upsilon/g or not, depending on how easy it is
        self.upsilon_g_lu_interpolator = PchipInterpolator(
            upsilon_temperatures,
            self.upsilon_lu_data.values,
            axis=1,
            extrapolate=True,
        )

    def solve(self, t_electrons):
        return pd.DataFrame(
            self.upsilon_g_lu_interpolator(t_electrons),
            index=self.upsilon_lu_data.index,
        )


class UpsilonChiantiSolver:
    """Solver for Upsilon / g_i for Chianti data."""

    def __init__(
        self,
        upsilon_data,
    ):
        self.upsilon_lu_data = upsilon_data

    def upsilon_scaling(self, row, t_electrons):
        """Scales Upsilon from Chianti data using equations
        23-38 from Burgess & Tully 1992 - A&A 254, 436B.

        Parameters
        ----------
        row : pd.Series
            DataFrame row of Chianti collisional data
        t_electrons : np.ndarray
            1D array of electron temperatures to interpolate over

        Returns
        -------
        pd.Series
            Scaled Upsilon / g_lower

        Raises
        ------
        ValueError
            Incorrect scaling type provided
        """
        scaling_constant = row["cups"]
        x_knots = np.linspace(0, 1, len(row["btemp"]))
        y_knots = row["bscups"]
        delta_energy = row["delta_e"]
        g_lower = row["g_l"]

        scaling_type = row["ttype"]
        if scaling_type > 5:
            scaling_type -= 5

        kt = K_B_EV * t_electrons

        spline_tck = splrep(x_knots, y_knots)

        if scaling_type == 1:
            x = 1 - np.log(scaling_constant) / np.log(
                kt / delta_energy + scaling_constant
            )
            y_func = splev(x, spline_tck)
            upsilon = y_func * np.log(kt / delta_energy + np.exp(1))

        elif scaling_type == 2:
            x = (kt / delta_energy) / (kt / delta_energy + scaling_constant)
            y_func = splev(x, spline_tck)
            upsilon = y_func

        elif scaling_type == 3:
            x = (kt / delta_energy) / (kt / delta_energy + scaling_constant)
            y_func = splev(x, spline_tck)
            upsilon = y_func / (kt / delta_energy + 1)

        elif scaling_type == 4:
            x = 1 - np.log(scaling_constant) / np.log(
                kt / delta_energy + scaling_constant
            )
            y_func = splev(x, spline_tck)
            upsilon = y_func * np.log(kt / delta_energy + scaling_constant)

        elif scaling_type > 4:
            raise ValueError(
                "Not sure what to do with scaling type greater than 4"
            )

        upsilon_g_lu = upsilon / g_lower
        return pd.Series(data=upsilon_g_lu, name="upsilon_g")

    def solve(self, t_electrons):
        """Solve the Upsilon / g_lower collisional values for arbitrary temperatures.

        Parameters
        ----------
        t_electrons : np.ndarray
            1D array of electron temperatures to interpolate over

        Returns
        -------
        pd.DataFrame
            DataFrame with columns of Upsilon / g_lower per transition and temperature.
        """
        upsilon_g_lu = self.upsilon_lu_data.apply(
            self.upsilon_scaling,
            axis=1,
            args=(t_electrons.value,),
        )
        return pd.DataFrame(
            upsilon_g_lu,
            index=self.upsilon_lu_data.index,
        )


class UpsilonRegemorterSolver:
    def __init__(self, transition_data, g_bar=0.2) -> None:
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
        self.g_bar = g_bar

    def solve(self, t_electrons):
        """
        Calculate collision strengths in the van Regemorter approximation.

        This function calculates thermally averaged effective collision
        strengths (divided by the statistical weight of the lower level)
        Y_ij / g_i using the van Regemorter approximation. A very good description can be found in
        Mihalas Chapter on collisional rates

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
        upsilon_g_lu = (
            self.transition_data.f_lu.values
            * (
                HYDROGEN_IONIZATION_ENERGY
                / (const.h * self.transition_data.nu.values * u.Hz)
            )
            ** 2
        )

        upsilon_g_lu = (
            14.5
            * REGEMORTER_CONSTANT
            * t_electrons.value
            * upsilon_g_lu[:, np.newaxis]
        )

        u0 = (
            const.h.cgs.value * self.transition_data.nu.values[np.newaxis].T
        ) / (t_electrons.value * const.k_B.cgs.value)
        gamma_component = 0.276 * exp1_times_exp(u0)  # Eq 9.59 in Mihalas
        # choice of transitions between principal quantum numbers g_bar = 0.2, otherwise gbar = 0.7
        # NOTE currently we assume all transitions have changes in principal quantum numbers which is wrong
        gamma = np.maximum(self.g_bar, gamma_component)
        upsilon_g_lu *= u0 * gamma / BETA_COLL
        upsilon_g_lu = pd.DataFrame(
            upsilon_g_lu.cgs.value,
            index=self.transition_data.index,
        )
        return upsilon_g_lu


class CollExcRateCoeff:
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
