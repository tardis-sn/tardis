import numpy as np
import pandas as pd
from astropy import units as u
from scipy.special import exp1

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


class CollisionalCrossSectionYG(CollisionalCrossSections):
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

    def calculate(self, atomic_data, continuum_interaction_species):
        yg_data = atomic_data.yg_data
        if yg_data is None:
            raise ValueError(
                "Tardis does not support continuum interactions for atomic data sources that do not contain yg_data"
            )

        mask_selected_species = yg_data.index.droplevel(
            ["level_number_lower", "level_number_upper"]
        ).isin(continuum_interaction_species)
        yg_data = yg_data[mask_selected_species]

        t_yg = atomic_data.collision_data_temperatures
        yg_data.columns = t_yg
        approximate_yg_data = self.calculate_yg_van_regemorter(
            atomic_data, t_yg, continuum_interaction_species
        )

        yg_data = yg_data.combine_first(approximate_yg_data)

        energies = atomic_data.levels.energy
        index = yg_data.index
        lu_index = index.droplevel("level_number_lower")
        ll_index = index.droplevel("level_number_upper")
        delta_E = energies.loc[lu_index].values - energies.loc[ll_index].values
        delta_E = pd.Series(delta_E, index=index)

        source_idx = atomic_data.macro_atom_references.loc[
            ll_index
        ].references_idx
        destination_idx = atomic_data.macro_atom_references.loc[
            lu_index
        ].references_idx
        yg_idx = pd.DataFrame(
            {
                "source_level_idx": source_idx.values,
                "destination_level_idx": destination_idx.values,
            },
            index=index,
        )
        return yg_data, t_yg, index, delta_E, yg_idx

    @classmethod
    def calculate_yg_van_regemorter(
        cls, atomic_data, t_electrons, continuum_interaction_species
    ):
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
        HYDROGEN_IONIZATION_ENERGY = atomic_data.ionization_data.loc[(1, 1)]

        mask_selected_species = atomic_data.lines.index.droplevel(
            ["level_number_lower", "level_number_upper"]
        ).isin(continuum_interaction_species)
        lines_filtered = atomic_data.lines[mask_selected_species]
        f_lu = lines_filtered.f_lu.values
        nu_lines = lines_filtered.nu.values

        coll_const = A0**2 * np.pi * np.sqrt(8 * K_B / (np.pi * M_E))
        yg = 14.5 * coll_const * t_electrons * yg[:, np.newaxis]

        u0 = nu_lines[np.newaxis].T / t_electrons * (H / K_B)
        gamma = 0.276 * cls.exp1_times_exp(u0)
        gamma[gamma < 0.2] = 0.2
        yg *= u0 * gamma / BETA_COLL
        yg = pd.DataFrame(yg, index=lines_filtered.index, columns=t_electrons)

        return yg

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


class CollisionCrossSectionRegemorter:
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
            * t_electrons
            * collision_cross_section[:, np.newaxis]
        )

        u0 = (
            self.transition_data.nu[np.newaxis].T
            / t_electrons
            * (const.h / const.k_B)
        )
        gamma = 0.276 * exp1_times_exp(u0)
        gamma[gamma < 0.2] = 0.2
        collision_cross_section *= u0 * gamma / BETA_COLL
        collision_cross_section = pd.DataFrame(
            collision_cross_section,
            index=self.transition_data.index,
            columns=t_electrons,
        )

        return collision_cross_section
