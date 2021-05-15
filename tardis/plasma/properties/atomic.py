import logging

import numpy as np
import pandas as pd
from scipy.special import expn
from scipy.interpolate import PchipInterpolator
from collections import Counter as counter
from tardis import constants as const

from tardis.plasma.properties.base import (
    ProcessingPlasmaProperty,
    HiddenPlasmaProperty,
    BaseAtomicDataProperty,
)
from tardis.plasma.exceptions import IncompleteAtomicData
from tardis.plasma.properties.continuum_processes import (
    get_ground_state_multi_index,
    K_B,
    BETA_COLL,
    H,
    A0,
    M_E,
    C,
)

logger = logging.getLogger(__name__)

__all__ = [
    "Levels",
    "Lines",
    "LinesLowerLevelIndex",
    "LinesUpperLevelIndex",
    "AtomicMass",
    "IonizationData",
    "ZetaData",
    "NLTEData",
    "PhotoIonizationData",
    "YgData",
    "YgInterpolator",
    "LevelIdxs2LineIdx",
    "TwoPhotonData",
]


class Levels(BaseAtomicDataProperty):
    """
    Attributes
    ----------
    levels : pandas.MultiIndex
        (atomic_number, ion_number, level_number)
        Index of filtered atomic data. Index used for all other attribute dataframes for this class
    excitation_energy : pandas.DataFrame, dtype float
        Excitation energies of atomic levels.
        Index is levels.
    metastability : pandas.DataFrame, dtype bool
        Records whether atomic levels are metastable.
        Index is levels.
    g : pandas.DataFrame (index=levels), dtype float
        Statistical weights of atomic levels.
    """

    outputs = ("levels", "excitation_energy", "metastability", "g")
    latex_name = (
        r"\textrm{levels}",
        r"\epsilon_{\textrm{k}}",
        r"\textrm{metastability}",
        "g",
    )

    def _filter_atomic_property(self, levels, selected_atoms):
        return levels
        # return levels[levels.atomic_number.isin(selected_atoms)]

    def _set_index(self, levels):
        # levels = levels.set_index(['atomic_number', 'ion_number',
        #                          'level_number'])
        return (
            levels.index,
            levels["energy"],
            levels["metastable"],
            levels["g"],
        )


class Lines(BaseAtomicDataProperty):
    """
    Attributes
    ----------
    lines : pandas.DataFrame
        Atomic lines data. Columns are wavelength, atomic_number,ion_number,
        f_ul, f_lu, level_number_lower, level_number_upper, nu, B_lu, B_ul, A_ul,
        wavelength. Index is line_id.
    nu : pandas.DataFrame, dtype float
        Line frequency data. Index is line_id.
    f_lu : pandas.DataFrame, dtype float
        Transition probability data. Index is line_id.
    wavelength_cm : pandas.DataFrame, dtype float
        Line wavelengths in cm. Index is line_id.
    """

    # Would like for lines to just be the line_id values
    outputs = ("lines", "nu", "f_lu", "wavelength_cm")

    def _filter_atomic_property(self, lines, selected_atoms):
        # return lines[lines.atomic_number.isin(selected_atoms)]
        return lines

    def _set_index(self, lines):
        # lines.set_index('line_id', inplace=True)
        return lines, lines["nu"], lines["f_lu"], lines["wavelength_cm"]


class PhotoIonizationData(ProcessingPlasmaProperty):
    """
    Attributes
    ----------
    photo_ion_cross_sections : pandas.DataFrame, dtype float
        Photoionization cross sections as a function of frequency.
        Columns are nu, x_sect, index=('atomic_number','ion_number','level_number')
    photo_ion_block_references : numpy.ndarray, dtype int
        Indices where the photoionization data for
        a given level starts. Needed for calculation
        of recombination rates.
    nu_i : pandas.Series, dtype float
        Threshold frequencies for ionization
    energy_i : pandas.Series, dtype float
        Energies of levels with bound-free transitions. Needed to calculate
        for example internal transition probabilities in the macro atom scheme.
    photo_ion_index : pandas.MultiIndex, dtype int
        Atomic, ion and level numbers for which photoionization data exists.
    level2continuum_idx : pandas.Series, dtype int
        Maps a level MultiIndex (atomic_number, ion_number, level_number) to
        the continuum_idx of the corresponding bound-free continuum (which are
        sorted by decreasing frequency).
    """

    outputs = (
        "photo_ion_cross_sections",
        "photo_ion_block_references",
        "photo_ion_index",
        "nu_i",
        "energy_i",
        "photo_ion_idx",
        "level2continuum_idx",
    )
    latex_name = (
        r"\xi_{\textrm{i}}(\nu)",
        "",
        "",
        r"\nu_i",
        r"\epsilon_i",
        "",
    )

    def calculate(self, atomic_data, continuum_interaction_species):
        photoionization_data = atomic_data.photoionization_data.set_index(
            ["atomic_number", "ion_number", "level_number"]
        )
        mask_selected_species = photoionization_data.index.droplevel(
            "level_number"
        ).isin(continuum_interaction_species)
        photoionization_data = photoionization_data[mask_selected_species]
        phot_nus = photoionization_data["nu"]
        block_references = np.pad(
            phot_nus.groupby(level=[0, 1, 2]).count().values.cumsum(), [1, 0]
        )
        photo_ion_index = photoionization_data.index.unique()
        nu_i = photoionization_data.groupby(level=[0, 1, 2]).first().nu
        energy_i = atomic_data.levels.loc[photo_ion_index].energy

        source_idx = atomic_data.macro_atom_references.loc[
            photo_ion_index
        ].references_idx
        destination_idx = atomic_data.macro_atom_references.loc[
            get_ground_state_multi_index(photo_ion_index)
        ].references_idx
        photo_ion_idx = pd.DataFrame(
            {
                "source_level_idx": source_idx.values,
                "destination_level_idx": destination_idx.values,
            },
            index=photo_ion_index,
        )

        level2continuum_edge_idx = pd.Series(
            np.arange(len(nu_i)),
            nu_i.sort_values(ascending=False).index,
            name="continuum_idx",
        )
        return (
            photoionization_data,
            block_references,
            photo_ion_index,
            nu_i,
            energy_i,
            photo_ion_idx,
            level2continuum_edge_idx,
        )


class TwoPhotonData(ProcessingPlasmaProperty):
    outputs = ("two_photon_data", "two_photon_idx")
    """
    Attributes:
    two_photon_data : pandas.DataFrame, dtype float
        A DataFrame containing the *two photon decay data* with:
            index: atomic_number, ion_number, level_number_lower, level_number_upper
            columns: A_ul[1/s], nu0[Hz], alpha, beta, gamma
        alpha, beta, gamma are fit coefficients for the frequency dependent
        transition probability A(y) of the two photon decay. See Eq. 2 in
        Nussbaumer & Schmutz (1984).
    two_photon_idx : pandas.DataFrame, dtype int
    """

    def calculate(self, atomic_data, continuum_interaction_species):
        two_photon_data = atomic_data.two_photon_data
        mask_selected_species = two_photon_data.index.droplevel(
            ["level_number_lower", "level_number_upper"]
        ).isin(continuum_interaction_species)
        if not mask_selected_species.sum():
            raise IncompleteAtomicData(
                "two photon transition data for the requested "
                "continuum_interactions species: {}".format(
                    continuum_interaction_species.values.tolist()
                )
            )
        two_photon_data = two_photon_data[mask_selected_species]
        index_lower = two_photon_data.index.droplevel("level_number_upper")
        index_upper = two_photon_data.index.droplevel("level_number_lower")
        source_idx = atomic_data.macro_atom_references.loc[
            index_upper
        ].references_idx
        destination_idx = atomic_data.macro_atom_references.loc[
            index_lower
        ].references_idx
        two_photon_idx = pd.DataFrame(
            {
                "source_level_idx": source_idx.values,
                "destination_level_idx": destination_idx.values,
            },
            index=two_photon_data.index,
        )
        return two_photon_data, two_photon_idx


class LinesLowerLevelIndex(HiddenPlasmaProperty):
    """
    Attributes
    ----------
    lines_lower_level_index : numpy.ndrarray, dtype int
        Levels data for lower levels of particular lines
    """

    outputs = ("lines_lower_level_index",)

    def calculate(self, levels, lines):
        levels_index = pd.Series(
            np.arange(len(levels), dtype=np.int64), index=levels
        )
        lines_index = lines.index.droplevel("level_number_upper")
        return np.array(levels_index.loc[lines_index])


class LinesUpperLevelIndex(HiddenPlasmaProperty):
    """
    Attributes
    ----------
    lines_upper_level_index : numpy.ndarray, dtype int
        Levels data for upper levels of particular lines
    """

    outputs = ("lines_upper_level_index",)

    def calculate(self, levels, lines):
        levels_index = pd.Series(
            np.arange(len(levels), dtype=np.int64), index=levels
        )
        lines_index = lines.index.droplevel("level_number_lower")
        return np.array(levels_index.loc[lines_index])


class LevelIdxs2LineIdx(HiddenPlasmaProperty):
    """
    Attributes:
    level_idxs2line_idx : pandas.Series, dtype int
       Maps a source_level_idx destination_level_idx pair to a line_idx.
    """

    outputs = ("level_idxs2line_idx",)

    def calculate(self, atomic_data):
        index = pd.MultiIndex.from_arrays(
            [
                atomic_data.lines_upper2level_idx,
                atomic_data.lines_lower2level_idx,
            ],
            names=["source_level_idx", "destination_level_idx"],
        )
        level_idxs2line_idx = pd.Series(
            np.arange(len(index)), index=index, name="lines_idx"
        )
        return level_idxs2line_idx


class AtomicMass(ProcessingPlasmaProperty):
    """
    Attributes
    ----------
    atomic_mass : pandas.Series
        Atomic masses of the elements used. Indexed by atomic number.
    """

    outputs = ("atomic_mass",)

    def calculate(self, atomic_data, selected_atoms):
        if getattr(self, self.outputs[0]) is not None:
            return (getattr(self, self.outputs[0]),)
        else:
            return atomic_data.atom_data.loc[selected_atoms].mass


class IonizationData(BaseAtomicDataProperty):
    """
    Attributes
    ----------
    ionization_data : pandas.Series
        Holding ionization energies
        Indexed by atomic number, ion number.
    """

    outputs = ("ionization_data",)

    def _filter_atomic_property(self, ionization_data, selected_atoms):
        mask = ionization_data.index.isin(selected_atoms, level="atomic_number")
        ionization_data = ionization_data[mask]
        counts = ionization_data.groupby(level="atomic_number").count()

        if np.alltrue(counts.index == counts):
            return ionization_data
        else:
            raise IncompleteAtomicData(
                f"ionization data for the ion ({str(counts.index[counts.index != counts])}, {str(counts[counts.index != counts])})"
            )

    def _set_index(self, ionization_data):
        return ionization_data


class ZetaData(BaseAtomicDataProperty):
    """
    Attributes
    ----------
    zeta_data : pandas.DataFrame, dtype float
        Zeta data for the elements used. Indexed by atomic number, ion number.
        Columns are temperature values up to 40,000 K in iterations of 2,000 K.
        The zeta value represents the fraction of recombination events
        from the ionized state that go directly to the ground state.
    """

    outputs = ("zeta_data",)

    def _filter_atomic_property(self, zeta_data, selected_atoms):
        zeta_data["atomic_number"] = zeta_data.index.codes[0] + 1
        zeta_data["ion_number"] = zeta_data.index.codes[1] + 1
        zeta_data = zeta_data[zeta_data.atomic_number.isin(selected_atoms)]
        zeta_data_check = counter(zeta_data.atomic_number.values)
        keys = np.array(list(zeta_data_check.keys()))
        values = np.array(zeta_data_check.values())
        if np.alltrue(keys + 1 == values):
            return zeta_data
        else:
            #            raise IncompleteAtomicData('zeta data')
            # This currently replaces missing zeta data with 1, which is necessary with
            # the present atomic data. Will replace with the error above when I have
            # complete atomic data.
            missing_ions = []
            updated_index = []
            for atom in selected_atoms:
                for ion in range(1, atom + 2):
                    if (atom, ion) not in zeta_data.index:
                        missing_ions.append((atom, ion))
                    updated_index.append([atom, ion])
            logger.warn(
                f"Zeta_data missing - replaced with 1s. Missing ions: {missing_ions}"
            )
            updated_index = np.array(updated_index)
            updated_dataframe = pd.DataFrame(
                index=pd.MultiIndex.from_arrays(
                    updated_index.transpose().astype(int)
                ),
                columns=zeta_data.columns,
            )
            for value in range(len(zeta_data)):
                updated_dataframe.loc[
                    zeta_data.atomic_number.values[value],
                    zeta_data.ion_number.values[value],
                ] = zeta_data.loc[
                    zeta_data.atomic_number.values[value],
                    zeta_data.ion_number.values[value],
                ]
            updated_dataframe = updated_dataframe.astype(float)
            updated_index = pd.DataFrame(updated_index)
            updated_dataframe["atomic_number"] = np.array(updated_index[0])
            updated_dataframe["ion_number"] = np.array(updated_index[1])
            updated_dataframe.fillna(1.0, inplace=True)
            return updated_dataframe

    def _set_index(self, zeta_data):
        return zeta_data.set_index(["atomic_number", "ion_number"])


class NLTEData(ProcessingPlasmaProperty):
    """
    Attributes
    ----------
    nlte_data :
        #Finish later (need atomic dataset with NLTE data).
    """

    outputs = ("nlte_data",)

    def calculate(self, atomic_data):
        if getattr(self, self.outputs[0]) is not None:
            return (getattr(self, self.outputs[0]),)
        else:
            return atomic_data.nlte_data


class YgData(ProcessingPlasmaProperty):
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

    outputs = ("yg_data", "t_yg", "yg_index", "delta_E_yg", "yg_idx")
    latex_name = (
        r"\frac{Y_{ij}}{g_i}",
        r"T_\textrm{Yg}",
        r"\textrm{yg_index}",
        r"\delta E_{ij}",
        r"\textrm{yg_idx}",
    )

    def calculate(self, atomic_data, continuum_interaction_species):
        yg_data = atomic_data.yg_data

        mask_selected_species = yg_data.index.droplevel(
            ["level_number_lower", "level_number_upper"]
        ).isin(continuum_interaction_species)
        yg_data = yg_data[mask_selected_species]

        t_yg = yg_data.columns.values.astype(float)
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

    @staticmethod
    def calculate_yg_van_regemorter(
        atomic_data, t_electrons, continuum_interaction_species
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
        I_H = atomic_data.ionization_data.loc[(1, 1)]

        mask_selected_species = atomic_data.lines.index.droplevel(
            ["level_number_lower", "level_number_upper"]
        ).isin(continuum_interaction_species)
        lines_filtered = atomic_data.lines[mask_selected_species]
        f_lu = lines_filtered.f_lu.values
        nu_lines = lines_filtered.nu.values

        yg = f_lu * (I_H / (H * nu_lines)) ** 2
        coll_const = A0 ** 2 * np.pi * np.sqrt(8 * K_B / (np.pi * M_E))
        yg = 14.5 * coll_const * t_electrons * yg[:, np.newaxis]

        u0 = nu_lines[np.newaxis].T / t_electrons * (H / K_B)
        gamma = 0.276 * np.exp(u0) * expn(1, u0)
        gamma[gamma < 0.2] = 0.2
        yg *= u0 * gamma / BETA_COLL
        yg = pd.DataFrame(yg, index=lines_filtered.index, columns=t_electrons)

        return yg


class YgInterpolator(ProcessingPlasmaProperty):
    """
    Attributes
    ----------
    yg_interp : scipy.interpolate.PchipInterpolator
        Interpolates the thermally averaged effective collision strengths
        (divided by the statistical weight of the lower level) Y_ij / g_i as
        a function of electron temperature.
    """

    outputs = ("yg_interp",)
    latex_name = ("\\frac{Y_ij}{g_i}_{\\textrm{interp}}",)

    def calculate(self, yg_data, t_yg):
        yg_interp = PchipInterpolator(t_yg, yg_data, axis=1, extrapolate=True)
        return yg_interp
