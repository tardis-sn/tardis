import logging

import numpy as np
import pandas as pd
from numba import njit
from scipy.special import exp1
from scipy.interpolate import PchipInterpolator
from collections import Counter as counter
import radioactivedecay as rd
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
    "IsotopeMass",
    "IonizationData",
    "ZetaData",
    "NLTEData",
    "MacroAtomData",
    "PhotoIonizationData",
    "YgData",
    "YgInterpolator",
    "LevelIdxs2LineIdx",
    "LevelIdxs2TransitionIdx",
    "TwoPhotonData",
    "ContinuumInteractionHandler",
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
        return levels[levels.index.isin(selected_atoms, level="atomic_number")]

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

    latex_name = (
        r"\textrm{lines}",
        r"\nu",
        r"f_lu",
        r"\lambda_{cm}",
    )

    def _filter_atomic_property(self, lines, selected_atoms):
        # return lines[lines.atomic_number.isin(selected_atoms)]
        return lines

    def _set_index(self, lines):
        # lines.set_index('line_id', inplace=True)
        return lines, lines["nu"], lines["f_lu"], lines["wavelength_cm"]


class MacroAtomData(BaseAtomicDataProperty):
    outputs = ("macro_atom_data",)

    def _filter_atomic_property(self, macro_atom_data, selected_atoms):
        return macro_atom_data

    def _set_index(self, macro_atom_data):
        return macro_atom_data


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
    level_idxs2continuum_idx : pandas.DataFrame, dtype int
        Maps a source_level_idx destination_level_idx pair to a continuum_idx.
    """

    outputs = (
        "photo_ion_cross_sections",
        "photo_ion_block_references",
        "photo_ion_index",
        "nu_i",
        "energy_i",
        "photo_ion_idx",
        "level2continuum_idx",
        "level_idxs2continuum_idx",
    )
    latex_name = (
        r"\xi_{\textrm{i}}(\nu)",
        "",
        "",
        r"\nu_i",
        r"\epsilon_i",
        "",
        "",
    )

    def calculate(self, atomic_data, continuum_interaction_species):
        # photoionization_data = atomic_data.photoionization_data.set_index(
        #    ["atomic_number", "ion_number", "level_number"]
        # )
        photoionization_data = atomic_data.photoionization_data
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

        level_idxs2continuum_idx = photo_ion_idx.copy()
        level_idxs2continuum_idx["continuum_idx"] = level2continuum_edge_idx
        level_idxs2continuum_idx = level_idxs2continuum_idx.set_index(
            ["source_level_idx", "destination_level_idx"]
        )
        return (
            photoionization_data,
            block_references,
            photo_ion_index,
            nu_i,
            energy_i,
            photo_ion_idx,
            level2continuum_edge_idx,
            level_idxs2continuum_idx,
        )


class ContinuumInteractionHandler(ProcessingPlasmaProperty):
    outputs = (
        "get_current_bound_free_continua",
        "determine_bf_macro_activation_idx",
        "determine_continuum_macro_activation_idx",
    )

    def calculate(
        self,
        photo_ion_cross_sections,
        level2continuum_idx,
        photo_ion_idx,
        k_packet_idx,
    ):
        nus = photo_ion_cross_sections.nu.loc[
            level2continuum_idx.index
        ]  # Sort by descending frequency
        nu_mins = nus.groupby(level=[0, 1, 2], sort=False).first().values
        nu_maxs = nus.groupby(level=[0, 1, 2], sort=False).last().values

        @njit(error_model="numpy", fastmath=True)
        def get_current_bound_free_continua(nu):
            """
            Determine bound-free continua for which absorption is possible.

            Parameters
            ----------
            nu : float
                Comoving frequency of the r-packet.

            Returns
            -------
            numpy.ndarray, dtype int
                Continuum ids for which absorption is possible for frequency `nu`.
            """
            # searchsorted would be faster but would need stricter format for photoionization data
            current_continua = np.where(
                np.logical_and(nu >= nu_mins, nu <= nu_maxs)
            )[0]
            return current_continua

        destination_level_idxs = photo_ion_idx.loc[
            level2continuum_idx.index, "destination_level_idx"
        ].values

        @njit(error_model="numpy", fastmath=True)
        def determine_bf_macro_activation_idx(
            nu, chi_bf_contributions, active_continua
        ):
            """
            Determine the macro atom activation level after bound-free absorption.

            Parameters
            ----------
            nu : float
                Comoving frequency of the r-packet.
            chi_bf_contributions : numpy.ndarray, dtype float
                Cumulative distribution of bound-free opacities at frequency
                `nu`.
            active_continua : numpy.ndarray, dtype int
                Continuum ids for which absorption is possible for frequency `nu`.

            Returns
            -------
            float
                Macro atom activation idx.
            """
            # Perform a MC experiment to determine the continuum for absorption
            index = np.searchsorted(chi_bf_contributions, np.random.random())
            continuum_id = active_continua[index]

            # Perform a MC experiment to determine whether thermal or
            # ionization energy is created
            nu_threshold = nu_mins[continuum_id]
            fraction_ionization = nu_threshold / nu
            if (
                np.random.random() < fraction_ionization
            ):  # Create ionization energy (i-packet)
                destination_level_idx = destination_level_idxs[continuum_id]
            else:  # Create thermal energy (k-packet)
                destination_level_idx = k_packet_idx
            return destination_level_idx

        @njit(error_model="numpy", fastmath=True)
        def determine_continuum_macro_activation_idx(
            nu, chi_bf, chi_ff, chi_bf_contributions, active_continua
        ):
            """
            Determine the macro atom activation level after a continuum absorption.

            Parameters
            ----------
            nu : float
                Comoving frequency of the r-packet.
            chi_bf : numpy.ndarray, dtype float
                Bound-free opacity.
            chi_bf : numpy.ndarray, dtype float
                Free-free opacity.
            chi_bf_contributions : numpy.ndarray, dtype float
                Cumulative distribution of bound-free opacities at frequency
                `nu`.
            active_continua : numpy.ndarray, dtype int
                Continuum ids for which absorption is possible for frequency `nu`.

            Returns
            -------
            float
                Macro atom activation idx.
            """
            fraction_bf = chi_bf / (chi_bf + chi_ff)
            # TODO: In principle, we can also decide here whether a Thomson
            # scattering event happens and need one less RNG call.
            if np.random.random() < fraction_bf:  # Bound-free absorption
                destination_level_idx = determine_bf_macro_activation_idx(
                    nu, chi_bf_contributions, active_continua
                )
            else:  # Free-free absorption (i.e. k-packet creation)
                destination_level_idx = k_packet_idx
            return destination_level_idx

        return (
            get_current_bound_free_continua,
            determine_bf_macro_activation_idx,
            determine_continuum_macro_activation_idx,
        )


class TwoPhotonData(ProcessingPlasmaProperty):
    outputs = ("two_photon_data", "two_photon_idx")
    """
    Attributes
    ----------
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
        if len(two_photon_data) != 1:
            raise NotImplementedError(
                "Currently only one two-photon decay is supported but there "
                f"are {len(two_photon_data)} in the atomic data."
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
    Attributes
    ----------
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


class LevelIdxs2TransitionIdx(HiddenPlasmaProperty):
    """
    Attributes
    ----------
    level_idxs2transition_idx : pandas.DataFrame, dtype int
       Maps a source_level_idx destination_level_idx pair to a transition_idx
       and transition type.
    """

    outputs = ("level_idxs2transition_idx",)

    def calculate(self, level_idxs2line_idx, level_idxs2continuum_idx):
        level_idxs2line_idx = level_idxs2line_idx.to_frame()
        level_idxs2line_idx.insert(1, "transition_type", -1)

        level_idxs2continuum_idx = level_idxs2continuum_idx.copy()
        level_idxs2continuum_idx.insert(1, "transition_type", -2)
        level_idxs2continuum_idx = level_idxs2continuum_idx.rename(
            columns=({"continuum_idx": "lines_idx"})
        )

        names = level_idxs2continuum_idx.index.names
        level_idxs2continuum_idx = level_idxs2continuum_idx.swaplevel()
        level_idxs2continuum_idx.index.names = names

        # TODO: This should probably be defined somewhere else.
        # One possibility would be to attach it to the cooling properties as
        # a class attribute.
        index_cooling = pd.MultiIndex.from_product(
            [["k"], ["ff", "adiabatic", "bf"]], names=names
        )
        num_cool = len(index_cooling)
        level_idxs2cooling_idx = pd.DataFrame(
            {
                "lines_idx": np.ones(num_cool, dtype=int) * -1,
                "transition_type": np.arange(-3, -3 - num_cool, -1),
            },
            index=index_cooling,
        )
        level_idxs2transition_idx = pd.concat(
            [
                level_idxs2continuum_idx,
                level_idxs2line_idx,
                level_idxs2cooling_idx,
            ]
        )

        # TODO: Add two-photon processes

        return level_idxs2transition_idx


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


class IsotopeMass(ProcessingPlasmaProperty):
    """
    Attributes
    ----------
    isotope_mass : pandas.Series
        Masses of the isotopes used. Indexed by isotope name e.g. 'Ni56'.
    """

    outputs = ("isotope_mass",)

    def calculate(self, isotope_abundance):
        """
        Determine mass of each isotope.

        Parameters
        ----------
        isotope_abundance : pandas.DataFrame
            Fractional abundance of isotopes.

        Returns
        -------
        pandas.DataFrame
            Masses of the isotopes used. Indexed by isotope name e.g. 'Ni56'.
        """
        if getattr(self, self.outputs[0]) is not None:
            return (getattr(self, self.outputs[0]),)
        else:
            if isotope_abundance.empty:
                return None
            isotope_mass_dict = {}
            for Z, A in isotope_abundance.index:
                element_name = rd.utils.Z_to_elem(Z)
                isotope_name = element_name + str(A)

                isotope_mass_dict[(Z, A)] = rd.Nuclide(isotope_name).atomic_mass

            isotope_mass_df = pd.DataFrame.from_dict(
                isotope_mass_dict, orient="index", columns=["mass"]
            )
            isotope_mass_df.index = pd.MultiIndex.from_tuples(
                isotope_mass_df.index
            )
            isotope_mass_df.index.names = ["atomic_number", "mass_number"]
            return isotope_mass_df / const.N_A


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
        if np.alltrue(keys + 1 == values) and keys:
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
        I_H = atomic_data.ionization_data.loc[(1, 1)]

        mask_selected_species = atomic_data.lines.index.droplevel(
            ["level_number_lower", "level_number_upper"]
        ).isin(continuum_interaction_species)
        lines_filtered = atomic_data.lines[mask_selected_species]
        f_lu = lines_filtered.f_lu.values
        nu_lines = lines_filtered.nu.values

        yg = f_lu * (I_H / (H * nu_lines)) ** 2
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
