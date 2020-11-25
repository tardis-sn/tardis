import logging

import numpy as np
import pandas as pd
from collections import Counter as counter

from tardis.plasma.properties.base import (
    ProcessingPlasmaProperty,
    HiddenPlasmaProperty,
    BaseAtomicDataProperty,
)
from tardis.plasma.exceptions import IncompleteAtomicData

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
]


class Levels(BaseAtomicDataProperty):
    """
    Attributes
    ----------
    levels : Pandas MultiIndex (atomic_number, ion_number, level_number)
             Index of filtered atomic data. Index used for all other attribute dataframes for this class
    excitation_energy : Pandas DataFrame (index=levels), dtype float
             Excitation energies of atomic levels
    metastability : Pandas DataFrame (index=levels), dtype bool
             Records whether atomic levels are metastable
    g : Pandas DataFrame (index=levels), dtype float
             Statistical weights of atomic levels
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
    lines : Pandas DataFrame (wavelength, atomic_number, ion_number, f_ul, f_lu, level_number_lower,
                              level_number_upper, nu, B_lu, B_ul, A_ul, wavelength)
            All atomic lines data. Index = line_id.
    nu : Pandas DataFrame (index=line_id), dtype float
            Line frequency data
    f_lu : Pandas DataFrame (index=line_id), dtype float
            Transition probability data
    wavelength_cm: Pandas DataFrame (index=line_id), dtype float
            Line wavelengths in cm
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
    photo_ion_cross_sections: Pandas DataFrame (nu, x_sect,
                                                index=['atomic_number',
                                                       'ion_number',
                                                       'level_number']),
                                                dtype float)
                              Table of photoionization cross sections as a
                              function of frequency.
    photo_ion_block_references: One-dimensional Numpy Array, dtype int
                              Indices where the photoionization data for
                              a given level starts. Needed for calculation
                              of recombination rates.
    photo_ion_index: Pandas MultiIndex, dtype int
                              Atomic, ion and level numbers for which
                              photoionization data exists.
    """

    outputs = (
        "photo_ion_cross_sections",
        "photo_ion_block_references",
        "photo_ion_index",
    )
    latex_name = (r"\xi_{\textrm{i}}(\nu)", "", "")

    def calculate(self, atomic_data, continuum_interaction_species):
        photoionization_data = atomic_data.photoionization_data.set_index(
            ["atomic_number", "ion_number", "level_number"]
        )
        selected_species_idx = pd.IndexSlice[
            continuum_interaction_species.get_level_values("atomic_number"),
            continuum_interaction_species.get_level_values("ion_number"),
            slice(None),
        ]
        photoionization_data = photoionization_data.loc[selected_species_idx]
        phot_nus = photoionization_data["nu"]
        block_references = np.hstack(
            [[0], phot_nus.groupby(level=[0, 1, 2]).count().values.cumsum()]
        )
        photo_ion_index = photoionization_data.index.unique()
        return photoionization_data, block_references, photo_ion_index


class LinesLowerLevelIndex(HiddenPlasmaProperty):
    """
    Attributes:
    lines_lower_level_index : One-dimensional Numpy Array, dtype int
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
    Attributes:
    lines_upper_level_index : One-dimensional Numpy Array, dtype int
        Levels data for upper levels of particular lines
    """

    outputs = ("lines_upper_level_index",)

    def calculate(self, levels, lines):
        levels_index = pd.Series(
            np.arange(len(levels), dtype=np.int64), index=levels
        )
        lines_index = lines.index.droplevel("level_number_lower")
        return np.array(levels_index.loc[lines_index])


class AtomicMass(ProcessingPlasmaProperty):
    """
    Attributes:
    atomic_mass : Pandas Series
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
    Attributes:
    ionization_data : Pandas Series holding ionization energies
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
                "ionization data for the ion ({}, {})".format(
                    str(counts.index[counts.index != counts]),
                    str(counts[counts.index != counts]),
                )
            )

    def _set_index(self, ionization_data):
        return ionization_data


class ZetaData(BaseAtomicDataProperty):
    """
    Attributes:
    zeta_data : Pandas DataFrame, dtype float
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
                "Zeta_data missing - replaced with 1s. Missing ions: {}".format(
                    missing_ions
                )
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
        Attributes:
        nlte_data :
    #Finish later (need atomic dataset with NLTE data).
    """

    outputs = ("nlte_data",)

    def calculate(self, atomic_data):
        if getattr(self, self.outputs[0]) is not None:
            return (getattr(self, self.outputs[0]),)
        else:
            return atomic_data.nlte_data
