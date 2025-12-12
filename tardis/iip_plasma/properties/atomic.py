import logging
from collections import Counter as counter

import numpy as np
import pandas as pd

from tardis.iip_plasma.exceptions import IncompleteAtomicData
from tardis.iip_plasma.properties.base import (
    BaseAtomicDataProperty,
    HiddenPlasmaProperty,
    ProcessingPlasmaProperty,
)

logger = logging.getLogger(__name__)

__all__ = [
    "AtomicMass",
    "ContinuumData",
    "IonizationData",
    "Levels",
    "Lines",
    "LinesLowerLevelIndex",
    "LinesUpperLevelIndex",
    "NLTEData",
    "YgData",
    "ZetaData",
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
        "\\textrm{levels}",
        "\\epsilon_{\\textrm{k}}",
        "\\textrm{metastability}",
        "g",
    )

    def _filter_atomic_property(self, levels, selected_atoms):
        return levels.loc[(selected_atoms), :]

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

    # Would like for lines to just be the line_id values, but this is necessary because of the filtering
    # error noted below.
    outputs = ("lines", "nu", "f_lu", "wavelength_cm", "lines_multi_index")

    def _filter_atomic_property(self, lines, selected_atoms):
        # return lines[lines.atomic_number.isin(selected_atoms)]
        return lines

    def _set_index(self, lines):
        # Filtering process re-arranges the index. This re-orders it.
        # It seems to be important that the lines stay indexed in the correct order
        # so that the tau_sobolevs values are in the right order for the montecarlo
        # code, or it returns the wrong answer.
        # try:
        #     reindexed = lines.reindex(lines.index)
        # except:
        #     reindexed = lines.reindex(lines.index)
        # lines = reindexed.dropna(subset=['atomic_number'])
        return (
            lines,
            lines["nu"],
            lines["f_lu"],
            lines["wavelength_cm"],
            lines.index,
        )


class LinesLowerLevelIndex(HiddenPlasmaProperty):
    """
    Attributes
    ----------
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
    Attributes
    ----------
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


class IonCXData(BaseAtomicDataProperty):
    outputs = ("ion_cx_data",)

    def _filter_atomic_property(self, ion_cx_data, selected_atoms):
        return ion_cx_data[
            ion_cx_data.atomic_number.isin(
                [selected_atoms]
                if np.isscalar(selected_atoms)
                else selected_atoms
            )
        ]

    def _set_index(self, ion_cx_data):
        return ion_cx_data.set_index(
            ["atomic_number", "ion_number", "level_number"]
        )


class AtomicMass(ProcessingPlasmaProperty):
    """
    Attributes
    ----------
    atomic_mass : Pandas Series
        Atomic masses of the elements used. Indexed by atomic number.
    """

    outputs = ("atomic_mass",)

    def calculate(self, atomic_data, selected_atoms):
        if getattr(self, self.outputs[0]) is not None:
            return (getattr(self, self.outputs[0]),)
        return atomic_data.atom_data.loc[selected_atoms].mass


class IonizationData(BaseAtomicDataProperty):
    """
    Attributes
    ----------
    ionization_data : Pandas DataFrame
        Ionization energies. Indexed by atomic number, ion number.
    """

    outputs = ("ionization_data",)

    def _filter_atomic_property(self, ionization_data, selected_atoms):
        mask = ionization_data.index.isin(selected_atoms, level="atomic_number")
        ionization_data = ionization_data[mask]
        counts = ionization_data.groupby(level="atomic_number").count()

        if np.all(counts.index == counts):
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
        if np.all(keys + 1 == values) and keys:
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
            logger.warning(
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
            updated_dataframe = updated_dataframe.fillna(1.0)
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
        return atomic_data.nlte_data


class ContinuumData(BaseAtomicDataProperty):
    """
    Attributes
    ----------
    continuum_data : `tardis.iip_plasma.continuum.base_continuum_data.ContinuumData`-object
             Container object for continuum data.
    photo_ion_cross_sections : Pandas DataFrame (nu, x_sect, index=['atomic_number',
                                                'ion_number','level_number']), dtype float
             Table of photoionization cross sections as a function of frequency.
    photo_ion_index_sorted : Pandas MultiIndex (atomic_number, ion_number, level_number)
             Index of levels sorted with respect to the threshold frequency of photoionization (decreasing).
    """

    outputs = (
        "continuum_data",
        "photo_ion_cross_sections",
        "photo_ion_index_sorted",
    )
    latex_name = ("\\omega_{\\textrm{i}}", "\\textrm{photo_ion_index_sorted}")

    def _filter_atomic_property(self, continuum_data, selected_atoms):
        return continuum_data

    def _set_index(self, continuum_data):
        return (
            continuum_data,
            continuum_data.photoionization_data,
            continuum_data.multi_index_nu_sorted,
        )


class YgData(BaseAtomicDataProperty):
    """
    Attributes
    ----------
    yg_data : Pandas DataFrame
        Table of effective collision strengths (divided by g_i). Columns are temperatures.
    T_Yg: Numpy Array
        Temperatures at which collision strengths are tabulated.
    Yg_index: Pandas MultiIndex
    """

    outputs = ("yg_data", "T_Yg", "Yg_index")
    latex_name = ("Y__{\\textrm{ij} / g_i}", "T_\\textrm{Yg}", "Yg_index")

    def _filter_atomic_property(self, yg_data, selected_atoms):
        return yg_data

    def _set_index(self, yg_data):
        return (yg_data, yg_data.columns.values.astype(float), yg_data.index)
