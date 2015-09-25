import logging

import numpy as np
import pandas as pd
from collections import Counter as counter

from tardis.plasma.properties.base import (ProcessingPlasmaProperty,
    HiddenPlasmaProperty, BaseAtomicDataProperty)
from tardis.plasma.exceptions import IncompleteAtomicData

logger = logging.getLogger(__name__)

__all__ = ['Levels', 'Lines', 'LinesLowerLevelIndex', 'LinesUpperLevelIndex',
           'AtomicMass', 'IonizationData', 'ZetaData', 'NLTEData', 'Chi0']

class Levels(BaseAtomicDataProperty):
    """
    Outputs:
    levels : Pandas DataFrame
        Levels data needed for particular simulation
    """
    outputs = ('levels', 'excitation_energy', 'metastability', 'g')
    latex_name = ('\\textrm{levels}', '\\epsilon_{\\textrm{k}}', '\\textrm{metastability}',
        'g')

    def _filter_atomic_property(self, levels, selected_atoms):
        return levels[levels.atomic_number.isin(selected_atoms)]

    def _set_index(self, levels):
        levels = levels.set_index(['atomic_number', 'ion_number',
                                 'level_number'])
        return (levels.index, levels['energy'], levels['metastable'],
            levels['g'])

class Lines(BaseAtomicDataProperty):
    """
    Outputs:
    lines : Pandas DataFrame
        Lines data needed for particular simulation
    """
    outputs = ('lines', 'nu', 'f_lu', 'wavelength_cm')

    def _filter_atomic_property(self, lines, selected_atoms):
        return lines[lines.atomic_number.isin(selected_atoms)]

    def _set_index(self, lines):
# Filtering process re-arranges the index. This re-orders it.
# It seems to be important that the lines stay indexed in the correct order
# so that the tau_sobolevs values are in the right order for the montecarlo
# code, or it returns the wrong answer.
        try:
            reindexed = lines.reindex(lines.index)
        except:
            reindexed = lines.reindex(lines.index)
        lines = reindexed.dropna(subset=['atomic_number'])
        return lines, lines['nu'], lines['f_lu'], lines['wavelength_cm']

class LinesLowerLevelIndex(HiddenPlasmaProperty):
    """
    Outputs:
    lines_lower_level_index : One-dimensional Numpy Array
        Levels data for lower levels of particular lines
        Usage: levels.ix[lines_lower_level_index]
    """
    outputs = ('lines_lower_level_index',)
    def calculate(self, levels, lines):
        levels_index = pd.Series(np.arange(len(levels), dtype=np.int64),
                                 index=levels)
        lines_index = lines.set_index(
            ['atomic_number', 'ion_number',
             'level_number_lower']).index
        return np.array(levels_index.ix[lines_index])

class LinesUpperLevelIndex(HiddenPlasmaProperty):
    """
    Outputs:
    lines_upper_level_index : One-dimensional Numpy Array
        Levels data for upper levels of particular lines
        Usage: levels.ix[lines_upper_level_index]
    """
    outputs = ('lines_upper_level_index',)

    def calculate(self, levels, lines):
        levels_index = pd.Series(np.arange(len(levels), dtype=np.int64),
                                 index=levels)
        lines_index = lines.set_index(
            ['atomic_number', 'ion_number',
             'level_number_upper']).index
        return np.array(levels_index.ix[lines_index])

class IonCXData(BaseAtomicDataProperty):
    outputs = ('ion_cx_data',)

    def _filter_atomic_property(self, ion_cx_data, selected_atoms):
        return ion_cx_data[ion_cx_data.atomic_number.isin([selected_atoms]
                                                          if np.isscalar(
            selected_atoms)
                                                          else selected_atoms)]

    def _set_index(self, ion_cx_data):
        return ion_cx_data.set_index(['atomic_number', 'ion_number',
                                      'level_number'])

class AtomicMass(ProcessingPlasmaProperty):
    """
    Outputs:
    atomic_mass : Pandas Series
        Atomic masses of the elements used, indexed by atomic number
    """
    outputs = ('atomic_mass',)

    def calculate(self, atomic_data, selected_atoms):
        if getattr(self, self.outputs[0]) is not None:
            return getattr(self, self.outputs[0]),
        else:
            return atomic_data.atom_data.ix[selected_atoms].mass

class IonizationData(BaseAtomicDataProperty):
    """
    Outputs:
    ionization_data : Pandas DataFrame
        Ionization energies of the elements used
    """
    outputs = ('ionization_data',)

    def _filter_atomic_property(self, ionization_data, selected_atoms):
        ionization_data['atomic_number'] = ionization_data.index.labels[0] + 1
        ionization_data['ion_number'] = ionization_data.index.labels[1] + 1
        ionization_data = ionization_data[ionization_data.atomic_number.isin(
            selected_atoms)]
        ion_data_check = counter(ionization_data.atomic_number.values)
        keys = np.array(ion_data_check.keys())
        values = np.array(ion_data_check.values())
        if np.alltrue(keys == values):
            return ionization_data
        else:
            raise IncompleteAtomicData('ionization data for the ion (' +
                                       str(keys[keys != values]) +
                                       str(values[keys != values]) + ')')

    def _set_index(self, ionization_data):
        return ionization_data.set_index(['atomic_number', 'ion_number'])

class ZetaData(BaseAtomicDataProperty):
    """
    Outputs:
    zeta_data : Pandas DataFrame
        Zeta data for the elements used
        Required for the nebular ionization scheme.
        The zeta value represents the fraction of recombination events
        from the ionized state that go directly to the ground state.
    """
    outputs = ('zeta_data',)

    def _filter_atomic_property(self, zeta_data, selected_atoms):
        for temperature in range(42000, 82000, 2000):
            zeta_data[temperature] = zeta_data[40000]
        zeta_data['atomic_number'] = zeta_data.index.labels[0] + 1
        zeta_data['ion_number'] = zeta_data.index.labels[1] + 1
        zeta_data = zeta_data[zeta_data.atomic_number.isin(selected_atoms)]
        zeta_data_check = counter(zeta_data.atomic_number.values)
        keys = np.array(zeta_data_check.keys())
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
                        missing_ions.append((atom,ion))
                    updated_index.append([atom, ion])
            logger.warn('Zeta_data missing - replaced with 1s. Missing ions: {}'.format(missing_ions))
            updated_index = np.array(updated_index)
            updated_dataframe = pd.DataFrame(index=pd.MultiIndex.from_arrays(
                updated_index.transpose().astype(int)),
                columns=zeta_data.columns)
            for value in range(len(zeta_data)):
                updated_dataframe.ix[zeta_data.atomic_number.values[value]].ix[
                    zeta_data.ion_number.values[value]] = \
                    zeta_data.ix[zeta_data.atomic_number.values[value]].ix[
                        zeta_data.ion_number.values[value]]
            updated_dataframe = updated_dataframe.astype(float)
            updated_index = pd.DataFrame(updated_index)
            updated_dataframe['atomic_number'] = np.array(updated_index[0])
            updated_dataframe['ion_number'] = np.array(updated_index[1])
            updated_dataframe.fillna(1.0, inplace=True)
            return updated_dataframe

    def _set_index(self, zeta_data):
        return zeta_data.set_index(['atomic_number', 'ion_number'])

class NLTEData(ProcessingPlasmaProperty):
    outputs = ('nlte_data',)

    def calculate(self, atomic_data):
        if getattr(self, self.outputs[0]) is not None:
            return (getattr(self, self.outputs[0]),)
        else:
            return atomic_data.nlte_data

class Chi0(ProcessingPlasmaProperty):
    outputs = ('chi_0',)

    def calculate(self, atomic_data):
        return atomic_data.ionization_data.ionization_energy.ix[20].ix[2]
