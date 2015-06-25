from abc import ABCMeta, abstractmethod
import numpy as np
import pandas as pd
from collections import Counter as counter
import logging

from tardis.plasma.properties.base import ProcessingPlasmaProperty
from tardis.plasma.exceptions import IncompleteAtomicData

logger = logging.getLogger(__name__)

__all__ = ['Levels', 'Lines', 'LinesLowerLevelIndex', 'LinesUpperLevelIndex',
           'AtomicMass', 'IonizationData', 'ZetaData']

class BaseAtomicDataProperty(ProcessingPlasmaProperty):
    __metaclass__ = ABCMeta

    inputs = ['atomic_data', 'selected_atoms']

    def __init__(self, plasma_parent):
        super(BaseAtomicDataProperty, self).__init__(plasma_parent)
        self.value = None

    @abstractmethod
    def _set_index(self, raw_atomic_property, atomic_data):
        raise NotImplementedError('Needs to be implemented in subclasses')

    @abstractmethod
    def _filter_atomic_property(self, raw_atomic_property):
        raise NotImplementedError('Needs to be implemented in subclasses')


    def calculate(self, atomic_data, selected_atoms):
        if self.value is not None:
            return self.value
        else:
            try:
                raw_atomic_property = getattr(atomic_data, '_' + self.name)
                return self._set_index(self._filter_atomic_property(
                    raw_atomic_property, selected_atoms), atomic_data)
            except:
                raw_atomic_property = getattr(atomic_data, self.name)
                return self._set_index(self._filter_atomic_property(
                    raw_atomic_property, selected_atoms), atomic_data)


class Levels(BaseAtomicDataProperty):
    name = 'levels'

    def _filter_atomic_property(self, levels, selected_atoms):
        return levels[levels.atomic_number.isin(selected_atoms)]

    def _set_index(self, levels, atomic_data):
        return levels.set_index(['atomic_number', 'ion_number',
            'level_number'])

class Lines(BaseAtomicDataProperty):
    name = 'lines'

    def _filter_atomic_property(self, lines, selected_atoms):
        return lines[lines.atomic_number.isin(selected_atoms)]

    def _set_index(self, lines, atomic_data):
        try:
            reindexed = lines.reindex(atomic_data.lines.index)
        except:
            reindexed = lines.reindex(atomic_data._lines.index)
        return reindexed

class LinesLowerLevelIndex(ProcessingPlasmaProperty):
    name = 'lines_lower_level_index'

    def calculate(self, levels, lines):
        levels_index = pd.Series(np.arange(len(levels), dtype=np.int64),
                                 index=levels.index)
        lines_index = lines.set_index(
            ['atomic_number', 'ion_number',
             'level_number_lower']).index
        return np.array(levels_index.ix[lines_index])

class LinesUpperLevelIndex(ProcessingPlasmaProperty):
    name = 'lines_upper_level_index'

    def calculate(self, levels, lines):
        levels_index = pd.Series(np.arange(len(levels), dtype=np.int64),
                                 index=levels.index)
        lines_index = lines.set_index(
            ['atomic_number', 'ion_number',
             'level_number_upper']).index
        return np.array(levels_index.ix[lines_index])


class IonCXData(BaseAtomicDataProperty):
    name = 'ion_cx_data'

    def _filter_atomic_property(self, ion_cx_data, selected_atoms):
        return filtered_ion_cx_data

    def _set_index(self, ion_cx_data, atomic_data):
        return levels.set_index(['atomic_number', 'ion_number',
                                 'level_number'])


class AtomicMass(ProcessingPlasmaProperty):
    name = 'atomic_mass'

    def calculate(self, atomic_data, selected_atoms):
        if self.value is not None:
            return self.value
        else:
            return atomic_data.atom_data.ix[selected_atoms].mass

class IonizationData(BaseAtomicDataProperty):
    name = 'ionization_data'

    def _filter_atomic_property(self, ionization_data, selected_atoms):
        ionization_data['atomic_number'] = ionization_data.index.labels[0]+1
        ionization_data['ion_number'] = ionization_data.index.labels[1]+1
        ionization_data = ionization_data[ionization_data.atomic_number.isin(
            selected_atoms)]
        ion_data_check = counter(ionization_data.atomic_number.values)
        keys = np.array(ion_data_check.keys())
        values = np.array(ion_data_check.values())
        if np.alltrue(keys==values):
            return ionization_data
        else:
            raise IncompleteAtomicData('ionization data for the ion (' +
                            str(keys[keys!=values]) +
                            str(values[keys!=values]) + ')')

    def _set_index(self, ionization_data, atomic_data):
        return ionization_data.set_index(['atomic_number', 'ion_number'])

class ZetaData(BaseAtomicDataProperty):
    name = 'zeta_data'

    def _filter_atomic_property(self, zeta_data, selected_atoms):
        zeta_data['atomic_number'] = zeta_data.index.labels[0]+1
        zeta_data['ion_number'] = zeta_data.index.labels[1]+1
        zeta_data =  zeta_data[zeta_data.atomic_number.isin(selected_atoms)]
        zeta_data_check = counter(zeta_data.atomic_number.values)
        keys = np.array(zeta_data_check.keys())
        values = np.array(zeta_data_check.values())
        if np.alltrue(keys+1==values):
            return zeta_data
        else:
            logger.warn('Zeta_data missing - replaced with 1s')
            updated_index = []
            for atom in selected_atoms:
                for ion in range(1, atom+2):
                    updated_index.append([atom,ion])
            updated_index = np.array(updated_index)
            updated_dataframe = pd.DataFrame(index=pd.MultiIndex.from_arrays(
                updated_index.transpose().astype(int)),
                columns = zeta_data.columns)
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

    def _set_index(self, zeta_data, atomic_data):
        return zeta_data.set_index(['atomic_number', 'ion_number'])

