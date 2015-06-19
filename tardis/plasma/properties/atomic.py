from abc import ABCMeta, abstractmethod
import numpy as np
import pandas as pd

from tardis.plasma.properties.base import ProcessingPlasmaProperty
from tardis.plasma.exceptions import IncompleteAtomicData

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
            if not getattr(atomic_data, 'has_{0}'.format(
                    self.name)):
                raise IncompleteAtomicData(self.name)
            else:
                raw_atomic_property = getattr(atomic_data, '_' + self.name)
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

class IonizationData(ProcessingPlasmaProperty):
    name = 'ionization_data'

    def calculate(self, atomic_data, levels):
        if self.value is not None:
            return self.value
        else:
            return atomic_data.ionization_data

class ZetaData(ProcessingPlasmaProperty):
    name = 'zeta_data'

    def calculate(self, atomic_data, selected_atoms):
        if self.value is not None:
            return self.value
        else:
            return atomic_data.zeta_data
