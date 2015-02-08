from abc import ABCMeta

from tardis.plasma.plasma_properties import ProcessingPlasmaProperty
from tardis.plasma.exceptions import IncompleteAtomicData

class BaseAtomicDataProperty(ProcessingPlasmaProperty):
    __metaclass__ = ABCMeta

    inputs = ['atomic_data', 'selected_atoms']

    def __init__(self, plasma_parent):
        super(BaseAtomicDataProperty, self).__init__(plasma_parent)
        self.value = None

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
                    raw_atomic_property, selected_atoms))



class AtomicLevels(BaseAtomicDataProperty):
    name = 'levels'
    type_str = 'pandas.DataFrame'

    def _filter_atomic_property(self, levels, selected_atoms):
        return levels[levels.atomic_number.isin(selected_atoms)]

    def _set_index(self, levels):
        return levels.set_index(['atomic_number', 'ion_number', 'level_number'])

class AtomicLines(BaseAtomicDataProperty):
    name = 'lines'
    type_str = 'pandas.DataFrame'

    def _filter_atomic_property(self, lines, selected_atoms):
        return lines[lines.atomic_number.isin(selected_atoms)]

    def _set_index(self, lines):
        return lines

class AtomicMass(BaseAtomicDataProperty):
    name = 'atomic_mass'
    type_str = 'pandas.DataFrame'

    def calculate(self, atomic_data, selected_atoms):
        if self.value is not None:
            return self.value
        else:
            return atomic_data.atom_data.ix[selected_atoms].mass
