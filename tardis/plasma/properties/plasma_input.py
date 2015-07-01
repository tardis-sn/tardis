import numpy as np
import pandas as pd
from tardis.plasma.properties.base import BasePlasmaProperty

__all__ = ['TRadiative', 'DilutionFactor', 'AtomicData', 'Abundance', 'Density',
           'TimeExplosion', 'JBlues', 'LinkTRadTElectron',
           'RadiationFieldCorrectionInput']


class Input(BasePlasmaProperty):

    def set_value(self, value):
        self.value = value


class StaticInput(Input):
    pass

class DynamicInput(Input):
    pass


class ArrayInput(DynamicInput):
    def set_value(self, value):
        self.value = np.array(value, copy=False)

class DataFrameInput(DynamicInput):
    def set_value(self, value):
        self.value = pd.DataFrame(value)

class TRadiative(ArrayInput):
    outputs = ('t_rad',)
    latex_name = r'$T_\textrm{rad}$'


class DilutionFactor(ArrayInput):
    outputs = ('w',)
    latex_name = r'$W$'

class AtomicData(StaticInput):
    outputs = ('atomic_data',)


class Abundance(DynamicInput):
    outputs = ('abundance',)


class RadiationFieldCorrectionInput(StaticInput):
    outputs = ('delta_input',)

class Density(ArrayInput):
    outputs = ('density',)
    latex_name = r'$\rho$'

class TimeExplosion(DynamicInput):
    outputs = ('time_explosion',)

class JBlues(DataFrameInput):
    outputs = ('j_blues',)

class LinkTRadTElectron(StaticInput):
    outputs = ('link_t_rad_t_electron',)
