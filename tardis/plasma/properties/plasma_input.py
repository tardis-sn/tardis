import numpy as np
import pandas as pd
from tardis.plasma.properties.base import BasePlasmaProperty

__all__ = ['TRadiative', 'DilutionFactor', 'AtomicData', 'Abundance', 'Density',
           'TimeExplosion', 'JBlues', 'LinkTRadTElectron',
           'RadiationFieldCorrectionInput']


class Input(BasePlasmaProperty):

    def _set_output_value(self, output, value):
        setattr(self, output, value)

    def set_value(self, value):
        assert len(self.outputs) == 1
        self._set_output_value(self.outputs[0], value)


class StaticInput(Input):
    pass

class DynamicInput(Input):
    pass


class ArrayInput(DynamicInput):
    def _set_output_value(self, output, value):
        setattr(self, output, np.array(value, copy=False))


class DataFrameInput(DynamicInput):
    def _set_output_value(self, output, value):
        setattr(self, output, np.array(pd.DataFrame(value), copy=False))

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
