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
    """
    Outputs:
    t_rad : Numpy Array
    """
    outputs = ('t_rad',)
    latex_name = r'$T_{\\textrm{rad}}$'

class DilutionFactor(ArrayInput):
    """
    Outputs:
    w : Numpy Array
        Factor used in nebular ionisation / dilute excitation calculations
        to account for the dilution of the radiation field.
    """
    outputs = ('w',)
    latex_name = r'$W$'
    latex_formula = r'$\frac{J_{\nu}}{B_{\nu}(T_\textrm{rad})}$'

class AtomicData(StaticInput):
    outputs = ('atomic_data',)

class Abundance(DynamicInput):
    outputs = ('abundance',)

class RadiationFieldCorrectionInput(StaticInput):
    """
    Outputs:
    delta_input : Numpy Array
        Used to adjust the ionisation balance to account for greater line
        blanketing in the blue.
    """
    outputs = ('delta_input',)
    latex_name = r'$\\delta$'

class Density(ArrayInput):
    outputs = ('density',)
    latex_name = r'$\\rho$'

class TimeExplosion(DynamicInput):
    outputs = ('time_explosion',)
    latex_name = r'$t_{\\textrm{exp}}$'

class JBlues(DataFrameInput):
    """
    Outputs:
    j_blues : Pandas DataFrame
        Mean intensity in the blue wing of each line.
    """
    outputs = ('j_blues',)
    latex_name = r'$J^{b}_{lu}$'

class LinkTRadTElectron(StaticInput):
    outputs = ('link_t_rad_t_electron',)
    latex_formula = r'$\\frac{T_{\\textrm{electron}}}{T_{\\textrm{rad}}}$'
