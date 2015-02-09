from abc import ABCMeta, abstractmethod, abstractproperty
import logging

from astropy import constants as const
import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)

class BasePlasmaProperty(object):
    __metaclass__ = ABCMeta

    @abstractproperty
    def name(self):
        pass

    def __init__(self):
        self.value = None


    def get_label(self):
        return "Name: {0}\nType: {1}\n{2}".format(self.name.replace('_', r'\\_'),
                                                  self.type_str,
                                                  getattr(self,
                                                          'latex_str', ''))
    def _update_type_str(self):
        self.type_str = repr(type(self.value))


    def get_latex_label(self):
        latex_template = r"""\textbf{{Name}} {name}
\textbf{{Formula}} {formula}
{description}
"""
        name = self.name.replace('_', r'\_')
        latex_name = getattr(self, 'latex_name', '')
        if latex_name != '':
            complete_name = '{0} [{1}]'.format(name, self.latex_name)
        else:
            complete_name = name


        latex_label = latex_template.format(name=complete_name,
                                     formula=getattr(self,
                                                     'latex_formula', '--'),
                                     description=getattr(self,
                                                         'latex_description',
                                                         ''))
        return latex_label.replace('\\', r'\\')

class ProcessingPlasmaProperty(BasePlasmaProperty):
    __metaclass__ = ABCMeta

    def __init__(self, plasma_parent):
        super(ProcessingPlasmaProperty, self).__init__()
        self.plasma_parent = plasma_parent
        self._update_inputs()
        self._update_inputs()

    def _update_inputs(self):
        """
        This function uses the CPython API to read the variable names from the
        `calculate`-function and makes the plasma routines easily programmable.
        """
        calculate_call_signature = self.calculate.func_code.co_varnames[
                                   :self.calculate.func_code.co_argcount]
        self.inputs = [item for item in calculate_call_signature if
                      item != 'self']


    def update(self):
        args = [getattr(self.plasma_parent, item) for item in self.inputs]
        self.value = self.calculate(*args)

    @abstractmethod
    def calculate(self, *args, **kwargs):
        raise NotImplementedError('This method needs to be implemented by ')


class BetaRadiation(ProcessingPlasmaProperty):

    name = 'beta_rad'
    latex_name = r'$\beta_\textrm{rad}$'

    latex_formula = r'$\frac{1}{K_B T_\textrm{rad}}$'

    def __init__(self, plasma_parent):
        super(BetaRadiation, self).__init__(plasma_parent)
        self.k_B_cgs = const.k_B.cgs.value


    def calculate(self, t_rad):
        return (1 / (self.k_B_cgs * t_rad))

class GElectron(ProcessingPlasmaProperty):

    name = 'g_electron'
    latex_name = r'$g_\textrm{electron}$'

    latex_formula = (r'\left(\\frac{2\pi m_\textrm{e} '
                     r'\beta_\textrm{rad}}{h^2}\right)^{3/2}')

    @staticmethod
    def calculate(beta_rad):
        return ((2 * np.pi * const.m_e.cgs.value / beta_rad) /
                (const.h.cgs.value ** 2)) ** 1.5


class NumberDensity(ProcessingPlasmaProperty):
    name = 'number_density'

    def calculate(self, atomic_mass, abundance, density):
        number_densities = (abundance * density)
        return number_densities.div(atomic_mass.ix[abundance.index], axis=0)

class SelectedAtoms(ProcessingPlasmaProperty):
    name = 'selected_atoms'

    def calculate(self, abundance):
        return self.plasma_parent.abundance.index

#### Importing properties from other modules ########
from tardis.plasma.partition_function import (LTEPartitionFunction,
                                              LevelBoltzmannFactor,
                                              DiluteLTEPartitionFunction)

from tardis.plasma.level_population import (LevelPopulationLTE,
                                            LevelNumberDensity)

from tardis.plasma.ion_population import (IonNumberDensity, PhiSahaLTE,
                                          PhiSahaNebular,
                                          RadiationFieldCorrection)
from tardis.plasma.radiative_properties import TauSobolev
from tardis.plasma.atomic_properties import (AtomicMass, AtomicLevels,
                                             AtomicLines, IonizationData)
######################################################
