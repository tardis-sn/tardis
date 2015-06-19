import logging

import numpy as np
from astropy import constants as const

from tardis.plasma.properties.base import ProcessingPlasmaProperty

logger = logging.getLogger(__name__)

__all__ = ['BetaRadiation', 'GElectron', 'NumberDensity', 'SelectedAtoms',
           'ElectronTemperature', 'BetaElectron']

class BetaRadiation(ProcessingPlasmaProperty):
    name = 'beta_rad'
    latex_name = r'$\beta_\textrm{rad}$'
    latex_formula = r'$\frac{1}{K_B T_\textrm{rad}}$'

    def __init__(self, plasma_parent):
        super(BetaRadiation, self).__init__(plasma_parent)
        self.k_B_cgs = const.k_B.cgs.value


    def calculate(self, t_rad):
        return 1 / (self.k_B_cgs * t_rad)

class GElectron(ProcessingPlasmaProperty):

    name = 'g_electron'
    latex_name = r'$g_\textrm{electron}$'

    latex_formula = (r'$\left(\frac{2\pi m_\textrm{e} '
                     r'\beta_\textrm{rad}}{h^2}\right)^{3/2}$')

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
        return abundance.index

class ElectronTemperature(ProcessingPlasmaProperty):
    name = 't_electron'

    def calculate(self, t_rad, link_t_rad_t_electron):
        return t_rad * link_t_rad_t_electron

class BetaElectron(ProcessingPlasmaProperty):
    name = 'beta_electron'

    def __init__(self, plasma_parent):
        super(BetaElectron, self).__init__(plasma_parent)
        self.k_B_cgs = const.k_B.cgs.value

    def calculate(self, t_electron):
        return 1 / (self.k_B_cgs * t_electron)
