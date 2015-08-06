import logging

import numpy as np
from astropy import constants as const

from tardis.plasma.properties.base import ProcessingPlasmaProperty

logger = logging.getLogger(__name__)

__all__ = ['BetaRadiation', 'GElectron', 'NumberDensity', 'SelectedAtoms',
           'ElectronTemperature', 'BetaElectron']

class BetaRadiation(ProcessingPlasmaProperty):
    """
    Outputs:
    beta_rad : Numpy Array
    """
    outputs = ('beta_rad',)
    latex_name = ('\\beta_{\\textrm{rad}}',)
    latex_formula = ('\\dfrac{1}{k_{B} T_{\\textrm{rad}}}',)

    def __init__(self, plasma_parent):
        super(BetaRadiation, self).__init__(plasma_parent)
        self.k_B_cgs = const.k_B.cgs.value

    def calculate(self, t_rad):
        return 1 / (self.k_B_cgs * t_rad)

class GElectron(ProcessingPlasmaProperty):
    """
    Outputs:
    g_electron : Numpy Array
    """
    outputs = ('g_electron',)
    latex_name = ('g_{\\textrm{electron}}',)
    latex_formula = ('\\Big(\\dfrac{2\\pi m_{e}/\
                     \\beta_{\\textrm{rad}}}{h^2}\\Big)^{3/2}',)

    def calculate(self, beta_rad):
        return ((2 * np.pi * const.m_e.cgs.value / beta_rad) /
                (const.h.cgs.value ** 2)) ** 1.5

class NumberDensity(ProcessingPlasmaProperty):
    """
    Outputs:
    number_density : Pandas DataFrame
    """
    outputs = ('number_density',)
    latex_name = ('N_{i}',)

    def calculate(self, atomic_mass, abundance, density):
        number_densities = (abundance * density)
        return number_densities.div(atomic_mass.ix[abundance.index], axis=0)

class SelectedAtoms(ProcessingPlasmaProperty):
    """
    Outputs:
    selected_atoms : Numpy Array
        Elements required for particular simulation
    """
    outputs = ('selected_atoms',)

    def calculate(self, abundance):
        return abundance.index

class ElectronTemperature(ProcessingPlasmaProperty):
    """
    Outputs:
    t_electron : Numpy Array
    """
    outputs = ('t_electrons',)
    latex_name = ('T_{\\textrm{electron}}',)
    latex_formula = ('\\textrm{const.}\\times T_{\\textrm{rad}}',)

    def calculate(self, t_rad, link_t_rad_t_electron):
        return t_rad * link_t_rad_t_electron

class BetaElectron(ProcessingPlasmaProperty):
    """
    Outputs:
    beta_electron : Numpy Array
    """
    outputs = ('beta_electron',)
    latex_name = ('\\beta_{\\textrm{electron}}',)
    latex_formula = ('\\frac{1}{K_{B} T_{\\textrm{electron}}}',)

    def __init__(self, plasma_parent):
        super(BetaElectron, self).__init__(plasma_parent)
        self.k_B_cgs = const.k_B.cgs.value

    def calculate(self, t_electron):
        return 1 / (self.k_B_cgs * t_electron)
