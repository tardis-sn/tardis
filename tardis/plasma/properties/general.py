import logging

import numpy as np
import pandas as pd
from astropy import units as u
from tardis import constants as const

from tardis.plasma.properties.base import ProcessingPlasmaProperty

logger = logging.getLogger(__name__)

__all__ = ['BetaRadiation', 'GElectron', 'NumberDensity', 'SelectedAtoms',
           'ElectronTemperature', 'BetaElectron', 'LuminosityInner',
           'TimeSimulation', 'ThermalGElectron']

class BetaRadiation(ProcessingPlasmaProperty):
    """
    Attributes
    ----------
    beta_rad : Numpy Array, dtype float
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
    Attributes
    ----------
    g_electron : Numpy Array, dtype float
    """
    outputs = ('g_electron',)
    latex_name = ('g_{\\textrm{electron}}',)
    latex_formula = ('\\Big(\\dfrac{2\\pi m_{e}/\
                     \\beta_{\\textrm{rad}}}{h^2}\\Big)^{3/2}',)

    def calculate(self, beta_rad):
        return ((2 * np.pi * const.m_e.cgs.value / beta_rad) /
                (const.h.cgs.value ** 2)) ** 1.5


class ThermalGElectron(GElectron):
    """
    Attributes
    ----------
    thermal_g_electron : Numpy Array, dtype float
    """
    outputs = ('thermal_g_electron',)
    latex_name = ('g_{\\textrm{electron_thermal}}',)
    latex_formula = ('\\Big(\\dfrac{2\\pi m_{e}/\
                     \\beta_{\\textrm{electron}}}{h^2}\\Big)^{3/2}',)

    def calculate(self, beta_electron):
        return super(ThermalGElectron, self).calculate(beta_electron)


class NumberDensity(ProcessingPlasmaProperty):
    """
    Attributes
    ----------
    number_density : Pandas DataFrame, dtype float
                     Indexed by atomic number, columns corresponding to zones
    """
    outputs = ('number_density',)
    latex_name = ('N_{i}',)

    @staticmethod
    def calculate(atomic_mass, abundance, density):
        number_densities = (abundance * density)
        return number_densities.div(atomic_mass.loc[abundance.index], axis=0)

class SelectedAtoms(ProcessingPlasmaProperty):
    """
    Attributes
    ----------
    selected_atoms : Pandas Int64Index, dtype int
                     Atomic numbers of elements required for particular simulation
    """
    outputs = ('selected_atoms',)

    def calculate(self, abundance):
        return abundance.index

class ElectronTemperature(ProcessingPlasmaProperty):
    """
    Attributes
    ----------
    t_electron : Numpy Array, dtype float
    """
    outputs = ('t_electrons',)
    latex_name = ('T_{\\textrm{electron}}',)
    latex_formula = ('\\textrm{const.}\\times T_{\\textrm{rad}}',)

    def calculate(self, t_rad, link_t_rad_t_electron):
        return t_rad * link_t_rad_t_electron

class BetaElectron(ProcessingPlasmaProperty):
    """
    Attributes
    ----------
    beta_electron : Numpy Array, dtype float
    """
    outputs = ('beta_electron',)
    latex_name = ('\\beta_{\\textrm{electron}}',)
    latex_formula = ('\\frac{1}{K_{B} T_{\\textrm{electron}}}',)

    def __init__(self, plasma_parent):
        super(BetaElectron, self).__init__(plasma_parent)
        self.k_B_cgs = const.k_B.cgs.value

    def calculate(self, t_electrons):
        return 1 / (self.k_B_cgs * t_electrons)

class LuminosityInner(ProcessingPlasmaProperty):
    outputs = ('luminosity_inner',)

    @staticmethod
    def calculate(r_inner, t_inner):
        return (4 * np.pi * const.sigma_sb.cgs * r_inner[0] ** 2
                * t_inner ** 4).to('erg/s')

class TimeSimulation(ProcessingPlasmaProperty):
    outputs = ('time_simulation',)

    @staticmethod
    def calculate(luminosity_inner):
        return 1.0 * u.erg / luminosity_inner
