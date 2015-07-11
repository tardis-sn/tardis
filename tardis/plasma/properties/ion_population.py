import logging

import numpy as np
import pandas as pd
from scipy import interpolate

from tardis.plasma.properties.base import ProcessingPlasmaProperty
from tardis.plasma.exceptions import PlasmaIonizationError

logger = logging.getLogger(__name__)

__all__ = ['PhiSahaNebular', 'PhiSahaLTE', 'RadiationFieldCorrection',
           'IonNumberDensity', 'PhiGeneral']

class PhiGeneral(ProcessingPlasmaProperty):
    """
    Outputs:
    general_phi : Pandas DataFrame
        Used as basis for PhiSahaLTE and PhiSahaNebular. Identical output to
        PhiSahaLTE. Separate property required as PhiSahaNebular is based
        on PhiSahaLTE, but the code cannot deal with the inclusion of two
        properties that generate a property called 'phi'.
    """
    outputs = ('general_phi',)
    latex_name = r'$\Phi_{i,j}(LTE)$'
    latex_formula = (r'$\frac{N_{i, j+1} n_e}{N_{i, j}} \\'
                     r' \Phi_{i, j} = g_e \times \frac{Z_{i, j+1}}{Z_{i, j}} '
                     r'e^{-\chi_{j\rightarrow j+1}/k_\textrm{B}T}$')

    @staticmethod
    def calculate(g_electron, beta_rad, partition_function, ionization_data):
        def calculate_phis(group):
            return group[1:] / group[:-1].values

        phis = partition_function.groupby(level='atomic_number').apply(
            calculate_phis)
        phis = pd.DataFrame(phis.values, index=phis.index.droplevel(0))
        phi_coefficient = (2 * g_electron * np.exp(np.outer(
            ionization_data.ionization_energy.ix[phis.index].values,
            -beta_rad)))
        return phis * phi_coefficient

class PhiSahaNebular(ProcessingPlasmaProperty):
    """
    Outputs:
    phi_saha_nebular: Pandas DataFrame
        The ionization equilibrium as calculated using a modified version of
        the Saha equation that accounts for dilution of the radiation field.
    """
    outputs = ('phi',)
    latex_name = r'$\Phi_{i,j}$'
    latex_formula = r'$\delta W[\zeta+W(1-\zeta)]\Phi_{i,j}_{\textrm{(LTE)}}$'

    def calculate(self, general_phi, t_rad, w, zeta_data, t_electron, delta):
        try:
            zeta = interpolate.interp1d(zeta_data.columns.values, zeta_data.ix[
                general_phi.index].values)(t_rad)
            zeta = zeta.astype(float)
        except ValueError:
            raise ValueError('t_rads outside of zeta factor interpolation'
                             ' zeta_min={0:.2f} zeta_max={1:.2f} '
                             '- requested {2}'.format(
                zeta_data.columns.values.min(), zeta_data.columns.values.max(),
                t_rad))
        phis = general_phi * delta * w * (zeta + w * (1 - zeta)) * \
               (t_electron/t_rad) ** .5
        return phis

class PhiSahaLTE(ProcessingPlasmaProperty):
    """
    Outputs:
    phi_saha_lte: Pandas DataFrame
        The ionization equilibrium as calculated using the Saha equation.
    """
    outputs = ('phi',)
    latex_name = r'$\Phi_{i,j}$'
    latex_formula = r'$\Phi_{i,j}_{\textrm{(LTE)}}$'

    def calculate(self, general_phi):
        return general_phi

class RadiationFieldCorrection(ProcessingPlasmaProperty):
    """
    Outputs:
    delta: Pandas DataFrame
        Calculates the radiation field correction (see Mazzali & Lucy, 1993) if
        not given as input in the config. file. The default chi_0_species is
        Ca II, which is good for type Ia supernovae. For type II supernovae,
        (1, 1) should be used.
    """
    outputs = ('delta',)
    latex_name = r'$\delta$'
    latex_formula = r'\chi>\chi_0:\frac{T_{e}}{b_1 W T_{R}}\exp\left ('\
                    r'\frac{\chi}{kT_{R}} - \frac{\chi}{kT_{e}} \right )'\
                    r'\chi<\chi_0: 1 - \exp\left ( \frac{h\nu}{kT_R} -'\
                    r'\frac{h\nu_0}{kT_R} \right ) + \frac{T_e}{b_1 WT_R}'\
                    r'\exp\left ( \frac{h\nu}{kT_R} -\frac{h\nu_0}{kT_e}'\
                    r'\right )'

    def __init__(self, plasma_parent, departure_coefficient=None,
        chi_0_species=(20,2)):
        super(RadiationFieldCorrection, self).__init__(plasma_parent)
        self.departure_coefficient = departure_coefficient
        self.chi_0_species = chi_0_species

    def calculate(self, w, ionization_data, beta_rad, t_electron, t_rad,
        beta_electron, levels, delta_input):
        if delta_input is None:
            if self.departure_coefficient is None:
                departure_coefficient = 1. / w
            else:
                departure_coefficient = self.departure_coefficient
            chi_0_species=self.chi_0_species
            chi_0 = ionization_data.ionization_energy.ix[chi_0_species]
            radiation_field_correction = -np.ones((len(ionization_data), len(
                beta_rad)))
            less_than_chi_0 = (
                ionization_data.ionization_energy < chi_0).values
            factor_a = (t_electron / (departure_coefficient * w * t_rad))
            radiation_field_correction[~less_than_chi_0] = factor_a * \
                np.exp(np.outer(ionization_data.ionization_energy.values[
                ~less_than_chi_0], beta_rad - beta_electron))
            radiation_field_correction[less_than_chi_0] = 1 - np.exp(np.outer(
                ionization_data.ionization_energy.values[less_than_chi_0],
                beta_rad) - beta_rad * chi_0)
            radiation_field_correction[less_than_chi_0] += factor_a * np.exp(
                np.outer(ionization_data.ionization_energy.values[
                less_than_chi_0],beta_rad) - chi_0 * beta_electron)
        else:
            radiation_field_correction = np.ones((len(ionization_data),
                len(beta_rad))) * delta_input
        delta = pd.DataFrame(radiation_field_correction,
            columns=np.arange(len(t_rad)), index=ionization_data.index)
        return delta

class IonNumberDensity(ProcessingPlasmaProperty):
    """
    Calculate the ionization balance

    .. math::
        N(X) = N_1 + N_2 + N_3 + \\dots
        N(X) = (N_2/N_1) \\times N_1 + (N3/N2) \\times (N_2/N_1) \\times N_1 + \\dots
        N(X) = N_1(1 + N_2/N_1 + (N_3/N_2) \\times (N_2/N_1) + \\dots
        N(X) = N_1(1+ \\Phi_{i,j}/N_e + \\Phi_{i, j}/N_e \\times \\Phi_{i, j+1}/N_e + \\dots)

    """
    outputs = ('ion_number_density', 'electron_densities')
    latex_name = (r'$N_{i,j}$', r'$n_e$')
    latex_formula = (r'$N(X) = N_1 + N_2 + N_3 + \dots \\ '
                     r'N(X) = (N_2/N_1) \times N_1 + (N3/N2) '
                     r'\times (N_2/N_1) \times N_1 + \dots \\'
                     r'N(X) = N_1(1 + N_2/N_1 + (N_3/N_2) \times (N_2/N_1) '
                     r'+ \dots \\'
                     r'N(X) = N_1(1+ \Phi_{i,j}/N_e + \Phi_{i, j}/N_e '
                     r'\times \Phi_{i, j+1}/N_e + \dots)$')

    outputs = ('ion_number_density', 'electron_densities')

    def __init__(self, plasma_parent, ion_zero_threshold=1e-20):
        super(IonNumberDensity, self).__init__(plasma_parent)
        self.ion_zero_threshold = ion_zero_threshold

    def calculate_with_n_electron(self, phi, partition_function,
                                  number_density, n_electron):
        ion_populations = pd.DataFrame(data=0.0,
            index=partition_function.index.copy(),
            columns=partition_function.columns.copy(), dtype=np.float64)

        for atomic_number, groups in phi.groupby(level='atomic_number'):
            current_phis = (groups / n_electron).replace(np.nan, 0.0).values
            phis_product = np.cumproduct(current_phis, axis=0)
            neutral_atom_density = (number_density.ix[atomic_number] /
                                    (1 + np.sum(phis_product, axis=0)))
            ion_populations.ix[atomic_number, 0] = (
                neutral_atom_density.values)
            ion_populations.ix[atomic_number].values[1:] = (
                neutral_atom_density.values * phis_product)
            ion_populations[ion_populations < self.ion_zero_threshold] = 0.0
        return ion_populations

    def calculate(self, phi, partition_function, number_density):
        n_e_convergence_threshold = 0.05
        n_electron = number_density.sum(axis=0)
        n_electron_iterations = 0
        while True:
            ion_number_density = self.calculate_with_n_electron(
                phi, partition_function, number_density, n_electron)
            ion_numbers = ion_number_density.index.get_level_values(1).values
            ion_numbers = ion_numbers.reshape((ion_numbers.shape[0], 1))
            new_n_electron = (ion_number_density.values * ion_numbers).sum(
                axis=0)
            if np.any(np.isnan(new_n_electron)):
                raise PlasmaIonizationError('n_electron just turned "nan" -'
                                            ' aborting')
            n_electron_iterations += 1
            if n_electron_iterations > 100:
                logger.warn('n_electron iterations above 100 ({0}) -'
                            ' something is probably wrong'.format(
                    n_electron_iterations))
            if np.all(np.abs(new_n_electron - n_electron)
                              / n_electron < n_e_convergence_threshold):
                break
            n_electron = 0.5 * (new_n_electron + n_electron)
        return ion_number_density, n_electron
