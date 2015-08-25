import logging

import numpy as np
import pandas as pd
from scipy import interpolate
from astropy import constants, units
import math
from scipy.integrate import trapz

from tardis.plasma.properties.base import ProcessingPlasmaProperty
from tardis.plasma.exceptions import PlasmaIonizationError

logger = logging.getLogger(__name__)

__all__ = ['PhiSahaNebular', 'PhiSahaLTE', 'PhiSahaNLTE', 'PhiSahaNoNLTE',
           'RadiationFieldCorrection', 'IonNumberDensity']

class PhiSahaLTE(ProcessingPlasmaProperty):
    """
    Outputs:
    general_phi : Pandas DataFrame
        Used as basis for PhiSahaLTE and PhiSahaNebular. Identical output to
        PhiSahaLTE. Separate property required as PhiSahaNebular is based
        on PhiSahaLTE, but the code cannot deal with the inclusion of two
        properties that generate a property called 'phi'.
    """
    outputs = ('general_phi',)
    latex_name = ('\\Phi',)
    latex_formula = ('\\dfrac{2Z_{i,j+1}}{Z_{i,j}}\\Big(\
                     \\dfrac{2\\pi m_{e}/\\beta_{\\textrm{rad}}}{h^2}\
                     \\Big)^{3/2}e^{\\dfrac{-\\chi_{i,j}}{kT_{\
                     \\textrm{rad}}}}',)

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
    outputs = ('general_phi',)
    latex_name = ('\\Phi',)
    latex_formula = ('W(\\delta\\zeta_{i,j}+W(1-\\zeta_{i,j}))\\left(\
                     \\dfrac{T_{\\textrm{electron}}}{T_{\\textrm{rad}}}\
                     \\right)^{1/2}',)

    @staticmethod
    def calculate(t_rad, w, zeta_data, t_electrons, delta,
            g_electron, beta_rad, partition_function, ionization_data):
        phi_lte = PhiSahaLTE.calculate(g_electron, beta_rad,
            partition_function, ionization_data)
        zeta = PhiSahaNebular.get_zeta_values(zeta_data, phi_lte, t_rad)
        phis = phi_lte * w * ((zeta * delta) + w * (1 - zeta)) * \
               (t_electrons/t_rad) ** .5
        return phis

    @staticmethod
    def get_zeta_values(zeta_data, general_phi, t_rad):
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
        return zeta

class PhiSahaNLTE(ProcessingPlasmaProperty):
    outputs = ('phi',)
    def calculate(self, general_phi, g_electron, beta_rad, partition_function,
        ionization_data, nlte_ionization_species, beta_electron, t_rad,
        previous_electron_densities):
        for species in nlte_ionization_species:
            number_of_ions = species+1
            partition_function_species, ionization_data_species = \
                self.filter_phi_input_data(species, partition_function,
                    ionization_data)
            phi_lte = PhiSahaLTE.calculate(g_electron, beta_rad,
                partition_function_species, ionization_data_species)
            r_lu_matrix = np.zeros((species+1, species+1, len(beta_rad)))
            r_lu_matrix_reshaped = r_lu_matrix.reshape(((species+1)**2,
                len(beta_rad)))
            r_ul_matrix = np.zeros((species+1, species+1, len(beta_rad)))
            r_ul_matrix_reshaped = r_ul_matrix.reshape(((species+1)**2,
                len(beta_rad)))
            for ion in range(0, species):
                ionization_energy = ionization_data.ix[species].ix[
                    ion+1].ionization_energy * units.erg
                ionization_nu = ionization_energy / constants.h.cgs
                alpha_st, alpha_sp = self.calculate_rate_upper_lower(phi_lte,
                    ionization_nu, species, ion, beta_electron, t_rad)
                gamma = alpha_sp / phi_lte.ix[species].ix[ion+1]
                lower_index = (number_of_ions + 1) * ion
                upper_index = lower_index + number_of_ions
                r_lu_matrix_reshaped[lower_index] = -gamma
                r_lu_matrix_reshaped[upper_index] = gamma
                r_ul_matrix_reshaped[lower_index+1] = (alpha_st + alpha_sp) * \
                    previous_electron_densities
                r_ul_matrix_reshaped[upper_index+1] = -(alpha_st + alpha_sp) *\
                    previous_electron_densities
            rates_matrix = r_lu_matrix + r_ul_matrix
            rates_matrix[0, :, :] = 1.0
            x = np.zeros(rates_matrix.shape[0])
            x[0] = 1.0
            for i in xrange(len(beta_rad)):
                phi = \
                    np.linalg.solve(rates_matrix[:, :, i], x)
                general_phi[i].ix[species] = \
                    (phi[1:] / phi[:-1])
        print general_phi.ix[species]
        return general_phi

    @staticmethod
    def calculate_rate_upper_lower(phi_lte, ionization_nu, species, ion,
        beta_electron, t_rad):
#Assuming for now that a_{ik} = (v_{i}/v)^3
        sp_coefficient = 8 * np.pi * phi_lte.ix[species].ix[ion+1] * (
            constants.c.cgs.value)**(-2.0) * (ionization_nu)**(3.0)
        x_interval = 9*ionization_nu.value
        x_values = np.arange(ionization_nu.value, 10*ionization_nu.value, (
            x_interval/1000.0))
        exponential_coefficient = -constants.h.cgs.value * beta_electron
        alpha_sp_dataframe = pd.DataFrame(1.0, index=range(len(
            exponential_coefficient)),
            columns=range(len(x_values)))
        y_values_sp = (np.exp(alpha_sp_dataframe.mul(exponential_coefficient,
            axis=0).mul(x_values, axis=1))).mul(1/x_values, axis=1)
        alpha_sp = trapz(y_values_sp, x_values) * sp_coefficient
        alpha_st_dataframe = pd.DataFrame(1.0, index=range(len(
            exponential_coefficient)),
            columns=range(len(x_values)))
        st_coefficient = 4 * np.pi * phi_lte.ix[species].ix[ion+1] * (
            1 / constants.h.cgs.value) * (ionization_nu)**(3.0)
        J_nu_dataframe = pd.DataFrame(1.0, index=range(len(
            exponential_coefficient)), columns=range(len(x_values)))
        J_nu = (J_nu_dataframe * 2.0 * constants.h.cgs.value * \
            constants.c.cgs.value**(-2.0)).mul(x_values**(3.0), axis=1).mul((
            1 / np.exp(exponential_coefficient * t_rad) - 1.0), axis=0)
        y_values_st = (np.exp(alpha_st_dataframe.mul(exponential_coefficient,
            axis=0).mul(x_values, axis=1))).mul(x_values**-4.0, axis=1).mul(
            J_nu, axis=1)
        alpha_st = trapz(y_values_st, x_values) * st_coefficient
        return alpha_st, alpha_sp

    @staticmethod
    def filter_phi_input_data(atomic_number, partition_function,
        ionization_data):
        index_tuples = []
        for ion in range(atomic_number+1):
            index_tuples.append((atomic_number, ion))
        partition_function_index = pd.MultiIndex.from_tuples(index_tuples,
            names=['atomic_number', 'ion_number'])
        ionization_data_index = pd.MultiIndex.from_tuples(index_tuples[1:],
            names=['atomic_number', 'ion_number'])
        partition_function = pd.DataFrame(
            partition_function.ix[atomic_number].values,
            index=partition_function_index, columns=partition_function.columns)
        ionization_data = pd.DataFrame(
            ionization_data.ix[atomic_number].values,
            index=ionization_data_index, columns=['ionization_energy'])
        return partition_function, ionization_data

class PhiSahaNoNLTE(ProcessingPlasmaProperty):
    outputs = ('phi',)
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
    latex_name = ('\\delta',)

    def __init__(self, plasma_parent, departure_coefficient=None):
        super(RadiationFieldCorrection, self).__init__(plasma_parent)
        self.departure_coefficient = departure_coefficient

    def calculate(self, w, ionization_data, beta_rad, t_electrons, t_rad,
        beta_electron, delta_input, chi_0):
        if delta_input is None:
            if self.departure_coefficient is None:
                departure_coefficient = 1. / w
            else:
                departure_coefficient = self.departure_coefficient
            radiation_field_correction = -np.ones((len(ionization_data), len(
                beta_rad)))
            less_than_chi_0 = (
                ionization_data.ionization_energy < chi_0).values
            factor_a = (t_electrons / (departure_coefficient * w * t_rad))
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
    Outputs:
    ion_number_density: Pandas DataFrame
    electron_densities: Numpy Array
        Convergence process to find the correct solution. A trial value for
        the electron density is initiated in a particular zone. The ion
        number densities are then calculated using the Saha equation. The
        electron density is then re-calculated by using the ion number
        densities to sum over the number of free electrons. If the two values
        for the electron densities are not similar to within the threshold
        value, a new guess for the value of the electron density is chosen
        and the process is repeated.
    """
    outputs = ('ion_number_density', 'electron_densities')
    latex_name = ('N_{i,j}','n_{e}',)

    def __init__(self, plasma_parent, ion_zero_threshold=1e-20):
        super(IonNumberDensity, self).__init__(plasma_parent)
        self.ion_zero_threshold = ion_zero_threshold

    def update_helium_nlte(self, ion_number_density, number_density):
        ion_number_density.ix[2].ix[0] = 0.0
        ion_number_density.ix[2].ix[2] = 0.0
        ion_number_density.ix[2].ix[1].update(number_density.ix[2])
        return ion_number_density

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
            if hasattr(self.plasma_parent, 'plasma_properties_dict'):
                if 'HeliumNLTE' in \
                    self.plasma_parent.plasma_properties_dict.keys():
                    ion_number_density = \
                        self.update_helium_nlte(ion_number_density,
                        number_density)
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
