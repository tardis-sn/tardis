import logging
import warnings

import numpy as np
import pandas as pd
import numexpr as ne

from scipy import interpolate

from tardis.plasma.properties.base import ProcessingPlasmaProperty
from tardis.plasma.exceptions import PlasmaIonizationError


logger = logging.getLogger(__name__)

__all__ = ['PhiSahaNebular', 'PhiSahaLTE', 'RadiationFieldCorrection',
           'IonNumberDensity', 'IonNumberDensityHeNLTE']


def calculate_block_ids_from_dataframe(dataframe):
        block_start_id = np.where(np.diff(
            dataframe.index.get_level_values(0)) != 0.0)[0] + 1
        return np.hstack(([0], block_start_id, [len(dataframe)]))


class PhiSahaLTE(ProcessingPlasmaProperty):
    """
    Attributes:
    phi : Pandas DataFrame, dtype float
          Used for LTE ionization. Indexed by atomic number, ion number. Columns are zones.
    """
    outputs = ('phi',)
    latex_name = ('\\Phi',)
    latex_formula = ('\\dfrac{2Z_{i,j+1}}{Z_{i,j}}\\Big(\
                     \\dfrac{2\\pi m_{e}/\\beta_{\\textrm{rad}}}{h^2}\
                     \\Big)^{3/2}e^{\\dfrac{-\\chi_{i,j}}{kT_{\
                     \\textrm{rad}}}}',)

    broadcast_ionization_energy = None

    @staticmethod
    def calculate(g_electron, beta_rad, partition_function, ionization_data):

        phis = np.empty(
            (partition_function.shape[0] -
             partition_function.index.get_level_values(0).unique().size,
            partition_function.shape[1]))

        block_ids = calculate_block_ids_from_dataframe(partition_function)

        for i, start_id in enumerate(block_ids[:-1]):
            end_id = block_ids[i + 1]
            current_block = partition_function.values[start_id:end_id]
            current_phis = current_block[1:] / current_block[:-1]
            phis[start_id - i:end_id - i - 1] = current_phis

        broadcast_ionization_energy = (
            ionization_data.loc[
                partition_function.index].dropna())
        phi_index = broadcast_ionization_energy.index
        broadcast_ionization_energy = broadcast_ionization_energy.values

        phi_coefficient = (2 * g_electron * np.exp(
            np.outer(broadcast_ionization_energy, -beta_rad)))

        return pd.DataFrame(phis * phi_coefficient, index=phi_index)

    @staticmethod
    def _calculate_block_ids(partition_function):
        partition_function.index.get_level_values(0).unique()



class PhiSahaNebular(ProcessingPlasmaProperty):
    """
    Attributes:
    phi : Pandas DataFrame, dtype float
          Used for nebular ionization. Indexed by atomic number, ion number. Columns are zones.
    """
    outputs = ('phi',)
    latex_name = ('\\Phi',)
    latex_formula = ('W(\\delta\\zeta_{i,j}+W(1-\\zeta_{i,j}))\\left(\
                     \\dfrac{T_{\\textrm{electron}}}{T_{\\textrm{rad}}}\
                     \\right)^{1/2}',)
    @staticmethod
    def calculate(t_rad, w, zeta_data, t_electrons, delta,
            g_electron, beta_rad, partition_function, ionization_data):
        phi_lte = PhiSahaLTE.calculate(g_electron, beta_rad,
            partition_function, ionization_data)
        zeta = PhiSahaNebular.get_zeta_values(zeta_data, phi_lte.index, t_rad)
        phis = phi_lte * w * ((zeta * delta) + w * (1 - zeta)) * \
               (t_electrons/t_rad) ** .5
        return phis

    @staticmethod
    def get_zeta_values(zeta_data, ion_index, t_rad):
        zeta_t_rad = zeta_data.columns.values.astype(np.float64)
        zeta_values = zeta_data.loc[ion_index].values.astype(np.float64)
        zeta = interpolate.interp1d(zeta_t_rad, zeta_values, bounds_error=False,
                                    fill_value=np.nan)(t_rad)
        zeta = zeta.astype(float)

        if np.any(np.isnan(zeta)):
            warnings.warn('t_rads outside of zeta factor interpolation'
                                 ' zeta_min={0:.2f} zeta_max={1:.2f} '
                                 '- replacing with 1s'.format(
                zeta_data.columns.values.min(), zeta_data.columns.values.max(),
                t_rad))
            zeta[np.isnan(zeta)] = 1.0

        return zeta

class RadiationFieldCorrection(ProcessingPlasmaProperty):
    """
    Attributes:
    delta : Pandas DataFrame, dtype float
            Calculates the radiation field correction (see Mazzali & Lucy, 1993) if
            not given as input in the config. file. The default chi_0_species is
            Ca II, which is good for type Ia supernovae. For type II supernovae,
            (1, 1) should be used. Indexed by atomic number, ion number. The columns are zones.
    """
    outputs = ('delta',)
    latex_name = ('\\delta',)

    def __init__(self, plasma_parent=None, departure_coefficient=None,
        chi_0_species=(20,2), delta_treatment=None):
        super(RadiationFieldCorrection, self).__init__(plasma_parent)
        self.departure_coefficient = departure_coefficient
        self.delta_treatment = delta_treatment
        self.chi_0_species = chi_0_species

    def _set_chi_0(self, ionization_data):
        if self.chi_0_species == (20, 2):
            self.chi_0 = 1.9020591570241798e-11
        else:
            self.chi_0 = ionization_data.loc[self.chi_0_species]

    def calculate(self, w, ionization_data, beta_rad, t_electrons, t_rad,
        beta_electron):
        if getattr(self, 'chi_0', None) is None:
            self._set_chi_0(ionization_data)
        if self.delta_treatment is None:
            if self.departure_coefficient is None:
                departure_coefficient = 1. / w
            else:
                departure_coefficient = self.departure_coefficient
            radiation_field_correction = -np.ones((len(ionization_data), len(
                beta_rad)))
            less_than_chi_0 = (
                ionization_data < self.chi_0).values
            factor_a = (t_electrons / (departure_coefficient * w * t_rad))
            radiation_field_correction[~less_than_chi_0] = factor_a * \
                np.exp(np.outer(ionization_data.values[
                ~less_than_chi_0], beta_rad - beta_electron))
            radiation_field_correction[less_than_chi_0] = 1 - np.exp(np.outer(
                ionization_data.values[less_than_chi_0],
                beta_rad) - beta_rad * self.chi_0)
            radiation_field_correction[less_than_chi_0] += factor_a * np.exp(
                np.outer(ionization_data.values[
                less_than_chi_0],beta_rad) - self.chi_0 * beta_electron)
        else:
            radiation_field_correction = np.ones((len(ionization_data),
                len(beta_rad))) * self.delta_treatment
        delta = pd.DataFrame(radiation_field_correction,
            columns=np.arange(len(t_rad)), index=ionization_data.index)
        return delta

class IonNumberDensity(ProcessingPlasmaProperty):
    """
    Attributes:
    ion_number_density : Pandas DataFrame, dtype float
                         Index atom number, ion number. Columns zones.
    electron_densities : Numpy Array, dtype float

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

    def __init__(self, plasma_parent, ion_zero_threshold=1e-20, electron_densities=None):
        super(IonNumberDensity, self).__init__(plasma_parent)
        self.ion_zero_threshold = ion_zero_threshold
        self.block_ids = None
        self._electron_densities = electron_densities

    @staticmethod
    def calculate_with_n_electron(phi, partition_function,
                                  number_density, n_electron, block_ids,
                                  ion_zero_threshold):
        if block_ids is None:
            block_ids = IonNumberDensity._calculate_block_ids(phi)

        ion_populations = np.empty_like(partition_function.values)

        phi_electron = np.nan_to_num(phi.values / n_electron.values)

        for i, start_id in enumerate(block_ids[:-1]):
            end_id = block_ids[i + 1]
            current_phis = phi_electron[start_id:end_id]
            phis_product = np.cumprod(current_phis, 0)

            tmp_ion_populations = np.empty((current_phis.shape[0] + 1,
                                            current_phis.shape[1]))
            tmp_ion_populations[0] = (number_density.values[i] /
                                    (1 + np.sum(phis_product, axis=0)))
            tmp_ion_populations[1:] = tmp_ion_populations[0] * phis_product

            ion_populations[start_id + i:end_id + 1 + i] = tmp_ion_populations

        ion_populations[ion_populations < ion_zero_threshold] = 0.0

        return pd.DataFrame(data = ion_populations,
                            index=partition_function.index), block_ids

    @staticmethod
    def _calculate_block_ids(phi):
        return calculate_block_ids_from_dataframe(phi)

    def calculate(self, phi, partition_function, number_density):
        if self._electron_densities is None:
            n_e_convergence_threshold = 0.05
            n_electron = number_density.sum(axis=0)
            n_electron_iterations = 0

            while True:
                ion_number_density, self.block_ids = \
                    self.calculate_with_n_electron(
                    phi, partition_function, number_density, n_electron,
                        self.block_ids, self.ion_zero_threshold)
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
        else:
            n_electron = self._electron_densities
            ion_number_density, self.block_ids = \
                self.calculate_with_n_electron(
                    phi, partition_function, number_density, n_electron,
                    self.block_ids, self.ion_zero_threshold)

        return ion_number_density, n_electron

class IonNumberDensityHeNLTE(ProcessingPlasmaProperty):
    """
    Attributes:
    ion_number_density : Pandas DataFrame, dtype float
                         Index atom number, ion number. Columns zones.
    electron_densities : Numpy Array, dtype float

    Convergence process to find the correct solution. A trial value for
    the electron density is initiated in a particular zone. The ion
    number densities are then calculated using the Saha equation. The
    electron density is then re-calculated by using the ion number
    densities to sum over the number of free electrons. If the two values
    for the electron densities are not similar to within the threshold
    value, a new guess for the value of the electron density is chosen
    and the process is repeated.
    """
    outputs = ('ion_number_density', 'electron_densities',
               'helium_population_updated')
    latex_name = ('N_{i,j}','n_{e}',)

    def __init__(self, plasma_parent, ion_zero_threshold=1e-20, electron_densities=None):
        super(IonNumberDensityHeNLTE, self).__init__(plasma_parent)
        self.ion_zero_threshold = ion_zero_threshold
        self.block_ids = None
        self._electron_densities = electron_densities

    def update_he_population(self, helium_population, n_electron,
                             number_density):
        helium_population_updated = helium_population.copy()
        he_one_population = helium_population_updated.loc[0].mul(n_electron)
        he_three_population = helium_population_updated.loc[2].mul(
            1./n_electron)
        helium_population_updated.loc[0].update(he_one_population)
        helium_population_updated.loc[2].update(he_three_population)
        unnormalised = helium_population_updated.sum()
        normalised = helium_population_updated.mul(number_density.loc[2] /
            unnormalised)
        helium_population_updated.update(normalised)
        return helium_population_updated

    def calculate(self, phi, partition_function, number_density,
                  helium_population):
        if self._electron_densities is None:
            n_e_convergence_threshold = 0.05
            n_electron = number_density.sum(axis=0)
            n_electron_iterations = 0
            while True:
                ion_number_density, self.block_ids = \
                    IonNumberDensity.calculate_with_n_electron(
                    phi, partition_function, number_density, n_electron,
                        self.block_ids, self.ion_zero_threshold)
                helium_population_updated = self.update_he_population(
                    helium_population, n_electron, number_density)
                ion_number_density.loc[2, 0].update(helium_population_updated.loc[
                                                        0].sum(axis=0))
                ion_number_density.loc[2, 1].update(helium_population_updated.loc[
                                                        1].sum(axis=0))
                ion_number_density.loc[2, 2].update(helium_population_updated.loc[
                                                        2, 0])
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
        else:
            n_electron = self._electron_densities
            ion_number_density, self.block_ids = \
                IonNumberDensity.calculate_with_n_electron(
                    phi, partition_function, number_density, n_electron,
                    self.block_ids, self.ion_zero_threshold)

            helium_population_updated = self.update_he_population(
                helium_population, n_electron, number_density)
            ion_number_density.loc[2, 0].update(helium_population_updated.loc[
                                                      0].sum(axis=0))
            ion_number_density.loc[2, 1].update(helium_population_updated.loc[
                                                      1].sum(axis=0))
            ion_number_density.loc[2, 2].update(helium_population_updated.loc[
                                                      2, 0])
        return ion_number_density, n_electron, helium_population_updated
