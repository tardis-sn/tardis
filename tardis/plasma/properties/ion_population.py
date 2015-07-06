import logging

import numpy as np
import pandas as pd
from scipy import interpolate

from tardis.plasma.properties.base import ProcessingPlasmaProperty
from tardis.plasma.exceptions import PlasmaIonizationError

logger = logging.getLogger(__name__)

__all__ = ['PhiSahaNebular', 'PhiSahaLTE', 'RadiationFieldCorrection',
           'IonNumberDensity', 'PhiGeneral', 'IonNumberDensityNLTE']

class PhiGeneral(ProcessingPlasmaProperty):

    outputs = ('general_phi',)

    latex_formula = (r'$\Phi_{i,j} = \frac{N_{i, j+1} n_e}{N_{i, j}} \\'
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
    Calculating the ionization equilibrium using the Saha equation, where i is atomic number,
    j is the ion_number, :math:`n_e` is the electron density, :math:`Z_{i, j}` are the partition functions
    and :math:`\chi` is the ionization energy. For the `NebularPlasma` we first calculate the
    ionization balance assuming LTE conditions (:math:`\\Phi_{i, j}(\\textrm{LTE})`) and use factors to more accurately
    describe the plasma. The two important factors are :math:`\\zeta` - a correction factor to take into account
    ionizations from excited states. The second factor is :math:`\\delta` , adjusting the ionization balance for the fact that
    there's more line blanketing in the blue.

    The :math:`\\zeta` factor for different temperatures is read in to the `~tardis.atomic.NebularAtomData` and then
    interpolated for the current temperature.

    The :math:`\\delta` factor is calculated with :meth:`calculate_radiation_field_correction`.

    Finally the ionization balance is adjusted (as equation 14 in :cite:`1993A&A...279..447M`):

    .. math::


        \\Phi_{i,j} =& \\frac{N_{i, j+1} n_e}{N_{i, j}} \\\\

        \\Phi_{i, j} =& W \\times[\\delta \\zeta + W ( 1 - \\zeta)] \\left(\\frac{T_\\textrm{e}}{T_\\textrm{R}}\\right)^{1/2}
        \\Phi_{i, j}(\\textrm{LTE})

    """

    outputs = ('phi',)

    @staticmethod
    def calculate(general_phi, t_rad, w, zeta_data, t_electron, delta):
        logger.debug('Calculating Saha using Nebular approximation')

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

    outputs = ('phi',)

    @staticmethod
    def calculate(general_phi):
        return general_phi

class RadiationFieldCorrection(ProcessingPlasmaProperty):

    outputs = ('delta',)
    """
    Calculating radiation field correction factors according to Mazzali & Lucy 1993 (:cite:`1993A&A...279..447M`; henceforth ML93)


    In ML93 the radiation field correction factor is denoted as :math:`\\delta` and is calculated in Formula 15 & 20

    The radiation correction factor changes according to a ionization energy threshold :math:`\\chi_\\textrm{T}`
    and the species ionization threshold (from the ground state) :math:`\\chi_0`.

    For :math:`\\chi_\\textrm{T} \\ge \\chi_0`

    .. math::
        \\delta = \\frac{T_\\textrm{e}}{b_1 W T_\\textrm{R}} \\exp(\\frac{\\chi_\\textrm{T}}{k T_\\textrm{R}} -
        \\frac{\\chi_0}{k T_\\textrm{e}})

    For :math:`\\chi_\\textrm{T} < \\chi_0`

    .. math::self.beta_rads * chi_0
        \\delta = 1 - \\exp(\\frac{\\chi_\\textrm{T}}{k T_\\textrm{R}} - \\frac{\\chi_0}{k T_\\textrm{R}}) + \\frac{T_\\textrm{e}}{b_1 W T_\\textrm{R}} \\exp(\\frac{\\chi_\\textrm{T}}{k T_\\textrm{R}} -
        \\frac{\\chi_0}{k T_\\textrm{e}}),

    where :math:`T_\\textrm{R}` is the radiation field Temperature, :math:`T_\\textrm{e}` is the electron temperature and W is the
    dilution factor.

    Parameters
        if line_interaction_type in ('downbranch', 'macroat
    ----------
    phi_table : `~astropy.table.Table`
        a table containing the field 'atomic_number', 'ion_number', 'phi'

    departure_coefficient : `~float` or `~None`, optional
        departure coefficient (:math:`b_1` in ML93) For the default (`None`) it is set to 1/W.

    chi_0_species : `~tuple`, optional
        This describes which ionization energy to use for the threshold. Default is Calcium II
        (1044 Angstrom; useful for Type Ia)
        For Type II supernovae use Lyman break (912 Angstrom) or (1,1) as the tuple

    Returns
    -------

    This function adds a field 'delta' to the phi table given to the function

    """

    def __init__(self, plasma_parent, departure_coefficient=None,
        chi_0_species=(20,2)):
        super(RadiationFieldCorrection, self).__init__(plasma_parent)
        self.departure_coefficient = departure_coefficient
        self.chi_0_species = chi_0_species

    def calculate(self, w, ionization_data, beta_rad, t_electron, t_rad,
        beta_electron, levels, delta_input):
        # factor delta ML 1993
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

class IonNumberDensityNLTE(ProcessingPlasmaProperty):
    """
    Calculate the ionization balance

    .. math::
        N(X) = N_1 + N_2 + N_3 + \\dots
        N(X) = (N_2/N_1) \\times N_1 + (N3/N2) \\times (N_2/N_1) \\times N_1 + \\dots
        N(X) = N_1(1 + N_2/N_1 + (N_3/N_2) \\times (N_2/N_1) + \\dots
        N(X) = N_1(1+ \\Phi_{i,j}/N_e + \\Phi_{i, j}/N_e \\times \\Phi_{i, j+1}/N_e + \\dots)

    """

    latex_formula = (r'$N(X) = N_1 + N_2 + N_3 + \dots \\ '
                     r'N(X) = (N_2/N_1) \times N_1 + (N3/N2) '
                     r'\times (N_2/N_1) \times N_1 + \dots \\'
                     r'N(X) = N_1(1 + N_2/N_1 + (N_3/N_2) \times (N_2/N_1) '
                     r'+ \dots \\'
                     r'N(X) = N_1(1+ \Phi_{i,j}/N_e + \Phi_{i, j}/N_e '
                     r'\times \Phi_{i, j+1}/N_e + \dots)$')

    outputs = ('ion_number_density', 'electron_densities',
               'level_boltzmann_factor')

    def __init__(self, plasma_parent, ion_zero_threshold=1e-20):
        super(IonNumberDensityNLTE, self).__init__(plasma_parent)
        self.ion_zero_threshold = ion_zero_threshold

    def calculate_nlte(self, nlte_input, beta_sobolev, j_blues, levels,
        nlte_data, t_electron, electron_densities,
        general_level_boltzmann_factor):
        """
        Calculating the NLTE level populations for specific ions

        """
        if j_blues.empty:
            j_blues = np.zeros_like(nlte_data.lines.index)
        if nlte_input.get('coronal_approximation', False):
            beta_sobolev_nlte = np.ones_like(beta_sobolev)
            j_blues_nlte = np.zeros_like(j_blues)
            logger.info('using coronal approximation = setting beta_sobolevs to 1 AND j_blues to 0')
        else:
            beta_sobolev_nlte = beta_sobolev
            j_blues_nlte = j_blues

        if nlte_input.get('classical_nebular', False):
            logger.info('using Classical Nebular = setting beta_sobolevs to 1')
            beta_sobolev_nlte = np.ones_like(beta_sobolev)

        for species in nlte_input.species:
            logger.info('Calculating rates for species %s', species)
            number_of_levels = len(levels.ix[species])

            lnl = nlte_data.lines_level_number_lower[species]
            lnu = nlte_data.lines_level_number_upper[species]

            lines_index = nlte_data.lines_idx[species]
            A_uls = nlte_data.A_uls[species]
            B_uls = nlte_data.B_uls[species]
            B_lus = nlte_data.B_lus[species]

            r_lu_index = lnu * number_of_levels + lnl
            r_ul_index = lnl * number_of_levels + lnu

            j_blues_nlte = j_blues_nlte[lines_index]

            r_ul_matrix = np.zeros((number_of_levels, number_of_levels,
                len(t_electron)), dtype=np.float64)
            r_ul_matrix_reshaped = r_ul_matrix.reshape((number_of_levels**2,
                len(t_electron)))
            r_ul_matrix_reshaped[r_ul_index] = A_uls[np.newaxis].T + \
                B_uls[np.newaxis].T * j_blues_nlte[np.newaxis].T
            r_ul_matrix_reshaped[r_ul_index] *= beta_sobolev_nlte[lines_index]

            r_lu_matrix = np.zeros_like(r_ul_matrix)
            r_lu_matrix_reshaped = r_lu_matrix.reshape((number_of_levels**2,
                len(t_electron)))
            r_lu_matrix_reshaped[r_lu_index] = B_lus[np.newaxis].T * \
                j_blues_nlte[np.newaxis].T * \
                beta_sobolev_nlte[lines_index]

#            collision_matrix = nlte_data.get_collision_matrix(species,
#                t_electron) * electron_densities.values

            collision_matrix = r_ul_matrix.copy()
            collision_matrix.fill(0.0)

            rates_matrix = r_lu_matrix + r_ul_matrix + collision_matrix

            level_boltzmann_factor = general_level_boltzmann_factor

            for i in xrange(number_of_levels):
                rates_matrix[i, i] = -rates_matrix[:, i].sum(axis=0)

            rates_matrix[0, :, :] = 1.0

            x = np.zeros(rates_matrix.shape[0])
            x[0] = 1.0
            for i in xrange(len(t_electron)):
                nlte_level_boltzmann_factor = \
                    np.linalg.solve(rates_matrix[:, :, i], x)
                level_boltzmann_factor[i].ix[species] = \
                    nlte_level_boltzmann_factor
            return level_boltzmann_factor

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

    def update_dependent_properties(self, levels, level_boltzmann_factor,
        g_electron, beta_rad, ionization_data, w, t_rad, zeta_data, t_electron,
        delta):
        partition_function = self.plasma_parent.plasma_properties_dict[
            'PartitionFunction'].calculate(levels, level_boltzmann_factor)
        general_phi = PhiGeneral.calculate(g_electron, beta_rad,
            partition_function, ionization_data)
        if self.plasma_parent.ionization_mode == 'lte':
            phi = PhiSahaLTE.calculate(general_phi)
        elif self.plasma_parent.excitation_mode == 'dilute-lte':
            phi = PhiSahaNebular.calculate(general_phi, t_rad, w, zeta_data,
                t_electron, delta)
        level_population_fraction = LevelPopulation.calculate(levels,
            partition_function, level_boltzmann_factor)

        return partition_function, phi, level_population_fraction

    def calculate_level_boltzmann_factor(self, nlte_input, j_blues, levels,
        t_electron, n_electron, general_level_boltzmann_factor, nlte_data,
        g_electron, beta_rad, ionization_data, w, t_rad, zeta_data, delta):
        beta_sobolev = np.ones((len(nlte_data.lines), len(t_electron)))
        while True:
            level_boltzmann_factor = self.calculate_nlte(
                nlte_input, beta_sobolev, j_blues, levels, nlte_data,
                t_electron, n_electron, general_level_boltzmann_factor)
            partition_function, phi, level_population_fraction = \
                self.update_dependent_properties(levels,
                level_boltzmann_factor, g_electron, beta_rad, ionization_data,
                w, t_rad, zeta_data, delta)

            if np.all(np.abs(new_level_boltzmann_factor -
                level_boltzmann_factor) / n_electron < 0.05):
                break
            return

    def calculate(self, number_density, nlte_input,
        j_blues, levels, nlte_data, t_electron,
        g_electron, beta_rad,
        ionization_data, w, t_rad, zeta_data, delta, lines_lower_level_index,
        lines_upper_level_index, lines, time_explosion,
        general_level_boltzmann_factor, atomic_data):
        n_e_convergence_threshold = 0.05
        n_electron = number_density.sum(axis=0)
        n_electron_iterations = 0
        while True:
            self.calculate_level_boltzmann_factor(nlte_input, j_blues,
                levels, t_electron, n_electron, general_level_boltzmann_factor,
                nlte_data, g_electron, beta_rad, ionization_data, w, t_rad,
                zeta_data, delta)
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
        return (ion_number_density, n_electron, level_boltzmann_factor,
            partition_function, phi, level_population_fraction,
            level_number_density, stimulated_emission_factor, tau_sobolevs,
            beta_sobolev)
